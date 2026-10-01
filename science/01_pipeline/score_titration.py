#!/usr/bin/env python3
"""Score the H2 MSA-depth titration against the pre-registered criterion.

    python3 score_titration.py --work ~/h2 --kit ~/apheris_kit_rorgamma

Applies BOTH pre-registered conditions:

  1. fold-quality floor : mean pLDDT >= 70
  2. state call         : region RMSD(479-486) to 4ZJW < to 3KYT, margin >= 1.0 A

and the primary criterion: >= 2 of 5 seeds satisfying both at some depth.

This also fixes the gap found in H1 scoring, where structures 27 A from BOTH
references were reported as decisive state calls because the margin test is
purely relative. Here a structure failing the floor is `unscorable` — it cannot
count as a success, and it is printed rather than silently dropped.

Writes titration_results.csv and titration_summary.md next to --work.
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import re
import sys

import numpy as np

REGION = list(range(479, 487))
PLDDT_FLOOR = 70.0
MARGIN = 1.0
MIN_SEEDS = 2


def mean_plddt(cif: pathlib.Path):
    """Mean pLDDT from the sibling scores json, else from the CIF B-factors."""
    for cand in (
        cif.with_name(cif.name.replace("_model.cif", "_scores.json")),
        cif.with_name(cif.stem + "_scores.json"),
        cif.parent / "confidence.json",
    ):
        if cand.exists():
            try:
                s = json.load(open(cand))
            except Exception:
                continue
            for k, v in s.items():
                if "local_distance" in k and isinstance(v, list) and v:
                    m = float(np.mean(v))
                    return m * 100 if m <= 1.0 else m
            for k in ("mean_plddt", "plddt", "complex_plddt"):
                if k in s and isinstance(s[k], (int, float)):
                    m = float(s[k])
                    return m * 100 if m <= 1.0 else m
    return None


def plddt_from_bfactor(chain):
    vals = []
    for res in chain.residues:
        ca = res.find_atom("CA", "*")
        if ca is not None:
            vals.append(ca.b_iso)
    if not vals:
        return None
    m = float(np.mean(vals))
    return m * 100 if m <= 1.0 else m


def parse_tags(path: pathlib.Path):
    s = str(path)
    depth = None
    m = re.search(r"depth_(\d+)", s)
    if m:
        depth = int(m.group(1))
    seed = None
    m = re.search(r"seed[_-](\d+)", s)
    if m:
        seed = int(m.group(1))
    samp = None
    m = re.search(r"sample[_-](\d+)", s)
    if m:
        samp = int(m.group(1))
    return depth, seed, samp


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", default="~/h2")
    ap.add_argument("--kit", required=True,
                    help="path to apheris_kit_rorgamma (for refs_rorgamma/)")
    ap.add_argument("--rgkit", default=None,
                    help="path to the rgkit package dir, if not inside --kit")
    args = ap.parse_args()

    work = pathlib.Path(args.work).expanduser()
    kit = pathlib.Path(args.kit).expanduser()

    for p in filter(None, [args.rgkit, str(kit), str(kit.parent)]):
        sys.path.insert(0, p)
    try:
        from rgkit import state_recovery2 as SR
        from rgkit import structio as S
    except ImportError:
        print("ERROR: cannot import rgkit. Pass --rgkit /path/to/dir containing "
              "the rgkit package (from the cross-check bundle).", file=sys.stderr)
        return 2

    refs_dir = kit / "refs_rorgamma"
    a = S.load_chain(refs_dir / "3KYT_clean.cif", source="3KYT")   # H12-in
    b = S.load_chain(refs_dir / "4ZJW_clean.cif", source="4ZJW")   # H12-out
    refs = {"active": a, "inactive": b}

    cifs = sorted(p for p in (work / "out").rglob("*.cif") if "__MACOSX" not in str(p))
    if not cifs:
        print(f"no .cif under {work/'out'}", file=sys.stderr)
        return 2

    rows = []
    for p in cifs:
        depth, seed, samp = parse_tags(p)
        row = {"depth": depth, "seed": seed, "sample": samp,
               "file": str(p.relative_to(work)), "status": "", "plddt": "",
               "rmsd_H12in": "", "rmsd_H12out": "", "margin": "",
               "verdict": "", "counts_as_success": 0}
        try:
            c = S.load_chain(p, source=p.name)
        except Exception as exc:
            row["status"] = f"load-error: {exc}"[:80]
            rows.append(row)
            continue

        pl = mean_plddt(p)
        if pl is None:
            pl = plddt_from_bfactor(c)
        row["plddt"] = round(pl, 1) if pl is not None else ""

        r = SR.state_recovery_generic(c, refs, REGION, decisive_margin=MARGIN)
        d_in, d_out = r["region_rmsd_active"], r["region_rmsd_inactive"]
        row["rmsd_H12in"] = round(d_in, 3)
        row["rmsd_H12out"] = round(d_out, 3)
        row["margin"] = round(abs(d_in - d_out), 3)

        if pl is not None and pl < PLDDT_FLOOR:
            row["status"] = "unscorable"
            row["verdict"] = f"fold failed (pLDDT {pl:.1f} < {PLDDT_FLOOR:g})"
        else:
            row["status"] = "ok"
            if d_out < d_in and (d_in - d_out) >= MARGIN:
                row["verdict"] = "H12-out"
                row["counts_as_success"] = 1
            elif d_in < d_out and (d_out - d_in) >= MARGIN:
                row["verdict"] = "H12-in"
            else:
                row["verdict"] = "non-decisive"
        rows.append(row)

    out_csv = work / "titration_results.csv"
    with out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # ---- per-depth summary against the pre-registered criterion -------------
    depths = sorted({r["depth"] for r in rows if r["depth"] is not None})
    lines = ["# H2 titration — result against the pre-registered criterion", "",
             f"Criterion: at some depth, >= {MIN_SEEDS}/5 seeds give mean pLDDT "
             f">= {PLDDT_FLOOR:g} AND H12-out with margin >= {MARGIN:g} A.", "",
             "| depth | structures | folded (pLDDT>=70) | mean pLDDT | H12-out structures | seeds with >=1 H12-out | criterion |",
             "|---|---|---|---|---|---|---|"]
    met_any = False
    print(f"\n{'depth':>8} {'n':>5} {'folded':>7} {'meanpLDDT':>10} "
          f"{'H12-out':>8} {'seeds':>6}  criterion")
    for d in depths:
        sub = [r for r in rows if r["depth"] == d]
        folded = [r for r in sub if r["status"] == "ok"]
        succ = [r for r in sub if r["counts_as_success"]]
        seeds_ok = sorted({r["seed"] for r in succ if r["seed"] is not None})
        pls = [r["plddt"] for r in sub if isinstance(r["plddt"], float)]
        mp = f"{np.mean(pls):.1f}" if pls else "-"
        met = len(seeds_ok) >= MIN_SEEDS
        met_any |= met
        verdict = "MET" if met else "not met"
        print(f"{d:>8} {len(sub):>5} {len(folded):>7} {mp:>10} "
              f"{len(succ):>8} {len(seeds_ok):>6}  {verdict}")
        lines.append(f"| {d} | {len(sub)} | {len(folded)} | {mp} | {len(succ)} | "
                     f"{len(seeds_ok)} ({seeds_ok}) | **{verdict}** |")

    lines += ["", "## Outcome", ""]
    if met_any:
        lines.append("**H2 supported.** At least one MSA depth produced H12-out in "
                     f">= {MIN_SEEDS} independent seeds with the fold intact. Report the "
                     "depth window and the fold-quality curve beside it.")
    else:
        lines.append("**H2 not supported — declared negative applies.** On RORγ LBD, "
                     "fold stability and H12-state accessibility could not be separated "
                     "by MSA depth: the state was inaccessible at every depth that "
                     "preserved the fold. This closes the mechanism the literature "
                     "suggests should work, and is the stronger of the two negatives.")

    n_uns = sum(1 for r in rows if r["status"] == "unscorable")
    if n_uns:
        lines += ["", f"{n_uns} structure(s) failed the pLDDT floor and are recorded "
                      "as `unscorable` in titration_results.csv — listed, not dropped."]

    (work / "titration_summary.md").write_text("\n".join(lines) + "\n")
    print(f"\n-> {out_csv}\n-> {work/'titration_summary.md'}")
    print("\nOUTCOME:", "H2 SUPPORTED" if met_any else "H2 NOT SUPPORTED (declared negative)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
