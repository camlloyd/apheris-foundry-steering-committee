"""
Score a whole experiment grid and produce the table you put on the slide.

Two modes:

  1. Validate your references FIRST. Run this before anything else:
         python score_run.py validate --config targets.json --target abl1
     It classifies each reference structure and checks it comes back as the
     state you think it is. If your DFG-out reference does not classify as
     DFG-out, every number you produce afterwards is meaningless. This takes
     thirty seconds and has saved entire hackathon days.

  2. Score predictions:
         python score_run.py score --config targets.json --target abl1 \
             --runs runs/abl1 --out results.csv

     Walks each arm's output directory, scores every sample, and reports
     state recovery and pose accuracy as SEPARATE columns -- because
     KinConfBench showed they decouple, and conflating them is the mistake
     the literature has already made.

Headline metric: state recovery rate = fraction of samples per arm that land
in the target conformational state. Not mean RMSD, which hides a bimodal
distribution: an arm that gets it right half the time and badly wrong half the
time looks identical to an arm that is consistently mediocre.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
from collections import defaultdict

from structure_io import load_structure
from kinase_state import (
    classify_dfg, find_kinase_motifs, motifs_from_overrides, transfer_motifs,
)
from state_recovery import (
    contact_recovery, ligand_centroid_distance, state_recovery,
)


def _find_prediction_files(arm_dir: str) -> list:
    """Model outputs, whatever the tool decided to call them."""
    patterns = ["*.cif", "*.pdb", "*.mmcif",
                "**/*.cif", "**/*.pdb", "**/*.mmcif"]
    found = []
    for p in patterns:
        found.extend(glob.glob(os.path.join(arm_dir, p), recursive=True))
    return sorted(set(found))


# ------------------------------------------------------------- mode: validate

def validate_references_generic(cfg: dict) -> int:
    """
    Non-kinase validation ("order_parameter": "none"): the question is not
    DFG geometry but whether the two references actually define a state
    difference at all. derive_state_region answers that from the data.
    """
    from state_recovery2 import derive_state_region

    print("Validating reference structures (generic two-state check)\n"
          + "-" * 60)
    problems = 0
    paths = []
    for role, key in (("target state", "ref_target_state"),
                      ("other state", "ref_other_state")):
        path = cfg.get(key)
        paths.append(path)
        if not path or not os.path.exists(path):
            print(f"[MISSING] {role:<14s} {path}  -- download it first")
            problems += 1
            continue
        st = load_structure(path)
        ligs = st.ligands()
        print(f"[  ok ] {role:<14s} {os.path.basename(path)}")
        lig_note = "NONE -- pocket derivation needs a ligand here"
        print(f"          ligands: {[l.name for l in ligs] or lig_note}")
        if not ligs:
            problems += 1

    if problems:
        print(f"\n{problems} problem(s). Fix these before scoring anything.")
        return problems

    a = load_structure(paths[0])
    b = load_structure(paths[1])
    der = derive_state_region(a, b)
    maxd = max((d for _, d in der["profile"]), default=0.0)
    print(f"\n  two-state check: {len(der['region'])} region residues, "
          f"segments {der['segments']}, max CA displacement {maxd:.2f} A")
    for w in der["warnings"]:
        print(f"  ! {w}")
    if not der["region"] or maxd < 2.0:
        problems += 1
        print("  -> the two references do not define a usable state "
              "difference. Pick a different pair before scoring.")
    else:
        print("  references look good -- go ahead and score")
    return problems


def validate_references(cfg: dict) -> int:
    if cfg.get("order_parameter", "dfg") == "none":
        return validate_references_generic(cfg)

    print("Validating reference structures\n" + "-" * 60)
    problems = 0

    for role, key in (("target state", "ref_target_state"),
                      ("other state", "ref_other_state")):
        path = cfg.get(key)
        expected = cfg["target_state"] if key == "ref_target_state" else cfg["other_state"]

        if not path or not os.path.exists(path):
            print(f"[MISSING] {role:<14s} {path}  -- download it first")
            problems += 1
            continue

        st = load_structure(path)
        overrides = cfg.get("motif_overrides")
        motifs = (motifs_from_overrides(st, overrides) if overrides
                  else find_kinase_motifs(st))
        res = classify_dfg(st, motifs)

        ok = res.label == expected
        tag = "  ok  " if ok else " WRONG"
        print(f"[{tag}] {role:<14s} {os.path.basename(path)}")
        print(f"          expected {expected}, measured {res}")
        print(f"          motifs: {motifs.describe(st)}")

        ligs = st.ligands()
        if ligs:
            print(f"          ligands: {[l.name for l in ligs]}")
        else:
            print("          ligands: NONE -- pocket derivation needs a ligand here")
            problems += 1
        if not ok:
            problems += 1
            print("          -> check the motif assignment above before you "
                  "blame the thresholds in kinase_state.DFG_THRESHOLDS")
        print()

    if problems:
        print(f"{problems} problem(s). Fix these before scoring anything.")
    else:
        print("references look good -- go ahead and score")
    return problems


# --------------------------------------------------------------- mode: score

def score_grid(cfg: dict, runs_dir: str, out_csv: str) -> list:
    ref_target = load_structure(cfg["ref_target_state"])
    ref_other = load_structure(cfg["ref_other_state"])

    target_state = cfg["target_state"]
    other_state = cfg["other_state"]

    # Autodetect motifs ONCE, on the crystal structure, then transfer -- or
    # take them from the config's motif_overrides when set (needed for
    # DFG-out references with a broken salt bridge, e.g. 1IEP).
    # Targets with "order_parameter": "none" (non-kinases) skip this entirely.
    skip_dfg = cfg.get("order_parameter", "dfg") == "none"
    overrides = cfg.get("motif_overrides")
    if skip_dfg:
        ref_motifs = None
        print("order_parameter: none -- kinase DFG check disabled "
              "(state recovery is the whole call)\n")
    else:
        ref_motifs = (motifs_from_overrides(ref_target, overrides) if overrides
                      else find_kinase_motifs(ref_target))
        print(f"reference motifs: {ref_motifs.describe(ref_target)}\n")

    region = set(cfg.get("state_region_seqids") or [])
    if not region and ref_motifs.dfg_asp:
        # default: the DFG motif plus a window of the activation loop
        region = set(range(ref_motifs.dfg_asp - 2, ref_motifs.dfg_asp + 12))
        print(f"no state_region_seqids given; using DFG +/- window: "
              f"{sorted(region)}\n")

    back_pocket = set(cfg.get("back_pocket_seqids") or [])

    manifest_path = os.path.join(runs_dir, "manifest.json")
    arms = {}
    if os.path.exists(manifest_path):
        with open(manifest_path) as fh:
            arms = {a["id"]: a for a in json.load(fh)["arms"]}

    rows = []
    arm_dirs = sorted(
        d for d in glob.glob(os.path.join(runs_dir, "*"))
        if os.path.isdir(d)
    )
    if not arm_dirs:
        print(f"no arm subdirectories under {runs_dir}. Expected one directory "
              f"per arm, each holding that arm's model outputs.")
        return []

    for arm_dir in arm_dirs:
        arm_id = os.path.basename(arm_dir)
        files = _find_prediction_files(arm_dir)
        if not files:
            print(f"  {arm_id}: no structure files found, skipping")
            continue

        for path in files:
            try:
                pred = load_structure(path)
            except Exception as exc:
                print(f"  ! could not read {path}: {exc}")
                continue

            row = {
                "arm": arm_id,
                "sample": os.path.basename(path),
                "description": arms.get(arm_id, {}).get("desc", ""),
            }

            # --- state recovery (the headline)
            rec = state_recovery(
                pred,
                {target_state: ref_target, other_state: ref_other},
                region_seqids={target_state: region, other_state: region},
            )
            row["verdict"] = rec.verdict
            row["target_state_recovered"] = int(rec.verdict == target_state)
            row["margin_A"] = round(rec.margin, 3) if rec.margin != float("inf") else ""
            row["decisive"] = int(rec.decisive)
            for c in rec.comparisons:
                tag = c.ref_name.replace("-", "_")
                row[f"region_rmsd_{tag}"] = round(c.region_rmsd, 3) \
                    if c.region_rmsd is not None else ""
                row[f"tm_{tag}"] = round(c.tm, 4)

            # --- independent check: the DFG order parameter itself
            # (kinase-only; empty for "order_parameter": "none" targets)
            if skip_dfg:
                row["dfg_label"] = ""
                row["dfg_agrees_with_target"] = ""
                row["D1"] = ""
                row["D2"] = ""
            else:
                moved = transfer_motifs(ref_target, ref_motifs, pred)
                dfg = classify_dfg(pred, moved)
                row["dfg_label"] = dfg.label
                row["D1"] = round(dfg.d1, 2) if dfg.d1 is not None else ""
                row["D2"] = round(dfg.d2, 2) if dfg.d2 is not None else ""
                row["dfg_agrees_with_target"] = int(dfg.label == target_state)

            # --- pose accuracy, reported separately and never mixed in
            d = ligand_centroid_distance(pred, ref_target)
            row["ligand_centroid_dist_A"] = round(d, 2) if d is not None else ""

            cr = contact_recovery(pred, ref_target)
            row["contact_recall"] = round(cr["recall"], 3) \
                if cr["recall"] is not None else ""

            if back_pocket:
                bp = contact_recovery(pred, ref_target, restrict_to=back_pocket)
                row["back_pocket_recall"] = round(bp["recall"], 3) \
                    if bp["recall"] is not None else ""

            rows.append(row)

    if not rows:
        print("nothing scored")
        return []

    fieldnames = list({k for r in rows for k in r})
    order = ["arm", "sample", "verdict", "target_state_recovered", "margin_A",
             "decisive", "dfg_label", "dfg_agrees_with_target", "D1", "D2",
             "ligand_centroid_dist_A", "contact_recall", "back_pocket_recall"]
    fieldnames = [f for f in order if f in fieldnames] + \
                 [f for f in sorted(fieldnames) if f not in order]

    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {len(rows)} scored samples to {out_csv}")

    summarise(rows, target_state)
    return rows


def summarise(rows: list, target_state: str) -> None:
    """The table that goes on the slide."""
    by_arm = defaultdict(list)
    for r in rows:
        by_arm[r["arm"]].append(r)

    print("\n" + "=" * 92)
    print(f"STATE RECOVERY  (target state: {target_state})")
    print("=" * 92)
    print(f"{'arm':<26s} {'n':>3s} {'state':>7s} {'DFG':>7s} {'agree':>6s} "
          f"{'decisive':>9s} {'lig dist':>9s} {'contact':>8s}")
    print("-" * 92)

    for arm in sorted(by_arm):
        rs = by_arm[arm]
        n = len(rs)
        state_rate = sum(r["target_state_recovered"] for r in rs) / n
        dfg_vals = [r["dfg_agrees_with_target"] for r in rs
                    if r["dfg_agrees_with_target"] != ""]
        dfg_rate = (sum(dfg_vals) / len(dfg_vals)) if dfg_vals else None
        if dfg_rate is None:
            agree = None
        else:
            agree = sum(int(r["target_state_recovered"] ==
                            r["dfg_agrees_with_target"])
                        for r in rs) / n
        decisive = sum(r["decisive"] for r in rs) / n

        dists = [r["ligand_centroid_dist_A"] for r in rs
                 if r["ligand_centroid_dist_A"] != ""]
        recs = [r["contact_recall"] for r in rs if r["contact_recall"] != ""]
        md = f"{sum(dists)/len(dists):.2f}" if dists else "   -"
        mr = f"{sum(recs)/len(recs):.2f}" if recs else "   -"

        dfg_s = f"{dfg_rate:>6.0%}" if dfg_rate is not None else "     -"
        agree_s = f"{agree:>5.0%}" if agree is not None else "    -"
        print(f"{arm:<26s} {n:>3d} {state_rate:>6.0%} {dfg_s} "
              f"{agree_s} {decisive:>8.0%} {md:>9s} {mr:>8s}")

    print("-" * 92)
    print("state   = fraction landing in the target state by dual-reference RMSD")
    print("DFG     = fraction landing there by the independent D1/D2 order parameter")
    print("agree   = fraction where the two metrics agree. If this is low, say so on")
    print("          the slide -- it means the state call itself is shaky, and that")
    print("          is a finding, not something to hide")
    print("lig dist= mean ligand centroid distance, A. Report it NEXT TO state")
    print("          recovery, never instead of it")
    print("=" * 92)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["validate", "score"])
    ap.add_argument("--config", required=True)
    ap.add_argument("--target", required=True)
    ap.add_argument("--runs", default=None)
    ap.add_argument("--out", default="results.csv")
    args = ap.parse_args()

    with open(args.config) as fh:
        cfg = json.load(fh)[args.target]

    if args.mode == "validate":
        raise SystemExit(1 if validate_references(cfg) else 0)

    if not args.runs:
        ap.error("--runs is required for score mode")
    score_grid(cfg, args.runs, args.out)


if __name__ == "__main__":
    main()
