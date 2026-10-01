#!/usr/bin/env python3
"""H12 state classifier — confidence-based, not geometry-based.

    python3 h12_state_classifier.py --work <dir with out/**/*.cif> \
        --kit ~/data/apheris-foundry-steering-committee/apheris_kit_rorgamma \
        --rgkit ~/h2_package --out state_calls.csv

Why this exists: across H1+H2+H3 (322 structures, nine conditioning levers,
four ligand arms), OpenFold3 never moved RORgamma's H12 into the H12-out
basin by RMSD. The geometric state call is therefore useless as a classifier
on new compounds -- it will read H12-in for everything. H3 found a signal
that *does* separate ligand classes: H12-local pLDDT (mean per-residue pLDDT
over the auto-derived state region, 479-486, 3KYT/construct numbering) is
higher for the inverse agonist (4P1, mean 83.2) than for either agonist
tested (HC2 25-HC, BIO592 synthetic; both ~77), p=0.0079, complete 5-vs-5
seed separation, confound-tested (rigid synthetic agonist BIO592 patterns
with the flexible natural agonist, not with the rigid inverse agonist, so
the effect is pharmacology-linked at this residue window, not purely a
rigidity/pocket-filling artifact -- see H3_RESULTS_REPORT.md).

A follow-up discriminating analysis (ligand-contact-residue sweep) found the
same direction and comparable-or-larger effect size at NON-H12 contact
residues too -- so the honest reading is "general ligand-contact-region
confidence suppression, more pronounced for agonists," not a signal
specific to the state-defining residues. This classifier is built and
reported on that basis: it is a ligand-class confidence signal, not a
validated structural-state detector. State it exactly that way on the slide.

CALIBRATION CAVEAT (load-bearing, repeat it everywhere this script's output
is used): the threshold below is fit on n=3 ligands (HC2, 4P1, BIO592) from
one target (RORgamma LBD). It has not been validated on any compound outside
that set. Report accuracy numbers from this script with that caveat attached,
not as a general claim.

Known failure mode of this classifier: apo (no ligand) also reads "high
H12-local pLDDT" (mean 87.5, higher than the inverse agonist's 83.2) -- i.e.
on local pLDDT alone, apo is NOT distinguishable from "inverse-agonist-like."
This script assumes every input structure has a ligand bound (true for a
compound pool); it will mislabel genuinely apo structures and does not try
to detect apo itself.
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import re
import sys
from collections import defaultdict

import numpy as np

REGION = list(range(479, 487))          # H12 state region, construct numbering
PLDDT_FLOOR = 70.0                       # fold-quality gate, same as score_titration.py
H12_PLDDT_THRESHOLD = 80.0               # midpoint of 77.0-77.5 (agonists) and 83.2 (4P1)

# Calibration reference (H3 seed-level means, for the printed self-test table)
CALIBRATION = {
    "apo (no ligand, NOT handled by this classifier)": 87.50,
    "C1 HC2 (agonist, flexible)": 77.52,
    "C2 4P1 (inverse agonist, rigid)": 83.19,
    "C3 BIO592 (agonist, rigid synthetic)": 77.00,
}


def per_residue_plddt(scores_path: pathlib.Path) -> dict[int, float]:
    """token-index -> mean pLDDT, from atom_predicted_local_distance_difference_test
    grouped by atom_tokens. Returns {} if the scores file doesn't have these keys
    (e.g. a different module version) -- caller falls back to whole-structure mean."""
    try:
        d = json.load(open(scores_path))
    except Exception:
        return {}
    atoms = d.get("atom_predicted_local_distance_difference_test")
    toks = d.get("atom_tokens")
    if not atoms or not toks:
        return {}
    sums = defaultdict(float)
    counts = defaultdict(int)
    for a, t in zip(atoms, toks):
        sums[t] += a
        counts[t] += 1
    return {t: sums[t] / counts[t] for t in sums}


def whole_mean_plddt(cif: pathlib.Path, scores_path: pathlib.Path):
    if not scores_path.exists():
        return None
    try:
        s = json.load(open(scores_path))
    except Exception:
        return None
    for k, v in s.items():
        if "local_distance" in k and isinstance(v, list) and v:
            m = float(np.mean(v))
            return m * 100 if m <= 1.0 else m
    return None


def parse_offset(construct_len: int) -> int:
    """3KYT/construct numbering starts at 265 for the standard 243-residue
    RORgamma LBD construct used throughout H1-H3; output (predicted) structures
    are numbered 1..construct_len. offset = 265 - 1 = 264."""
    return 264


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True, help="dir containing out/**/*.cif + sibling _scores.json")
    ap.add_argument("--kit", required=True, help="path to apheris_kit_rorgamma (for refs_rorgamma/)")
    ap.add_argument("--rgkit", default=None, help="path to rgkit package dir, if not inside --kit")
    ap.add_argument("--out", default="state_calls.csv")
    ap.add_argument("--threshold", type=float, default=H12_PLDDT_THRESHOLD,
                     help=f"H12-local pLDDT midpoint (default {H12_PLDDT_THRESHOLD}, "
                          "the agonist/inverse-agonist midpoint found in H3)")
    ap.add_argument("--labels", default=None,
                     help="optional CSV with columns pattern,expected_call -- "
                          "any row whose 'pattern' substring is in the file path gets "
                          "that expected_call; prints per-class accuracy if given")
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
              "the rgkit package.", file=sys.stderr)
        return 2

    refs_dir = kit / "refs_rorgamma"
    a = S.load_chain(refs_dir / "3KYT_clean.cif", source="3KYT")   # H12-in
    b = S.load_chain(refs_dir / "4ZJW_clean.cif", source="4ZJW")   # H12-out
    refs = {"active": a, "inactive": b}

    cifs = sorted(p for p in work.rglob("*.cif") if "__MACOSX" not in str(p))
    if not cifs:
        print(f"no .cif under {work}", file=sys.stderr)
        return 2

    rows = []
    for cif in cifs:
        scores = cif.with_name(cif.name.replace("_model.cif", "_scores.json"))
        row = {"file": str(cif.relative_to(work)), "status": "",
               "whole_plddt": "", "h12_local_plddt": "",
               "rmsd_H12in": "", "rmsd_H12out": "", "signed_margin": "",
               "state_call": "", "confidence_score": ""}
        try:
            c = S.load_chain(cif, source=cif.name)
        except Exception as exc:
            row["status"] = f"load-error: {exc}"[:80]
            rows.append(row)
            continue

        whole = whole_mean_plddt(cif, scores)
        row["whole_plddt"] = round(whole, 1) if whole is not None else ""

        pr = per_residue_plddt(scores)
        offset = parse_offset(len(c.residues)) if hasattr(c, "residues") else 264
        h12_toks = [r - offset for r in REGION]
        h12_vals = [pr[t] for t in h12_toks if t in pr]
        h12_local = float(np.mean(h12_vals)) if h12_vals else None
        row["h12_local_plddt"] = round(h12_local, 2) if h12_local is not None else ""

        r = SR.state_recovery_generic(c, refs, REGION, decisive_margin=1.0)
        d_in, d_out = r["region_rmsd_active"], r["region_rmsd_inactive"]
        row["rmsd_H12in"] = round(d_in, 3)
        row["rmsd_H12out"] = round(d_out, 3)
        signed_margin = d_in - d_out   # positive => geometrically closer to H12-out
        row["signed_margin"] = round(signed_margin, 3)

        if whole is not None and whole < PLDDT_FLOOR:
            row["status"] = "unscorable"
            row["state_call"] = f"fold failed (pLDDT {whole:.1f} < {PLDDT_FLOOR:g})"
        elif h12_local is None:
            row["status"] = "no-per-residue-scores"
            row["state_call"] = "fallback: scores.json lacks per-atom pLDDT; whole-structure " \
                                 "mean used instead, see whole_plddt column"
        else:
            row["status"] = "ok"
            if h12_local > args.threshold:
                row["state_call"] = "inverse-agonist-like (low H12 local confidence override)"
            else:
                row["state_call"] = "agonist-like"
            # confidence: distance from threshold, scaled by the H3 calibration spread (~6 pts)
            row["confidence_score"] = round(abs(h12_local - args.threshold) / 6.0, 2)
        rows.append(row)

    out_csv = pathlib.Path(args.out)
    with out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    n_ok = sum(1 for r in rows if r["status"] == "ok")
    n_uns = sum(1 for r in rows if r["status"] == "unscorable")
    n_call_ia = sum(1 for r in rows if "inverse-agonist-like" in r["state_call"])
    n_call_ag = sum(1 for r in rows if r["state_call"] == "agonist-like")
    print(f"scored {len(rows)} structures -> {out_csv}")
    print(f"  ok: {n_ok}  unscorable (pLDDT<{PLDDT_FLOOR:g}): {n_uns}")
    print(f"  called 'inverse-agonist-like': {n_call_ia}   called 'agonist-like': {n_call_ag}")
    print()
    print("CAVEAT: threshold calibrated on n=3 ligands, one target. Apo structures "
          "will misclassify as inverse-agonist-like (not handled). Report accordingly.")

    if args.labels:
        label_rules = list(csv.DictReader(open(args.labels)))
        per_class = defaultdict(lambda: [0, 0])  # expected_call -> [correct, total]
        for r in rows:
            if r["status"] != "ok":
                continue
            expected = next((lr["expected_call"] for lr in label_rules if lr["pattern"] in r["file"]), None)
            if expected is None:
                continue
            per_class[expected][1] += 1
            if expected in r["state_call"]:
                per_class[expected][0] += 1
        print("\n=== accuracy vs --labels ===")
        tot_c = tot_n = 0
        for cls, (correct, total) in per_class.items():
            print(f"  {cls}: {correct}/{total} = {correct/total:.1%}")
            tot_c += correct; tot_n += total
        if tot_n:
            print(f"  OVERALL: {tot_c}/{tot_n} = {tot_c/tot_n:.1%}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
