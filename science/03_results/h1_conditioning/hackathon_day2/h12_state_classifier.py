#!/usr/bin/env python3
"""
h12_state_classifier.py — RORγ H12 state classifier for the Apheris Foundry hackathon pool.

Turns the H3 finding (H12-local confidence is ligand-dependent and partitions by
pharmacological class) into per-compound state calls for the deposited compound pool.

Two windows are tracked (see report_additions.md for the numbering resolution):
  * H12      = residues 501-507 — the challenge-slide definition; disordered (unmodeled)
               in the inverse-agonist reference 4ZJW. pLDDT only; no RMSD-to-4ZJW possible.
  * JUNCTION = residues 479-486 — the kit auto-derived window; modeled in BOTH references
               but in different positions, so it supports both pLDDT and RMSD-margin features.

Features per compound (median across seeds/samples):
  * h12_plddt      : mean pLDDT over 501-507          (primary; recalibrate tonight)
  * jx_plddt       : mean pLDDT over 479-486          (H3-calibrated: agonist ~77.0-77.5,
                     inverse agonist ~83.2, apo ~87.5)
  * jx_margin      : RMSD(jx->4ZJW) - RMSD(jx->3KYT)  (H3: agonist ~+7.1, inverse ~+4.4)

Decision rule (conservative, high-precision; calibration n=3 ligands):
  each feature votes agonist-like / antagonist-like; calls require >=2 agreeing votes,
  otherwise "uncertain". Structures below the global pLDDT floor are recorded as
  "unscorable" and listed, never silently dropped (project convention from H2).

Usage:
  python h12_state_classifier.py --scores pool_scores.csv --out classifier_out/
  python h12_state_classifier.py --scores pool_scores.csv --labels deposited_labels.csv --out classifier_out/
  python h12_state_classifier.py --calibrate h3_arm_scores.csv --out classifier_out/   # derive thresholds from H3 arms

scores CSV: one row per STRUCTURE (seed/sample), column names are matched case-insensitively;
adjust COLUMN_CANDIDATES below if the VM schema differs. The script fails loudly and prints
the columns it found rather than guessing silently.

labels CSV (optional): compound_id, true_state   with true_state in {agonist, antagonist}
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------- configuration ----------------

H12_WINDOW = (501, 507)        # challenge-slide H12; pLDDT feature only
JUNCTION_WINDOW = (479, 486)   # kit window; pLDDT + RMSD-margin features

GLOBAL_PLDDT_FLOOR = 70.0      # below this a structure is "unscorable" (H2 convention)

# Defaults from the H3 four-arm experiment, measured on the JUNCTION window.
# h12_plddt has NO calibration yet — recalibrate from the H3 arms (--calibrate)
# after re-extracting per-residue pLDDT over 501-507 from the existing outputs.
DEFAULT_THRESHOLDS = {
    "h12_plddt": None,   # set by --calibrate; if None, this feature is abstained
    "jx_plddt": 80.0,    # midpoint of agonist cluster (77.0-77.5) and inverse agonist (83.2)
    "jx_margin": 5.8,    # midpoint of agonist (+7.1) and inverse agonist (+4.4)
}

COLUMN_CANDIDATES = {
    "compound_id": ["compound_id", "compound", "ligand_id", "ligand", "mol_id", "name"],
    "structure_id": ["structure_id", "structure", "model_id", "sample_id", "pred_id"],
    "global_plddt": ["global_plddt", "plddt", "mean_plddt", "plddt_global"],
    "h12_plddt": ["h12_plddt", "plddt_501_507", "h12_local_plddt"],
    "jx_plddt": ["jx_plddt", "plddt_479_486", "h12_kit_plddt", "junction_plddt"],
    "jx_rmsd_in": ["jx_rmsd_in", "rmsd_479_486_to_3kyt", "rmsd_in", "rmsd_to_3KYT"],
    "jx_rmsd_out": ["jx_rmsd_out", "rmsd_479_486_to_4zjw", "rmsd_out", "rmsd_to_4ZJW"],
    "arm": ["arm", "condition", "group"],  # only used by --calibrate
}


def resolve_columns(df: pd.DataFrame) -> dict:
    lower = {c.lower(): c for c in df.columns}
    resolved = {}
    for key, candidates in COLUMN_CANDIDATES.items():
        for cand in candidates:
            if cand.lower() in lower:
                resolved[key] = lower[cand.lower()]
                break
    return resolved


def load_scores(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    cols = resolve_columns(df)
    required = ["compound_id", "global_plddt", "jx_plddt", "jx_rmsd_in", "jx_rmsd_out"]
    missing = [k for k in required if k not in cols]
    if missing:
        sys.exit(
            "ERROR: could not find columns for "
            f"{missing}.\nColumns present in {path}: {list(df.columns)}\n"
            "Edit COLUMN_CANDIDATES at the top of this script to match the VM schema."
        )
    out = pd.DataFrame({k: df[v] for k, v in cols.items() if k in cols})
    out["jx_margin"] = out["jx_rmsd_out"] - out["jx_rmsd_in"]
    return out


def aggregate_per_compound(df: pd.DataFrame) -> pd.DataFrame:
    """Median across seeds/samples per compound; unscorable structures listed, not dropped."""
    rows = []
    for comp, grp in df.groupby("compound_id"):
        scorable = grp[grp["global_plddt"] >= GLOBAL_PLDDT_FLOOR]
        n_unscorable = int((grp["global_plddt"] < GLOBAL_PLDDT_FLOOR).sum())
        if scorable.empty:
            rows.append({"compound_id": comp, "n_structures": len(grp),
                         "n_unscorable": n_unscorable, "call": "unscorable"})
            continue
        row = {"compound_id": comp, "n_structures": len(grp), "n_unscorable": n_unscorable}
        for feat in ["h12_plddt", "jx_plddt", "jx_margin", "global_plddt"]:
            if feat in scorable:
                row[feat] = float(scorable[feat].median())
        rows.append(row)
    return pd.DataFrame(rows)


def classify(df: pd.DataFrame, thresholds: dict) -> pd.DataFrame:
    calls, votes_detail = [], []
    for _, r in df.iterrows():
        if r.get("call") == "unscorable":
            calls.append("unscorable")
            votes_detail.append("")
            continue
        votes = []
        if thresholds.get("h12_plddt") is not None and not pd.isna(r.get("h12_plddt")):
            votes.append("antagonist" if r["h12_plddt"] >= thresholds["h12_plddt"] else "agonist")
        if not pd.isna(r.get("jx_plddt")):
            votes.append("antagonist" if r["jx_plddt"] >= thresholds["jx_plddt"] else "agonist")
        if not pd.isna(r.get("jx_margin")):
            votes.append("antagonist" if r["jx_margin"] <= thresholds["jx_margin"] else "agonist")
        n_ant = votes.count("antagonist")
        n_ago = votes.count("agonist")
        if n_ant >= 2 and n_ant > n_ago:
            call = "antagonist-like"
        elif n_ago >= 2 and n_ago > n_ant:
            call = "agonist-like"
        elif len(votes) == 1:           # only one feature available: low-confidence call
            call = f"{votes[0]}-like (single-feature)"
        else:
            call = "uncertain"
        calls.append(call)
        votes_detail.append("/".join(votes))
    df = df.copy()
    df["call"] = calls
    df["votes"] = votes_detail
    return df


def calibrate(h3_csv: str) -> dict:
    """Derive thresholds from the H3 arms (known pharmacology): C1/C3 = agonist, C2 = antagonist."""
    df = load_scores(h3_csv)
    if "arm" not in df.columns:
        sys.exit("ERROR: --calibrate needs an 'arm' column (apo/C1/C2/C3) in the H3 scores CSV.")
    per_arm = df.groupby("arm")[["h12_plddt", "jx_plddt", "jx_margin"]].median()
    ago = per_arm.loc[["C1", "C3"]]
    inv = per_arm.loc[["C2"]]
    thr = {}
    for feat, direction in [("h12_plddt", "high_is_ant"), ("jx_plddt", "high_is_ant"),
                            ("jx_margin", "low_is_ant")]:
        if per_arm[feat].isna().all():
            thr[feat] = None
            continue
        a, i = float(ago[feat].mean()), float(inv[feat].mean())
        thr[feat] = round((a + i) / 2.0, 3)
        if (direction == "high_is_ant" and i < a) or (direction == "low_is_ant" and i > a):
            print(f"WARNING: {feat} direction flipped vs H3 expectation "
                  f"(agonist={a}, inverse={i}) — inspect before trusting this threshold.")
    print("Calibrated thresholds from H3 arms:", json.dumps(thr, indent=2))
    return thr


def score_against_labels(calls: pd.DataFrame, labels_path: str) -> dict:
    labels = pd.read_csv(labels_path)
    labels.columns = [c.lower() for c in labels.columns]
    merged = calls.merge(labels, on="compound_id", how="inner")
    if merged.empty:
        sys.exit("ERROR: no compound_id overlap between calls and labels.")
    merged["true_simple"] = merged["true_state"].str.lower().str.contains("ant").map(
        {True: "antagonist", False: "agonist"})
    decided = merged[~merged["call"].isin(["unscorable", "uncertain"])].copy()
    decided["pred_simple"] = decided["call"].str.contains("ant").map(
        {True: "antagonist", False: "agonist"})
    report = {"n_labeled": len(merged), "n_decided": len(decided),
              "n_uncertain": int((merged["call"] == "uncertain").sum()),
              "n_unscorable": int((merged["call"] == "unscorable").sum())}
    if len(decided):
        report["accuracy_decided"] = float((decided["pred_simple"] == decided["true_simple"]).mean())
        for cls in ["agonist", "antagonist"]:
            sub = decided[decided["true_simple"] == cls]
            if len(sub):
                report[f"accuracy_{cls}"] = float((sub["pred_simple"] == cls).mean())
        cm = pd.crosstab(decided["true_simple"], decided["pred_simple"])
        report["confusion_matrix"] = cm.to_dict()
    return report


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", help="per-structure scores CSV from compute-metrics/state_recovery2")
    ap.add_argument("--labels", help="optional deposited-set labels CSV (compound_id, true_state)")
    ap.add_argument("--calibrate", help="H3 arm scores CSV to derive thresholds from (apo/C1/C2/C3)")
    ap.add_argument("--out", default="classifier_out")
    args = ap.parse_args()

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    thresholds = dict(DEFAULT_THRESHOLDS)
    if args.calibrate:
        thresholds.update(calibrate(args.calibrate))

    manifest = {"timestamp": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
                "thresholds": thresholds,
                "windows": {"H12": H12_WINDOW, "junction": JUNCTION_WINDOW},
                "global_plddt_floor": GLOBAL_PLDDT_FLOOR,
                "inputs": {"scores": args.scores, "labels": args.labels,
                           "calibrate": args.calibrate}}

    if args.scores:
        df = load_scores(args.scores)
        per_compound = aggregate_per_compound(df)
        calls = classify(per_compound, thresholds)
        calls_path = outdir / "pool_state_calls.csv"
        calls.to_csv(calls_path, index=False)
        print(f"\nWrote {calls_path}  ({len(calls)} compounds)")
        print(calls["call"].value_counts().to_string())
        if args.labels:
            report = score_against_labels(calls, args.labels)
            manifest["label_report"] = report
            print("\nLabel validation:")
            print(json.dumps(report, indent=2))
    manifest_path = outdir / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"\nWrote {manifest_path}")


if __name__ == "__main__":
    main()
