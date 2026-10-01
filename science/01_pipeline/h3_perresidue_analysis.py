#!/usr/bin/env python3
"""
h3_perresidue_analysis.py — discriminating analysis for the H3 mechanism question.

Question: is the agonist-associated H12-confidence suppression specific to the
state-defining region (pharmacology-linked), or generic across all ligand-contact
residues (contact-perturbation artifact)?

Reads the existing per-residue pLDDT extraction (h3_perresidue.pkl) for the four H3 arms
(apo, C1, C2, C3) and computes, per residue:
    suppression_vs_inverse = mean(C1, C3) - C2      (agonist mean minus inverse agonist)
    suppression_vs_apo     = mean(C1, C3) - apo
Negative values = agonists LOWER confidence than the comparison arm (the H3 pattern).

Regions:
  * H12       = 501-507 (challenge definition; disordered in 4ZJW)
  * junction  = 479-486 (kit auto-derived window; modeled in both references)
  * contact   = ligand-contact residues, supplied via --contacts (one residue number per line)
  * other     = everything else

Interpretation guide (printed with results):
  * suppression localized to H12/junction only      -> pharmacology-linked, state-adjacent
  * suppression at ALL contact residues similarly   -> contact-perturbation artifact
  * suppression everywhere including non-contacts   -> global effect, not ligand-specific

Usage:
  python h3_perresidue_analysis.py --pkl h3_perresidue.pkl --contacts pocket_contacts.txt --out perresidue_out/

The pkl is expected to be either a long DataFrame (columns: residue, arm, plddt — one row per
residue x arm x structure) or a wide DataFrame/dict (index=residue, one column per arm with
per-residue means). The loader prints what it finds and exits rather than guessing.
"""

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

H12_WINDOW = (501, 507)
JUNCTION_WINDOW = (479, 486)
ARMS = ["apo", "C1", "C2", "C3"]


def load_perresidue(path: str) -> pd.DataFrame:
    """Return long DataFrame with columns: residue (int), arm (str), plddt (float)."""
    with open(path, "rb") as fh:
        obj = pickle.load(fh)
    if isinstance(obj, dict):
        try:
            obj = pd.DataFrame(obj)
        except Exception:
            sys.exit(f"ERROR: pkl contains a dict with keys {list(obj.keys())} — "
                     "cannot coerce to DataFrame. Inspect and adapt the loader.")
    if not isinstance(obj, pd.DataFrame):
        sys.exit(f"ERROR: pkl contains {type(obj)} — expected DataFrame or dict. "
                 "Inspect and adapt the loader.")
    df = obj.copy()
    df.columns = [str(c) for c in df.columns]
    lower = {c.lower(): c for c in df.columns}
    if {"residue", "arm", "plddt"} <= set(lower):
        out = df[[lower["residue"], lower["arm"], lower["plddt"]]].copy()
        out.columns = ["residue", "arm", "plddt"]
        return out
    arm_cols = [c for c in df.columns if c in ARMS]
    if arm_cols:
        res_col = lower.get("residue") or lower.get("resid") or lower.get("res")
        if res_col is None:
            df = df.reset_index().rename(columns={df.index.name or "index": "residue"})
            res_col = "residue"
        out = df.melt(id_vars=[res_col], value_vars=arm_cols,
                      var_name="arm", value_name="plddt")
        out = out.rename(columns={res_col: "residue"})
        return out
    sys.exit(f"ERROR: could not parse pkl schema. Columns found: {list(df.columns)}. "
             "Expected long (residue/arm/plddt) or wide (one column per arm).")


def region_of(res: int, contacts: set) -> str:
    if H12_WINDOW[0] <= res <= H12_WINDOW[1]:
        return "H12"
    if JUNCTION_WINDOW[0] <= res <= JUNCTION_WINDOW[1]:
        return "junction"
    if res in contacts:
        return "contact"
    return "other"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pkl", required=True)
    ap.add_argument("--contacts", help="text file, one ligand-contact residue number per line")
    ap.add_argument("--out", default="perresidue_out")
    args = ap.parse_args()

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    contacts = set()
    if args.contacts:
        contacts = {int(x) for x in Path(args.contacts).read_text().split()}

    df = load_perresidue(args.pkl)
    per_arm = df.groupby(["residue", "arm"])["plddt"].mean().unstack("arm")
    missing = [a for a in ARMS if a not in per_arm.columns]
    if missing:
        sys.exit(f"ERROR: arms missing from pkl: {missing}. Arms found: {list(per_arm.columns)}")

    per_arm["suppression_vs_inverse"] = per_arm[["C1", "C3"]].mean(axis=1) - per_arm["C2"]
    per_arm["suppression_vs_apo"] = per_arm[["C1", "C3"]].mean(axis=1) - per_arm["apo"]
    per_arm["region"] = [region_of(r, contacts) for r in per_arm.index]

    table_path = outdir / "perresidue_suppression.csv"
    per_arm.to_csv(table_path)
    print(f"Wrote {table_path} ({len(per_arm)} residues)\n")

    summary = per_arm.groupby("region")[["suppression_vs_inverse", "suppression_vs_apo"]].agg(
        ["mean", "median", "count"])
    print("Suppression by region (negative = agonists lower than comparison arm):")
    print(summary.round(2).to_string())

    h12jx = per_arm[per_arm["region"].isin(["H12", "junction"])]["suppression_vs_inverse"].dropna()
    contact = per_arm[per_arm["region"] == "contact"]["suppression_vs_inverse"].dropna()
    other = per_arm[per_arm["region"] == "other"]["suppression_vs_inverse"].dropna()
    print("\n--- interpretation ---")
    print(f"H12+junction: mean {h12jx.mean():.2f} (n={len(h12jx)}) | "
          f"contact: {contact.mean() if len(contact) else float('nan'):.2f} (n={len(contact)}) | "
          f"other: {other.mean():.2f} (n={len(other)})")
    if len(contact) >= 3:
        try:
            from scipy.stats import mannwhitneyu
            u, p = mannwhitneyu(h12jx, contact, alternative="two-sided")
            print(f"Mann-Whitney H12+junction vs contact: U={u:.1f}, p={p:.4f} "
                  "(small n — treat as descriptive)")
        except ImportError:
            print("scipy not available — skipping Mann-Whitney (descriptive only).")
    else:
        print("Too few contact residues for a test — descriptive only.")

    # slide-ready strip plot
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
    colors = {"H12": "#D55E00", "junction": "#E69F00", "contact": "#0072B2", "other": "#999999"}
    fig, ax = plt.subplots(figsize=(9, 3.6))
    for region, grp in per_arm.groupby("region"):
        ax.scatter(grp.index, grp["suppression_vs_inverse"], s=14, alpha=0.75,
                   color=colors[region], label=region)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("Residue")
    ax.set_ylabel("Agonist suppression of pLDDT\n(mean(C1,C3) − C2)")
    ax.legend(frameon=False, ncol=4, fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "perresidue_suppression.png", dpi=200)
    fig.savefig(outdir / "perresidue_suppression.svg")
    print(f"\nWrote {outdir/'perresidue_suppression.png'} (+ .svg)")


if __name__ == "__main__":
    main()
