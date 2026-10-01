#!/usr/bin/env python3
"""
make_figure_day1.py — final hackathon figure (figure_day1.png).

Panels:
  a) H2 MSA-depth cliff: folded structures and mean pLDDT vs MSA depth (0 H12-out everywhere).
  b) H3 four-arm H12-junction pLDDT (seed means) with the classifier threshold.
     NOTE: panel b shows the junction window (479-486) pending tonight's re-extraction
     over the challenge H12 window (501-507); swap H3_SEED_MEANS via --h3-csv then.
  c) Pool classifier map (auto-added when --pool-csv pool_state_calls.csv exists):
     junction pLDDT x junction margin per compound, colored by state call.

Usage:
  python make_figure_day1.py --out figure_day1.png
  python make_figure_day1.py --pool-csv classifier_out/pool_state_calls.csv --out figure_day1.png
  python make_figure_day1.py --h3-csv h3_seed_means.csv --out figure_day1.png   # re-cut on 501-507
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
plt.rcParams["svg.fonttype"] = "none"

# Okabe-Ito colorblind-safe
C = {"ago": "#0072B2", "inv": "#D55E00", "apo": "#009E73", "unc": "#999999", "bar": "#56B4E9"}

# --- real data from the H2 titration (titration_summary.md) ---
H2_DEPTHS = [8, 32, 128, 512]
H2_FOLDED_OF_50 = [0, 0, 50, 50]
H2_MEAN_PLDDT = [46.5, 44.5, 89.6, 89.3]

# --- real data from H3 (seed-level mean pLDDT over the junction window 479-486) ---
H3_SEED_MEANS = {
    "apo": [86.9, 87.3, 88.3, 87.9, 87.1],
    "C2\ninverse agonist": [85.2, 81.2, 81.9, 84.2, 83.4],
    "C1\nagonist·sterol": [75.9, 79.0, 77.5, 78.8, 76.3],
    "C3\nagonist·synthetic": [76.8, 78.5, 76.5, 78.6, 74.6],
}
CLASSIFIER_THRESHOLD = 80.0


def panel_a(ax):
    x = np.arange(len(H2_DEPTHS))
    ax.bar(x, np.array(H2_FOLDED_OF_50) / 50 * 100, width=0.55, color=C["bar"],
           label="Folded structures (% of 50)")
    ax.set_xticks(x)
    ax.set_xticklabels([str(d) for d in H2_DEPTHS])
    ax.set_xlabel("MSA depth (sequences)")
    ax.set_ylabel("Folded (pLDDT ≥ 70), %")
    ax.set_ylim(0, 109)
    ax2 = ax.twinx()
    ax2.plot(x, H2_MEAN_PLDDT, "o-", color=C["inv"], lw=2, label="Mean pLDDT")
    ax2.axhline(70, color="grey", ls=":", lw=1)
    ax2.set_ylabel("Mean pLDDT", color=C["inv"])
    ax2.set_ylim(0, 100)
    ax2.tick_params(axis="y", labelcolor=C["inv"])
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8, loc="upper left")
    ax.set_title("a  MSA-depth cliff: fold or agonist, nothing between", fontsize=10, loc="left")


def panel_b(ax, h3_means):
    order = list(h3_means.keys())
    colors = [C["apo"] if k == "apo" else C["inv"] if "inverse" in k.lower() else C["ago"]
              for k in order]
    for i, k in enumerate(order):
        vals = h3_means[k]
        ax.scatter(np.full(len(vals), i) + np.linspace(-0.12, 0.12, len(vals)), vals,
                   s=28, color=colors[i], zorder=3)
        ax.hlines(np.mean(vals), i - 0.22, i + 0.22, color="black", lw=1.6, zorder=4)
    ax.axhline(CLASSIFIER_THRESHOLD, color="grey", ls="--", lw=1)
    ax.text(len(order) - 0.45, CLASSIFIER_THRESHOLD + 0.4, "classifier threshold",
            fontsize=7.5, color="grey", ha="right")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, fontsize=8)
    ax.set_ylabel("H12-junction pLDDT (479–486), seed means")
    ax.set_ylim(70, 92)
    ax.set_title("b  H3: confidence partitions by pharmacological class", fontsize=10, loc="left")


def panel_c(ax, pool_csv):
    df = pd.read_csv(pool_csv)
    df = df[df["call"] != "unscorable"]
    cmap = {"agonist-like": C["ago"], "antagonist-like": C["inv"], "uncertain": C["unc"]}
    for call, grp in df.groupby("call"):
        ax.scatter(grp["jx_margin"], grp["jx_plddt"], s=10, alpha=0.6,
                   color=cmap.get(grp["call"].iloc[0], C["unc"]), label=f"{call} (n={len(grp)})")
    ax.axhline(CLASSIFIER_THRESHOLD, color="grey", ls="--", lw=0.8)
    ax.axvline(5.8, color="grey", ls=":", lw=0.8)
    ax.set_xlabel("Junction RMSD margin (out − in), Å")
    ax.set_ylabel("Junction pLDDT (479–486)")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("c  Pool state calls from the confidence classifier", fontsize=10, loc="left")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool-csv", help="classifier_out/pool_state_calls.csv (adds panel c)")
    ap.add_argument("--h3-csv", help="CSV with columns arm,seed_mean to replace panel b data")
    ap.add_argument("--out", default="figure_day1.png")
    args = ap.parse_args()

    h3_means = H3_SEED_MEANS
    if args.h3_csv:
        h3df = pd.read_csv(args.h3_csv)
        h3_means = {k: g["seed_mean"].tolist() for k, g in h3df.groupby("arm")}

    n_panels = 3 if args.pool_csv else 2
    fig, axes = plt.subplots(1, n_panels, figsize=(4.2 * n_panels, 3.8))
    axes = np.atleast_1d(axes)
    panel_a(axes[0])
    panel_b(axes[1], h3_means)
    if args.pool_csv:
        panel_c(axes[2], args.pool_csv)
    fig.tight_layout()
    out = Path(args.out)
    fig.savefig(out, dpi=220)
    fig.savefig(out.with_suffix(".svg"))
    print(f"Wrote {out} (+ .svg)")


if __name__ == "__main__":
    main()
