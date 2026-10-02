#!/usr/bin/env python3
"""H4 two-arm strip/dot plot: per-structure peptide RMSD, D1 vs D2.
Reads results/h4_per_seed_rmsd.csv (already on disk), writes PNG+SVG.
"""
from __future__ import annotations

import csv
import random
import statistics as stats
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULTS = Path(__file__).parent.parent / "results"


def iqr(vals: list[float]) -> tuple[float, float, float]:
    vals = sorted(vals)
    q1 = stats.median(vals[: len(vals) // 2]) if len(vals) >= 2 else vals[0]
    q3 = stats.median(vals[(len(vals) + 1) // 2 :]) if len(vals) >= 2 else vals[0]
    return q1, stats.median(vals), q3


def main() -> int:
    rows = list(csv.DictReader(open(RESULTS / "h4_per_seed_rmsd.csv")))
    arms = ["D1_agonist", "D2_inverse_agonist"]
    labels = ["D1 (agonist, HC2)", "D2 (inverse agonist, 4P1)"]
    data = {a: [float(r["rmsd"]) for r in rows if r["arm"] == a] for a in arms}

    for a in arms:
        q1, med, q3 = iqr(data[a])
        print(f"{a}: n={len(data[a])} median={med:.2f} Å  IQR=[{q1:.2f}, {q3:.2f}]")

    needs_log = (max(data[arms[0]]) / min(data[arms[1]])) > 50
    print(f"log y-axis: {needs_log}")

    fig, ax = plt.subplots(figsize=(5.5, 4.5), dpi=300)
    rng = random.Random(20261002)
    colors = ["#4C72B0", "#C44E52"]
    for i, (a, lbl, color) in enumerate(zip(arms, labels, colors)):
        ys = data[a]
        xs = [i + 1 + rng.uniform(-0.12, 0.12) for _ in ys]
        ax.scatter(xs, ys, s=22, alpha=0.65, color=color, edgecolors="none", zorder=2)
        q1, med, q3 = iqr(ys)
        ax.plot([i + 1 - 0.22, i + 1 + 0.22], [med, med], color="black", lw=1.8, zorder=3)
        ax.plot([i + 1 - 0.14, i + 1 + 0.14], [q1, q1], color="black", lw=1.0, zorder=3)
        ax.plot([i + 1 - 0.14, i + 1 + 0.14], [q3, q3], color="black", lw=1.0, zorder=3)

    if needs_log:
        ax.set_yscale("log")
    ax.set_xticks([1, 2])
    ax.set_xticklabels(labels)
    ax.set_xlim(0.5, 2.5)
    ax.set_ylabel("Peptide heavy-atom RMSD to 3KYT pose (Å)")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")
    fig.tight_layout()

    png = RESULTS / "h4_dot_plot.png"
    svg = RESULTS / "h4_dot_plot.svg"
    fig.savefig(png, dpi=300, facecolor="white")
    fig.savefig(svg, facecolor="white")
    print(f"wrote {png}")
    print(f"wrote {svg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
