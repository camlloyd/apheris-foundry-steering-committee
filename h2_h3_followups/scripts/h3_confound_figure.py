#!/usr/bin/env python3
"""H3 confound-resolution figure: H12-local pLDDT, four arms (apo, C1, C2,
C3), colored by pharmacological class rather than rigidity -- the point
being that C3 (rigid, pocket-filling, like C2) lands with C1 (its own
pharmacological class) rather than with C2. Reads the already-committed
raw_runs/h3/state_calls_pool_c1-c6.csv, no new predictions or scoring.
"""
from __future__ import annotations

import csv
import random
import statistics as stats
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).parent.parent
SRC = ROOT / "raw_runs" / "h3" / "state_calls_pool_c1-c6.csv"
OUT = ROOT / "results"


def arm_of(path: str) -> str | None:
    for a in ("apo", "c1", "c2", "c3"):
        if f"/{a}/" in path:
            return a
    return None


def iqr(vals: list[float]) -> tuple[float, float, float]:
    vals = sorted(vals)
    q1 = stats.median(vals[: len(vals) // 2])
    q3 = stats.median(vals[(len(vals) + 1) // 2 :])
    return q1, stats.median(vals), q3


def main() -> int:
    rows = list(csv.DictReader(open(SRC)))
    data: dict[str, list[float]] = {a: [] for a in ("apo", "c1", "c2", "c3")}
    for r in rows:
        a = arm_of(r["file"])
        if a:
            data[a].append(float(r["h12_local_plddt"]))

    arms = ["apo", "c1", "c2", "c3"]
    labels = ["apo", "C1\n(agonist,\nflexible)", "C2\n(inverse agonist,\nrigid)", "C3\n(agonist,\nrigid)"]
    # color by pharmacological class, not rigidity -- the whole point of the figure
    colors = {"apo": "#888888", "c1": "#4C72B0", "c2": "#C44E52", "c3": "#4C72B0"}

    for a in arms:
        q1, med, q3 = iqr(data[a])
        print(f"{a}: n={len(data[a])} median={med:.2f} IQR=[{q1:.2f}, {q3:.2f}]")

    fig, ax = plt.subplots(figsize=(6.5, 4.5), dpi=300)
    rng = random.Random(20261002)
    for i, a in enumerate(arms):
        ys = data[a]
        xs = [i + 1 + rng.uniform(-0.12, 0.12) for _ in ys]
        ax.scatter(xs, ys, s=20, alpha=0.6, color=colors[a], edgecolors="none", zorder=2)
        q1, med, q3 = iqr(ys)
        ax.plot([i + 1 - 0.22, i + 1 + 0.22], [med, med], color="black", lw=1.8, zorder=3)
        ax.plot([i + 1 - 0.14, i + 1 + 0.14], [q1, q1], color="black", lw=1.0, zorder=3)
        ax.plot([i + 1 - 0.14, i + 1 + 0.14], [q3, q3], color="black", lw=1.0, zorder=3)

    # bracket highlighting C1 vs C3 (p=0.60, indistinguishable) vs C2 (p=0.0079 from both)
    y_bracket = max(max(v) for v in data.values()) + 2.5
    ax.plot([2, 2, 4, 4], [y_bracket, y_bracket + 0.6, y_bracket + 0.6, y_bracket],
             color="#4C72B0", lw=1.2)
    ax.text(3, y_bracket + 0.9, "p = 0.60 (indistinguishable)", ha="center", fontsize=8, color="#4C72B0")

    ax.set_xticks([1, 2, 3, 4])
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_xlim(0.5, 4.5)
    ax.set_ylabel("H12-local pLDDT (residues 479–486)")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")
    fig.tight_layout()

    png = OUT / "h3_confound_figure.png"
    svg = OUT / "h3_confound_figure.svg"
    fig.savefig(png, dpi=300, facecolor="white")
    fig.savefig(svg, facecolor="white")
    print(f"wrote {png}")
    print(f"wrote {svg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
