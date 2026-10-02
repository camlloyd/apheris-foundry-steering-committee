#!/usr/bin/env python3
"""AUROC for predicting call correctness from (a) decisive margin,
(b) mean pLDDT, on the existing 258-compound pool. No new predictions, no
classifier fit, no threshold tuning. Reuses task3_per_compound.csv.
"""
from __future__ import annotations

import csv
import random
from pathlib import Path

RESULTS = Path(__file__).parent.parent / "results"


def auroc(scores: list[float], labels: list[int]) -> float:
    """labels: 1 = correct, 0 = incorrect. Mann-Whitney U form, exact
    (handles ties by giving them 0.5 credit), no external deps."""
    pos = [s for s, y in zip(scores, labels) if y == 1]
    neg = [s for s, y in zip(scores, labels) if y == 0]
    if not pos or not neg:
        return float("nan")
    wins = 0.0
    for p in pos:
        for n in neg:
            if p > n:
                wins += 1
            elif p == n:
                wins += 0.5
    return wins / (len(pos) * len(neg))


def bootstrap_ci(scores: list[float], labels: list[int], n_boot: int = 2000, seed: int = 20261002) -> tuple[float, float, float]:
    rng = random.Random(seed)
    n = len(scores)
    point = auroc(scores, labels)
    boots = []
    for _ in range(n_boot):
        idx = [rng.randrange(n) for _ in range(n)]
        s = [scores[i] for i in idx]
        l = [labels[i] for i in idx]
        a = auroc(s, l)
        if a == a:  # not nan
            boots.append(a)
    boots.sort()
    lo = boots[int(0.025 * len(boots))]
    hi = boots[int(0.975 * len(boots))]
    return point, lo, hi


def pearson(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = sum((x - mx) ** 2 for x in xs) ** 0.5
    dy = sum((y - my) ** 2 for y in ys) ** 0.5
    return num / (dx * dy) if dx > 0 and dy > 0 else float("nan")


def main() -> int:
    rows = list(csv.DictReader(open(RESULTS / "task3_per_compound.csv")))
    usable = [r for r in rows if r["mean_margin"] and r["mean_plddt"]]
    n_dropped = len(rows) - len(usable)
    print(f"compounds: {len(rows)}, usable (have both margin and pLDDT): {len(usable)}, dropped: {n_dropped}")

    labels = [1 if r["correct"] == "True" else 0 for r in usable]
    margins = [float(r["mean_margin"]) for r in usable]
    plddts = [float(r["mean_plddt"]) for r in usable]

    auc_margin, lo_m, hi_m = bootstrap_ci(margins, labels)
    auc_plddt, lo_p, hi_p = bootstrap_ci(plddts, labels)
    corr = pearson(margins, plddts)

    lines = []
    lines.append("# Confidence-signal AUROC — margin vs. pLDDT\n")
    lines.append(f"**Scope:** existing 258-compound pool, no new predictions, no classifier fit.\n")
    lines.append(f"n = {len(usable)} compounds ({n_dropped} dropped: missing margin or pLDDT).\n")
    lines.append("## AUROC for predicting call correctness\n")
    lines.append("| signal | AUROC | 95% CI (bootstrap, n=2000) |")
    lines.append("|---|---|---|")
    lines.append(f"| decisive margin | {auc_margin:.3f} | [{lo_m:.3f}, {hi_m:.3f}] |")
    lines.append(f"| mean pLDDT | {auc_plddt:.3f} | [{lo_p:.3f}, {hi_p:.3f}] |")
    lines.append("")
    lines.append(f"## Margin vs. pLDDT correlation\n")
    lines.append(f"Pearson r = **{corr:.3f}** across the {len(usable)} compounds.")
    if abs(corr) < 0.2:
        lines.append(" Near-independent — consistent with margin and pLDDT carrying "
                      "different information about a call, not the same signal measured twice.")
    elif abs(corr) < 0.5:
        lines.append(" Moderately correlated, not independent — the two signals share real "
                      "information rather than contributing fully separate evidence.")
    else:
        lines.append(" Strongly correlated — treating them as separate evidence on a slide "
                      "would be double-counting.")
    lines.append("")
    lines.append("## Reading this against Task 3\n")
    if auc_plddt >= auc_margin:
        lines.append("pLDDT's AUROC here is **at or above** margin's, which is a weaker story "
                      "for margin-as-the-good-signal than Task 3's quartile table suggested on its "
                      "own (margin monotonic, pLDDT not) — AUROC and a 4-bin monotonicity check "
                      "are different questions, and they don't have to agree. Reported as measured.")
    else:
        lines.append("Consistent with Task 3: margin outperforms raw pLDDT as a correctness signal.")
    out = RESULTS / "confidence_auroc.md"
    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out}")
    print(f"AUROC margin: {auc_margin:.3f} [{lo_m:.3f}, {hi_m:.3f}]")
    print(f"AUROC pLDDT:  {auc_plddt:.3f} [{lo_p:.3f}, {hi_p:.3f}]")
    print(f"corr(margin, pLDDT) = {corr:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
