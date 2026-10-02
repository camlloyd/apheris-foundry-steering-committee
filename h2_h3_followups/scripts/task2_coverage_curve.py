#!/usr/bin/env python3
"""Task 2 — accuracy/coverage curve on the existing 258-compound pool.
No new predictions; sweeps the per-compound decisive-margin threshold
already used (>= 1.0 A) and reports accuracy vs. coverage.
"""
from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

CALLS_CSV = Path.home() / "h3" / "official_state_calls_all.csv"


def batch_of(path: str) -> str | None:
    if "pool120" in path:
        return "batch_A_curated_120"
    if "rand150" in path:
        return "batch_B_random_150"
    return None


def true_label(path: str) -> str | None:
    m = re.search(r"(?:POOL|RAND)_[A-Za-z0-9]+_(agonist|antagonist)", path)
    return m.group(1) if m else None


def cid_of(path: str) -> str | None:
    m = re.search(r"(?:POOL|RAND)_([A-Za-z0-9]+)_(?:agonist|antagonist)", path)
    return m.group(1) if m else None


def majority(states: list[str]) -> str | None:
    a, b = states.count("agonist"), states.count("antagonist")
    if a > b:
        return "agonist"
    if b > a:
        return "antagonist"
    return None


def main() -> int:
    rows = list(csv.DictReader(open(CALLS_CSV)))
    pool_rows = [r for r in rows if batch_of(r["file"])]
    by_compound: dict = defaultdict(list)
    for r in pool_rows:
        key = (batch_of(r["file"]), cid_of(r["file"]), true_label(r["file"]))
        by_compound[key].append(r)
    n_total = len(by_compound)
    print(f"pool compounds: {n_total}")

    out_rows = []
    for margin_thresh in [round(0.0 + 0.1 * i, 2) for i in range(0, 41)]:
        correct = 0
        scored = 0
        for (_, _, lbl), rs in by_compound.items():
            dec = [r["state"] for r in rs if r["margin"] and float(r["margin"]) >= margin_thresh]
            if not dec:
                continue  # abstains at this threshold
            call = majority(dec)
            scored += 1
            if call == lbl:
                correct += 1
        coverage = scored / n_total
        acc = correct / scored if scored else float("nan")
        out_rows.append({
            "margin_threshold": margin_thresh,
            "coverage": round(coverage, 4),
            "n_scored": scored,
            "accuracy_on_scored": round(acc, 4) if scored else "",
        })

    out = Path(__file__).parent.parent / "results" / "task2_accuracy_coverage_curve.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    print(f"wrote {out}")

    at_1 = next(r for r in out_rows if r["margin_threshold"] == 1.0)
    at_0 = next(r for r in out_rows if r["margin_threshold"] == 0.0)
    print(f"margin>=0.0 (full coverage): coverage={at_0['coverage']:.1%} accuracy={at_0['accuracy_on_scored']:.1%}")
    print(f"margin>=1.0 (current 'decisive' cut): coverage={at_1['coverage']:.1%} accuracy={at_1['accuracy_on_scored']:.1%}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
