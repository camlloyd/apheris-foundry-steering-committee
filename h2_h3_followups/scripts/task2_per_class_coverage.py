#!/usr/bin/env python3
"""Per-class accuracy/coverage curves (agonist vs. antagonist) on the
existing 258-compound pool. No new predictions, no threshold tuning.
"""
from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

CALLS_CSV = Path.home() / "h3" / "official_state_calls_all.csv"
RESULTS = Path(__file__).parent.parent / "results"


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

    n_by_class = {"agonist": 0, "antagonist": 0}
    for (_, _, lbl) in by_compound:
        n_by_class[lbl] += 1
    print(f"n agonist={n_by_class['agonist']}, n antagonist={n_by_class['antagonist']}")

    thresholds = [round(0.0 + 0.1 * i, 2) for i in range(0, 41)]
    out_rows = []
    crossing = {"agonist": None, "antagonist": None}
    for thresh in thresholds:
        row = {"margin_threshold": thresh}
        for cls in ("agonist", "antagonist"):
            compounds = [(k, rs) for k, rs in by_compound.items() if k[2] == cls]
            scored = 0
            correct = 0
            for key, rs in compounds:
                dec = [r["state"] for r in rs if r["margin"] and float(r["margin"]) >= thresh]
                if not dec:
                    continue
                call = majority(dec)
                scored += 1
                if call == cls:
                    correct += 1
            coverage = scored / n_by_class[cls]
            acc = correct / scored if scored else float("nan")
            row[f"{cls}_coverage"] = round(coverage, 4)
            row[f"{cls}_n_scored"] = scored
            row[f"{cls}_accuracy"] = round(acc, 4) if scored else ""
            if crossing[cls] is None and scored and acc >= 0.80:
                crossing[cls] = (thresh, coverage, acc)
        out_rows.append(row)

    out = RESULTS / "task2_per_class_coverage.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    print(f"wrote {out}")

    for cls in ("agonist", "antagonist"):
        c = crossing[cls]
        if c:
            print(f"{cls}: crosses 80% accuracy at margin>={c[0]} (coverage={c[1]:.1%})")
        else:
            print(f"{cls}: never reaches 80% accuracy in the swept range")

    if crossing["agonist"] and crossing["antagonist"]:
        delta_cov = crossing["antagonist"][1] - crossing["agonist"][1]
        print(f"coverage needed to hit 80%, antagonist minus agonist: {delta_cov:+.1%}")
        if delta_cov < 0:
            print("antagonist reaches 80% at LOWER coverage than agonist -- abstention rescues "
                  "antagonist calls more cheaply, in coverage terms, than agonist ones.")
        else:
            print("antagonist needs MORE coverage given up than agonist to reach 80% -- "
                  "abstention does not rescue antagonist more cheaply than agonist.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
