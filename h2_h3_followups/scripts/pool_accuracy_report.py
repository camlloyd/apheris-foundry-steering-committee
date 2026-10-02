#!/usr/bin/env python3
"""Combine the official-rescore CSV with the two ChEMBL pool batches (120
curated + 150 random) and report accuracy against real ground-truth labels,
per batch and pooled, both all-samples and decisive-margin-only.

    python3 pool_accuracy_report.py --calls h3/official_state_calls_all.csv
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict


def true_label(path: str) -> str | None:
    m = re.search(r"(?:POOL|RAND)_[A-Za-z0-9]+_(agonist|antagonist)", path)
    return m.group(1) if m else None


def batch_of(row_id: str) -> str | None:
    if "pool120" in row_id:
        return "batch_A_curated_120"
    if "rand150" in row_id:
        return "batch_B_random_150"
    return None


def majority_call(states: list[str]) -> str | None:
    """Majority vote over agonist/antagonist calls. Returns None on an exact
    tie rather than silently breaking it — tie-breaking via max(set(...))
    is non-deterministic across runs (Python's str hash randomization),
    which is not acceptable for a number going on a slide."""
    n_ago = states.count("agonist")
    n_anta = states.count("antagonist")
    if n_ago > n_anta:
        return "agonist"
    if n_anta > n_ago:
        return "antagonist"
    return None


def report(rows: list[dict]) -> None:
    by_batch = defaultdict(lambda: defaultdict(list))
    for r in rows:
        b = batch_of(r["file"])
        if b is None:
            continue
        lbl = true_label(r["file"])
        if lbl is None:
            continue
        cid = re.search(r"(?:POOL|RAND)_([A-Za-z0-9]+)_(?:agonist|antagonist)", r["file"]).group(1)
        by_batch[b][(cid, lbl)].append(r)

    all_compounds = defaultdict(list)
    for b, compounds in by_batch.items():
        print(f"\n=== {b} ===")
        correct = {"agonist": [0, 0], "antagonist": [0, 0]}
        correct_dec = {"agonist": [0, 0], "antagonist": [0, 0]}
        ties = tied_dec = 0
        for (cid, lbl), rs in compounds.items():
            call = majority_call([r["state"] for r in rs])
            correct[lbl][1] += 1
            if call is None:
                ties += 1
            elif call == lbl:
                correct[lbl][0] += 1
            dec = [r["state"] for r in rs if r["margin"] and float(r["margin"]) >= 1.0]
            if dec:
                call_dec = majority_call(dec)
                correct_dec[lbl][1] += 1
                if call_dec is None:
                    tied_dec += 1
                elif call_dec == lbl:
                    correct_dec[lbl][0] += 1
            all_compounds[(b, cid, lbl)] = rs

        for lbl, (c, t) in correct.items():
            print(f"  {lbl}: {c}/{t} = {c/t:.1%}" if t else f"  {lbl}: n/a")
        tc, tt = sum(c for c, t in correct.values()), sum(t for c, t in correct.values())
        print(f"  OVERALL (all samples): {tc}/{tt} = {tc/tt:.1%} ({ties} exact ties, counted as wrong)" if tt else "  n/a")

        for lbl, (c, t) in correct_dec.items():
            print(f"  {lbl} (decisive subset): {c}/{t} = {c/t:.1%}" if t else f"  {lbl}: no decisive compounds")
        tc2, tt2 = sum(c for c, t in correct_dec.values()), sum(t for c, t in correct_dec.values())
        print(f"  OVERALL (decisive subset): {tc2}/{tt2} = {tc2/tt2:.1%} ({tied_dec} exact ties, counted as wrong)" if tt2 else "  n/a")

    print("\n=== POOLED (both batches combined) ===")
    correct = {"agonist": [0, 0], "antagonist": [0, 0]}
    ties = 0
    for (b, cid, lbl), rs in all_compounds.items():
        call = majority_call([r["state"] for r in rs])
        correct[lbl][1] += 1
        if call is None:
            ties += 1
        elif call == lbl:
            correct[lbl][0] += 1
    for lbl, (c, t) in correct.items():
        print(f"  {lbl}: {c}/{t} = {c/t:.1%}" if t else f"  {lbl}: n/a")
    tc, tt = sum(c for c, t in correct.values()), sum(t for c, t in correct.values())
    print(f"  OVERALL: {tc}/{tt} = {tc/tt:.1%} ({ties} exact ties, counted as wrong)" if tt else "  n/a")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True)
    args = ap.parse_args()
    rows = list(csv.DictReader(open(args.calls)))
    report(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
