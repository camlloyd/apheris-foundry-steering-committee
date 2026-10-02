#!/usr/bin/env python3
"""Task 1 — decision-boundary recalibration on the existing 258-compound
pool. No new predictions; rescores official_state_calls_all.csv under the
pre-registered rule (a), then (clearly separate, post-hoc) sweeps the
boundary. See ../preregistrations/H5_TASK1_BOUNDARY_PREREGISTRATION.md.
"""
from __future__ import annotations

import csv
import re
import sys
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


def balanced_and_pooled(by_compound: dict, call_fn) -> tuple[float, float, dict]:
    correct = {"agonist": [0, 0], "antagonist": [0, 0]}
    for (_, _, lbl), rs in by_compound.items():
        call = call_fn(rs)
        correct[lbl][1] += 1
        if call == lbl:
            correct[lbl][0] += 1
    per_class = {lbl: (c / t if t else float("nan")) for lbl, (c, t) in correct.items()}
    balanced = sum(per_class.values()) / len(per_class)
    tc = sum(c for c, t in correct.values())
    tt = sum(t for c, t in correct.values())
    pooled = tc / tt
    return balanced, pooled, per_class


def main() -> int:
    rows = list(csv.DictReader(open(CALLS_CSV)))
    pool_rows = [r for r in rows if batch_of(r["file"])]
    by_compound: dict = defaultdict(list)
    for r in pool_rows:
        key = (batch_of(r["file"]), cid_of(r["file"]), true_label(r["file"]))
        by_compound[key].append(r)
    print(f"pool compounds: {len(by_compound)} ({len(pool_rows)} structures)")

    # --- D_ref: mutual H12 displacement between the two references ---
    # Measured once, mechanically, per the pre-registration (expected to be
    # a no-op given the symmetric-metric argument there).
    sys.path.insert(0, "/tmp/rorgt_kit/rorgt_candidate_kit")
    from lib import h12_state
    ag_path, anta_path = h12_state.reference_paths()
    work = Path("/tmp/h4_boundary_calib")
    work.mkdir(exist_ok=True)
    calls = h12_state.call_states(
        {"agonist_ref_as_input": ag_path, "antagonist_ref_as_input": anta_path}, work
    )
    d_ref_from_agonist = calls["agonist_ref_as_input"].rmsd_to_antagonist
    d_ref_from_antagonist = calls["antagonist_ref_as_input"].rmsd_to_agonist
    print(f"D_ref (agonist_ref -> antagonist_ref): {d_ref_from_agonist:.5f} A")
    print(f"D_ref (antagonist_ref -> agonist_ref): {d_ref_from_antagonist:.5f} A")
    print(f"References sit at delta = +-{d_ref_from_agonist:.5f} -> midpoint is delta = 0, "
          f"same as the current rule.")

    # --- current rule (delta > 0 -> agonist, delta threshold = 0) ---
    def call_current(rs):
        states = [r["state"] for r in rs]  # already computed with the current rule
        return majority(states)

    bal0, pooled0, per_class0 = balanced_and_pooled(by_compound, call_current)
    print(f"\ncurrent rule:  balanced={bal0:.1%} pooled={pooled0:.1%} per_class={per_class0}")

    # --- rule (a): recompute delta per-structure, threshold at 0 (== D_ref midpoint) ---
    def call_rule_a(rs):
        deltas = []
        for r in rs:
            if not r["rmsd_to_agonist"] or not r["rmsd_to_antagonist"]:
                continue  # unscored, same as the current rule's "state" column
            to_ag, to_an = float(r["rmsd_to_agonist"]), float(r["rmsd_to_antagonist"])
            deltas.append("agonist" if to_ag < to_an else "antagonist")
        return majority(deltas)

    bal_a, pooled_a, per_class_a = balanced_and_pooled(by_compound, call_rule_a)
    print(f"rule (a):      balanced={bal_a:.1%} pooled={pooled_a:.1%} per_class={per_class_a}")
    print(f"identical to current rule: {abs(bal_a - bal0) < 1e-12 and abs(pooled_a - pooled0) < 1e-12}")

    # --- post-hoc only: sweep the boundary on delta, report accuracy vs boundary ---
    print("\n--- post-hoc boundary sweep (sensitivity analysis, NOT a result) ---")
    sweep_out = Path(__file__).parent.parent / "results" / "task1_boundary_sweep.csv"
    sweep_rows = []
    lo, hi = -6.0, 6.0
    n = 49
    for i in range(n):
        thresh = lo + (hi - lo) * i / (n - 1)

        def call_at(rs, thresh=thresh):
            calls = []
            for r in rs:
                if not r["rmsd_to_agonist"] or not r["rmsd_to_antagonist"]:
                    continue
                to_ag, to_an = float(r["rmsd_to_agonist"]), float(r["rmsd_to_antagonist"])
                delta = to_an - to_ag
                calls.append("agonist" if delta > thresh else "antagonist")
            return majority(calls)

        bal, pooled, per_class = balanced_and_pooled(by_compound, call_at)
        sweep_rows.append({
            "delta_threshold": round(thresh, 4),
            "balanced_accuracy": round(bal, 4),
            "pooled_accuracy": round(pooled, 4),
            "agonist_accuracy": round(per_class["agonist"], 4),
            "antagonist_accuracy": round(per_class["antagonist"], 4),
        })
    best = max(sweep_rows, key=lambda r: r["balanced_accuracy"])
    print(f"optimum on this sweep: threshold={best['delta_threshold']}, "
          f"balanced_accuracy={best['balanced_accuracy']:.1%} "
          f"(current rule sits at threshold=0, balanced_accuracy={bal0:.1%})")
    with open(sweep_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sweep_rows[0].keys()))
        w.writeheader()
        w.writerows(sweep_rows)
    print(f"wrote {sweep_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
