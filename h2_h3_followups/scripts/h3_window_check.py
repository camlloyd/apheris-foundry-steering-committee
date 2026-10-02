#!/usr/bin/env python3
"""H12 window sensitivity check: 479-486 (auto-derived, primary) vs 501-507
(the challenge deck's window, through Tyr502/Phe506).

4ZJW does not resolve 502-507 (see INTEGRITY_NOTES.md), so no RMSD-to-4ZJW
state call is possible on that window -- but per-residue pLDDT only needs
the predicted structure, so the same seed-level permutation test from
H3_RESULTS_REPORT.md section 3.1 can be re-run on it directly.

    python3 h3_window_check.py --work ~/h3 --arms apo c1 c2 c3
"""
from __future__ import annotations

import argparse
import glob
import itertools
from collections import defaultdict

import numpy as np

from per_residue_discriminator import per_residue_plddt, seed_of, permutation_test, OFFSET

WINDOWS = {
    "H12_479-486_primary": set(range(479 - OFFSET, 487 - OFFSET)),
    "H12_501-507_deck": set(range(501 - OFFSET, 508 - OFFSET)),
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--work", required=True)
    ap.add_argument("--arms", nargs="+", required=True)
    args = ap.parse_args()

    data = {arm: [] for arm in args.arms}
    for arm in args.arms:
        for cif in sorted(glob.glob(f"{args.work}/out/{arm}/**/*.cif", recursive=True)):
            scores = cif.replace("_model.cif", "_scores.json")
            data[arm].append({"seed": seed_of(cif), "pr": per_residue_plddt(scores)})

    def seed_means(arm, residue_set):
        by_seed = defaultdict(list)
        for row in data[arm]:
            vals = [row["pr"][t] for t in residue_set if t in row["pr"]]
            if vals:
                by_seed[row["seed"]].append(np.mean(vals))
        seeds = sorted(by_seed)
        return np.array([np.mean(by_seed[s]) for s in seeds])

    for wname, wset in WINDOWS.items():
        print(f"\n=== {wname} ===")
        means = {}
        for arm in args.arms:
            m = seed_means(arm, wset)
            means[arm] = m
            print(f"  {arm}: {np.round(m, 2)}  mean={np.mean(m):.2f}")
        print("  -- pairwise --")
        for a, b in itertools.combinations(args.arms, 2):
            obs, p, d, cl = permutation_test(means[a], means[b])
            print(f"  {a} vs {b}: diff={obs:+.2f} p={p:.4f} d={d:+.2f} cliffs={cl:+.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
