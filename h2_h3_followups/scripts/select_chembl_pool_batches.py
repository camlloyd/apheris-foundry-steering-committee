#!/usr/bin/env python3
"""Reproduces the two ChEMBL pool batches scored in OFFICIAL_RESCORE_FINDINGS.md.

Batch A (120, diversity/potency-curated) is deterministic: most potent
compound per unique Bemis-Murcko scaffold, 60 agonist + 60 antagonist.

Batch B (150, true random draw, excluding batch A) uses Python's `random`
with a fixed seed (20261001, chosen as the run date) so the draw is
reproducible -- this is the number that was missing from the original run
(the selection was done ad hoc in a shell heredoc and not persisted as a
script). Re-running this file regenerates byte-identical
chembl_subset_120.json / chembl_subset_random150.json.

    python3 select_chembl_pool_batches.py \
        --chembl-csv /tmp/rorgt_kit/rorgt_candidate_kit/dataset/train/chembl/chembl_compounds.csv \
        --out-dir ~/h3
"""
from __future__ import annotations

import argparse
import csv
import json
import random
from pathlib import Path

RANDOM_SEED = 20261001  # the run date -- fixed so batch B is reconstructible


def pick_diverse(rows: list[dict], n: int) -> list[dict]:
    rows_sorted = sorted(rows, key=lambda r: r["potency_nm"])
    seen_scaffold: set[str] = set()
    out = []
    for r in rows_sorted:
        if r["scaffold"] in seen_scaffold:
            continue
        seen_scaffold.add(r["scaffold"])
        out.append(r)
        if len(out) >= n:
            break
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chembl-csv", required=True)
    ap.add_argument("--out-dir", default=str(Path.home() / "h3"))
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.chembl_csv)))
    for r in rows:
        r["potency_nm"] = float(r["potency_nm"])

    by_label: dict[str, list[dict]] = {"agonist": [], "antagonist": []}
    for r in rows:
        by_label[r["activity_label"]].append(r)

    batch_a = pick_diverse(by_label["agonist"], 60) + pick_diverse(by_label["antagonist"], 60)
    used = {r["molecule_chembl_id"] for r in batch_a}

    remaining = {"agonist": [], "antagonist": []}
    for r in rows:
        if r["molecule_chembl_id"] not in used:
            remaining[r["activity_label"]].append(r)

    random.seed(RANDOM_SEED)
    batch_b = random.sample(remaining["agonist"], 75) + random.sample(remaining["antagonist"], 75)
    random.shuffle(batch_b)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "chembl_subset_120.json").write_text(json.dumps(batch_a, indent=2))
    (out_dir / "chembl_subset_random150.json").write_text(json.dumps(batch_b, indent=2))
    print(f"batch A (curated): {len(batch_a)} -> {out_dir / 'chembl_subset_120.json'}")
    print(f"batch B (random, seed={RANDOM_SEED}): {len(batch_b)} -> {out_dir / 'chembl_subset_random150.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
