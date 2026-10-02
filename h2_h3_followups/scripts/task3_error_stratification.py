#!/usr/bin/env python3
"""Task 3 — error stratification on the existing 258-compound pool. No new
predictions, no new classifier, no touching h12_state_classifier.py.

Pool annotation does not carry an orthosteric/allosteric label (checked:
chembl_subset_120.json / chembl_subset_random150.json only have
molecule_chembl_id/smiles/activity_label/potency_nm/scaffold) -- per the
task's own fallback, stratify by scaffold class instead, plus mean pLDDT
and decisive margin.
"""
from __future__ import annotations

import csv
import json
import re
import statistics as stats
from collections import defaultdict
from pathlib import Path

H3 = Path.home() / "h3"
CALLS_CSV = H3 / "official_state_calls_all.csv"


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


def mean_plddt(cif_path: str) -> float | None:
    scores_path = Path(cif_path.replace("_model.cif", "_scores.json"))
    if not scores_path.exists():
        return None
    data = json.loads(scores_path.read_text())
    plddt = data.get("atom_predicted_local_distance_difference_test")
    if not plddt:
        return None
    return sum(plddt) / len(plddt)


def main() -> int:
    scaffold_by_cid: dict[str, str] = {}
    for fname in ("chembl_subset_120.json", "chembl_subset_random150.json"):
        for rec in json.loads((H3 / fname).read_text()):
            scaffold_by_cid[rec["molecule_chembl_id"]] = rec["scaffold"]

    rows = list(csv.DictReader(open(CALLS_CSV)))
    pool_rows = [r for r in rows if batch_of(r["file"])]
    by_compound: dict = defaultdict(list)
    for r in pool_rows:
        key = (batch_of(r["file"]), cid_of(r["file"]), true_label(r["file"]))
        by_compound[key].append(r)

    per_compound = []
    for (batch, cid, lbl), rs in by_compound.items():
        states = [r["state"] for r in rs]
        call = majority(states)
        margins = [float(r["margin"]) for r in rs if r["margin"]]
        plddts = [p for p in (mean_plddt(r["file"]) for r in rs) if p is not None]
        per_compound.append({
            "batch": batch, "cid": cid, "true_label": lbl, "call": call,
            "correct": call == lbl,
            "scaffold": scaffold_by_cid.get(cid, "?"),
            "mean_margin": sum(margins) / len(margins) if margins else None,
            "mean_plddt": sum(plddts) / len(plddts) if plddts else None,
        })

    out = Path(__file__).parent.parent / "results" / "task3_per_compound.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(per_compound[0].keys()))
        w.writeheader()
        w.writerows(per_compound)
    print(f"wrote {out} ({len(per_compound)} compounds)")

    n_errors = sum(1 for c in per_compound if not c["correct"])
    print(f"errors: {n_errors}/{len(per_compound)}")

    # --- stratify by scaffold ---
    print("\n=== by scaffold (scaffolds with >=3 compounds in the pool) ===")
    by_scaffold = defaultdict(list)
    for c in per_compound:
        by_scaffold[c["scaffold"]].append(c)
    scaffold_rows = []
    for sc, cs in by_scaffold.items():
        if len(cs) < 3:
            continue
        err = sum(1 for c in cs if not c["correct"])
        scaffold_rows.append((sc, len(cs), err, err / len(cs)))
    scaffold_rows.sort(key=lambda r: -r[3])
    for sc, n, err, rate in scaffold_rows:
        print(f"  n={n:3d}  errors={err:2d}  error_rate={rate:.1%}  scaffold={sc[:60]}")
    n_singleton_scaffolds = sum(1 for cs in by_scaffold.values() if len(cs) < 3)
    print(f"  ({n_singleton_scaffolds} scaffolds with <3 compounds, not shown individually)")

    # --- stratify by mean pLDDT quartile ---
    print("\n=== by mean pLDDT quartile ===")
    with_plddt = [c for c in per_compound if c["mean_plddt"] is not None]
    with_plddt.sort(key=lambda c: c["mean_plddt"])
    n = len(with_plddt)
    for q in range(4):
        lo, hi = n * q // 4, n * (q + 1) // 4
        chunk = with_plddt[lo:hi]
        err = sum(1 for c in chunk if not c["correct"])
        plddt_range = (chunk[0]["mean_plddt"], chunk[-1]["mean_plddt"])
        print(f"  Q{q+1} (pLDDT {plddt_range[0]:.1f}-{plddt_range[1]:.1f}): "
              f"n={len(chunk)} errors={err} error_rate={err/len(chunk):.1%}")
    if n < len(per_compound):
        print(f"  ({len(per_compound) - n} compounds missing a pLDDT: no scores.json found)")

    # --- stratify by decisive margin ---
    print("\n=== by mean decisive-margin band ===")
    bands = [(0, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 100)]
    with_margin = [c for c in per_compound if c["mean_margin"] is not None]
    for lo, hi in bands:
        chunk = [c for c in with_margin if lo <= c["mean_margin"] < hi]
        if not chunk:
            continue
        err = sum(1 for c in chunk if not c["correct"])
        print(f"  margin [{lo},{hi}): n={len(chunk)} errors={err} error_rate={err/len(chunk):.1%}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
