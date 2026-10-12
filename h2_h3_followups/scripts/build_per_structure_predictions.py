#!/usr/bin/env python3
"""Build results/per_structure_predictions.csv: one row per scored
structure, official frame only (5VB7/6T4I, window 484-507). No submission
schema exists anywhere in rorgt_candidate_kit (checked: README.md,
HACKATHON.md, dataset/ -- README.md says outright "the set submissions
are scored on is not in [this kit]"), so this is our own format.

Sources, no new GPU jobs:
  - ~/h3/official_state_calls_all.csv  -- rmsd_to_agonist/antagonist,
    margin, state, for every H1/H2/H2-posthoc/H3-core/H3-extended/pool
    structure already scored against the official reference (1,881 rows).
  - h12_state.call_states(), run fresh here, for H4's D1 arm only (50
    structures) -- the one round that postdates the official rescore.
  - Each structure's own *_scores.json for mean_plddt (whole-structure)
    and h12_local_plddt (official window 484-507, tokens 219-242 0-based
    -- see note below on the off-by-one bug this caught in the existing
    H3 classifier's own window math).
  - chembl_subset_120.json / chembl_subset_random150.json for pool
    true_label and the 12 compounds attempted but never scored.
"""
from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

H3_RAW = Path.home() / "h3"
H2_RAW = Path.home() / "h2"
H2P_RAW = Path.home() / "h2_posthoc"
APHERIS_RUNS = Path("/home/lyceum/data/apheris-foundry-steering-committee/apheris_kit_rorgamma/runs")
OUT = Path(__file__).parent.parent / "results" / "per_structure_predictions.csv"

OFFICIAL_WINDOW = (484, 507)          # UniProt/3KYT numbering, organizers' own window
CONSTRUCT_OFFSET = 264                # predicted residue i (1-based) -> construct residue i+264
# token (0-based, from atom_tokens) for construct residue r -> r - CONSTRUCT_OFFSET - 1
TOKEN_WINDOW = (OFFICIAL_WINDOW[0] - CONSTRUCT_OFFSET - 1, OFFICIAL_WINDOW[1] - CONSTRUCT_OFFSET - 1)


def h12_local_plddt(scores_path: Path) -> float | None:
    """Mean pLDDT over the official H12 window (484-507), receptor chain
    only. Chain A occupies the first len(chain A) tokens regardless of
    how many other chains follow (verified directly against gemmi atom
    order for both single-chain and H4's 3-chain structures), so the same
    fixed token range applies project-wide without a per-structure chain
    lookup.
    """
    try:
        d = json.loads(scores_path.read_text())
    except Exception:
        return None
    atoms = d.get("atom_predicted_local_distance_difference_test")
    toks = d.get("atom_tokens")
    if not atoms or not toks:
        return None
    lo, hi = TOKEN_WINDOW
    vals = [a for a, t in zip(atoms, toks) if lo <= t <= hi]
    return sum(vals) / len(vals) if vals else None


def mean_plddt(scores_path: Path) -> float | None:
    try:
        d = json.loads(scores_path.read_text())
    except Exception:
        return None
    atoms = d.get("atom_predicted_local_distance_difference_test")
    return sum(atoms) / len(atoms) if atoms else None


def scores_path_for(model_cif: str) -> Path:
    return Path(model_cif.replace("_model.cif", "_scores.json"))


ARM_TRUE_LABEL = {
    "apo": "", "c1": "agonist", "c2": "antagonist", "c3": "agonist",
    "c4": "antagonist", "c5": "antagonist", "c6": "antagonist",
}


def classify(path: str):
    """Returns (round, arm, compound_id, seed, sample, true_label) or None."""
    if "/apheris_kit_rorgamma/runs/" in path:
        m = re.search(r"/runs/(rorgamma_invago|rorgamma_25hc)/([^/]+)/\2_seed_(\d+)_sample_(\d+)_model\.cif$", path)
        if not m:
            return None
        target, arm, seed, sample = m.groups()
        true_label = "antagonist" if target == "rorgamma_invago" else "agonist"
        return "h1", arm, arm, seed, sample, true_label

    if "/h2/out/" in path:
        m = re.search(r"/out/(depth_\d+)/[^/]+/seed_(\d+)/[^/]+_sample_(\d+)_model\.cif$", path)
        if not m:
            return None
        arm, seed, sample = m.groups()
        return "h2", arm, arm, seed, sample, "antagonist"

    if "/h2_posthoc/out/" in path:
        m = re.search(r"/out/([^/]+)/[^/]+/seed_(\d+)/[^/]+_sample_(\d+)_model\.cif$", path)
        if not m:
            return None
        arm, seed, sample = m.groups()
        if arm == "coactivator_d1":
            return "h4", "D1_agonist", "D1_agonist", seed, sample, "agonist"
        # every remaining h2_posthoc arm (depth_48/64/96, template_3kyt/4zjw/5c4o,
        # coactivator/D2) uses 4P1 -- confirmed directly in run_msa_titration.sh's
        # hardcoded query and run_bisect.sh/run_template.sh reusing it unchanged.
        return "h2_posthoc", arm, arm, seed, sample, "antagonist"

    if "/h3/out/" in path:
        m = re.search(r"/out/(pool120|rand150)/((?:POOL|RAND)_([A-Za-z0-9]+)_(agonist|antagonist))/seed_(\d+)/[^/]+_sample_(\d+)_model\.cif$", path)
        if m:
            batch, _qid, cid, label, seed, sample = m.groups()
            arm = "batch_A_curated_120" if batch == "pool120" else "batch_B_random_150"
            return "pool", arm, cid, seed, sample, label
        m = re.search(r"/out/(apo|c1|c2|c3|c4|c5|c6)/[^/]+/seed_(\d+)/[^/]+_sample_(\d+)_model\.cif$", path)
        if m:
            arm, seed, sample = m.groups()
            return "h3", arm, arm, seed, sample, ARM_TRUE_LABEL[arm]
        return None

    return None


def main() -> int:
    rows_out: list[dict] = []

    # --- everything already in the official rescore ---
    official = list(csv.DictReader(open(H3_RAW / "official_state_calls_all.csv")))
    unclassified = 0
    for r in official:
        cls = classify(r["file"])
        if cls is None:
            unclassified += 1
            continue
        round_, arm, compound_id, seed, sample, true_label = cls
        scores_path = scores_path_for(r["file"])
        mp = mean_plddt(scores_path)
        h12p = h12_local_plddt(scores_path)
        to_ag = float(r["rmsd_to_agonist"]) if r["rmsd_to_agonist"] else None
        to_an = float(r["rmsd_to_antagonist"]) if r["rmsd_to_antagonist"] else None
        signed_margin = (to_ag - to_an) if (to_ag is not None and to_an is not None) else None
        if r["state"] == "unscored" or to_ag is None:
            call = "unscorable"
        else:
            call = r["state"]
        rows_out.append({
            "compound_id": compound_id, "round": round_, "arm": arm,
            "seed": seed, "sample_idx": sample,
            "mean_plddt": round(mp, 2) if mp is not None else "",
            "h12_local_plddt": round(h12p, 2) if h12p is not None else "",
            "rmsd_to_agonist_ref": round(to_ag, 3) if to_ag is not None else "",
            "rmsd_to_antagonist_ref": round(to_an, 3) if to_an is not None else "",
            "margin": round(signed_margin, 3) if signed_margin is not None else "",
            "quality_floor_pass": (mp >= 70) if mp is not None else "",
            "call": call,
            "true_label": true_label,
        })
    print(f"from official_state_calls_all.csv: {len(rows_out)} rows, {unclassified} unclassified (investigate if > 0)", file=sys.stderr)

    # --- H4 D1 (new, postdates the official rescore -- score it now) ---
    sys.path.insert(0, "/tmp/rorgt_kit/rorgt_candidate_kit")
    from lib import h12_state
    import glob
    d1_cifs = sorted(glob.glob(str(H2P_RAW / "out/coactivator_d1/**/*_model.cif"), recursive=True))
    structures = {}
    for f in d1_cifs:
        m = re.search(r"seed_(\d+)_sample_(\d+)_model\.cif$", f)
        structures[f"D1_{m.group(1)}_{m.group(2)}"] = Path(f)
    work = Path("/tmp/d1_official_rescore_work")
    work.mkdir(exist_ok=True)
    calls = h12_state.call_states(structures, work, default_offset=h12_state.OFFSET_TO_UNIPROT)
    for sid, path in structures.items():
        _, seed, sample = sid.split("_")
        c = calls[sid]
        scores_path = scores_path_for(str(path))
        mp = mean_plddt(scores_path)
        h12p = h12_local_plddt(scores_path)
        to_ag, to_an = c.rmsd_to_agonist, c.rmsd_to_antagonist
        signed_margin = (to_ag - to_an) if (to_ag is not None and to_an is not None) else None
        call = "unscorable" if c.state == h12_state.UNSCORED else c.state
        rows_out.append({
            "compound_id": "D1_agonist", "round": "h4", "arm": "D1_agonist",
            "seed": seed, "sample_idx": sample,
            "mean_plddt": round(mp, 2) if mp is not None else "",
            "h12_local_plddt": round(h12p, 2) if h12p is not None else "",
            "rmsd_to_agonist_ref": round(to_ag, 3) if to_ag is not None else "",
            "rmsd_to_antagonist_ref": round(to_an, 3) if to_an is not None else "",
            "margin": round(signed_margin, 3) if signed_margin is not None else "",
            "quality_floor_pass": (mp >= 70) if mp is not None else "",
            "call": call,
            "true_label": "agonist",
        })
    print(f"added {len(structures)} H4/D1 rows (freshly scored, official frame)", file=sys.stderr)

    # --- the 12 rand150 compounds attempted but never scored (abstentions, not omitted) ---
    labels_by_cid = {}
    for fname in ("chembl_subset_120.json", "chembl_subset_random150.json"):
        for rec in json.loads((H3_RAW / fname).read_text()):
            labels_by_cid[rec["molecule_chembl_id"]] = rec["activity_label"]
    scored_rand150_cids = {r["compound_id"] for r in rows_out if r["round"] == "pool" and r["arm"] == "batch_B_random_150"}
    all_rand150_cids = {rec["molecule_chembl_id"] for rec in json.loads((H3_RAW / "chembl_subset_random150.json").read_text())}
    missing = sorted(all_rand150_cids - scored_rand150_cids)
    for cid in missing:
        rows_out.append({
            "compound_id": cid, "round": "pool", "arm": "batch_B_random_150",
            "seed": "", "sample_idx": "",
            "mean_plddt": "", "h12_local_plddt": "",
            "rmsd_to_agonist_ref": "", "rmsd_to_antagonist_ref": "", "margin": "",
            "quality_floor_pass": "",
            "call": "no_call",
            "true_label": labels_by_cid.get(cid, ""),
        })
    print(f"added {len(missing)} no_call rows for rand150 compounds never scored "
          f"(documented RDKit featurizer failure): {missing}", file=sys.stderr)

    fieldnames = ["compound_id", "round", "arm", "seed", "sample_idx", "mean_plddt",
                  "h12_local_plddt", "rmsd_to_agonist_ref", "rmsd_to_antagonist_ref",
                  "margin", "quality_floor_pass", "call", "true_label"]
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows_out:
            w.writerow(r)
    print(f"\nwrote {OUT} ({len(rows_out)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
