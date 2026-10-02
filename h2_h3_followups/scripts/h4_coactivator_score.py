#!/usr/bin/env python3
"""H4 scoring: coactivator-peptide placement, D1 (agonist) vs D2 (inverse
agonist). Run this in the morning, after both arms have completed -- see
H4_PREREGISTRATION.md. Not run the night the arms were launched, on purpose.

For each structure:
  1. Superpose predicted chain A (receptor, CA atoms) onto 3KYT chain A
     (CA atoms, same residues only) via Kabsch -- fit on the RECEPTOR only,
     so the peptide's own placement can't inflate or deflate the fit.
  2. Apply that transform to predicted chain B (peptide, heavy atoms),
     compute RMSD to 3KYT chain C (peptide, 686-697 -> predicted 1-12,
     offset +685), same residue correspondence, heavy atoms only.
  3. Peptide-LBD interface contacts: predicted chain B heavy atoms within
     4.5 A of predicted chain A heavy atoms (same cutoff as
     per_residue_discriminator.py), reported as a contact-residue count.
  4. Official H12 state call on chain A, via the organizers' own
     lib.h12_state (5VB7/6T4I references, window 484-507) -- needs the
     rorgt_candidate_kit on disk and `apheris-data` wrapper on PATH, same
     as official_rescore.py.

    PYTHONPATH=/tmp/rorgt_kit/rorgt_candidate_kit PATH=$HOME/bin:$PATH \
        python3 h4_coactivator_score.py \
        --d1 ~/h2_posthoc/out/coactivator_d1 \
        --d2 ~/h2_posthoc/out/coactivator \
        --ref3kyt ~/data/apheris-foundry-steering-committee/apheris_kit_rorgamma/refs_rorgamma/3KYT.cif \
        --kit-root /tmp/rorgt_kit/rorgt_candidate_kit \
        --out h4_results.csv

Then run the pre-registered exact permutation test on D1 vs D2 seed-level
mean peptide RMSD -- same convention as H2/H3 (see per_residue_discriminator.py
for the test implementation, reusable here).
"""
from __future__ import annotations

import argparse
import csv
import glob
import re
import sys
from pathlib import Path

import gemmi
import numpy as np

PEPTIDE_OFFSET = 685  # predicted chain B residue 1 -> 3KYT chain C residue 686
RECEPTOR_OFFSET = 264  # predicted chain A residue 1 -> 3KYT/construct residue 265


def ca_coords(model, chain_id: str) -> dict[int, np.ndarray]:
    out = {}
    for res in model[chain_id]:
        for atom in res:
            if atom.name == "CA":
                out[res.seqid.num] = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
    return out


def heavy_coords(model, chain_id: str) -> dict[int, np.ndarray]:
    out: dict[int, list] = {}
    for res in model[chain_id]:
        pts = [
            [a.pos.x, a.pos.y, a.pos.z]
            for a in res
            if not a.is_hydrogen()
        ]
        if pts:
            out[res.seqid.num] = np.asarray(pts)
    return out


def kabsch(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Returns (R, t) such that mobile @ R + t ~= target, least squares."""
    mc, tc = mobile.mean(axis=0), target.mean(axis=0)
    mobile0, target0 = mobile - mc, target - tc
    u, _, vt = np.linalg.svd(mobile0.T @ target0)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = u @ np.diag([1, 1, d]) @ vt
    t = tc - mc @ r
    return r, t


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--d1", required=True, help="D1 (agonist) output dir")
    ap.add_argument("--d2", required=True, help="D2 (inverse agonist) output dir")
    ap.add_argument("--ref3kyt", required=True)
    ap.add_argument("--kit-root", required=True)
    ap.add_argument("--out", default="h4_results.csv")
    args = ap.parse_args()

    ref_st = gemmi.read_structure(args.ref3kyt)
    ref_model = ref_st[0]
    ref_ca_raw = ca_coords(ref_model, "A")  # keyed 265-507
    ref_peptide_heavy_raw = heavy_coords(ref_model, "C")  # keyed 686-697 (+ waters, filtered below)
    ref_peptide_heavy = {k: v for k, v in ref_peptide_heavy_raw.items() if 686 <= k <= 697}

    sys.path.insert(0, args.kit_root)
    from lib import h12_state  # noqa: E402

    rows = []
    for arm, root in (("D1_agonist", args.d1), ("D2_inverse_agonist", args.d2)):
        cifs = sorted(glob.glob(f"{root}/**/*_model.cif", recursive=True))
        for cif in cifs:
            seed = int(re.search(r"seed[_-](\d+)", cif).group(1))
            st = gemmi.read_structure(cif)
            model = st[0]
            pred_ca = ca_coords(model, "A")
            pred_peptide_heavy = heavy_coords(model, "B")
            pred_receptor_heavy = heavy_coords(model, "A")

            # map predicted chain A numbering (1-243) -> 3KYT numbering (265-507)
            common = sorted(
                (i, i + RECEPTOR_OFFSET) for i in pred_ca if (i + RECEPTOR_OFFSET) in ref_ca_raw
            )
            if len(common) < 20:
                rows.append({"arm": arm, "file": cif, "seed": seed, "status": "receptor_ca_mismatch"})
                continue
            mobile = np.array([pred_ca[i] for i, _ in common])
            target = np.array([ref_ca_raw[j] for _, j in common])
            r, t = kabsch(mobile, target)

            # map predicted chain B numbering (1-12) -> 3KYT chain C numbering (686-697)
            pep_common = sorted(
                (i, i + PEPTIDE_OFFSET) for i in pred_peptide_heavy if (i + PEPTIDE_OFFSET) in ref_peptide_heavy
            )
            if len(pep_common) < 8:
                rows.append({"arm": arm, "file": cif, "seed": seed, "status": "peptide_mismatch"})
                continue
            pred_pts, ref_pts = [], []
            for i, j in pep_common:
                p = pred_peptide_heavy[i] @ r + t
                n = min(len(p), len(ref_peptide_heavy[j]))
                pred_pts.append(p[:n])
                ref_pts.append(ref_peptide_heavy[j][:n])
            pred_pts = np.concatenate(pred_pts)
            ref_pts = np.concatenate(ref_pts)
            peptide_rmsd = float(np.sqrt(np.mean(np.sum((pred_pts - ref_pts) ** 2, axis=1))))

            # interface contacts: peptide heavy atoms within 4.5A of receptor heavy atoms
            receptor_pts = np.concatenate(list(pred_receptor_heavy.values())) if pred_receptor_heavy else np.zeros((0, 3))
            contact_residues = 0
            if len(receptor_pts):
                for pts in pred_peptide_heavy.values():
                    d = np.min(np.linalg.norm(pts[:, None, :] - receptor_pts[None, :, :], axis=-1))
                    if d <= 4.5:
                        contact_residues += 1

            rows.append({
                "arm": arm, "file": cif, "seed": seed, "status": "ok",
                "peptide_rmsd_to_3kyt": round(peptide_rmsd, 3),
                "peptide_contact_residues": contact_residues,
            })

    # official H12 state call, batched
    structures = {f"{r['arm']}__{Path(r['file']).stem}": Path(r["file"]) for r in rows if r.get("status") == "ok"}
    work = Path("h4_official_rescore_work")
    work.mkdir(exist_ok=True)
    calls = h12_state.call_states(structures, work, default_offset=h12_state.OFFSET_TO_UNIPROT)
    by_file = {str(path): calls.get(sid) for sid, path in structures.items()}
    for r in rows:
        c = by_file.get(r["file"])
        r["h12_state"] = c.state if c else ""
        r["h12_margin"] = c.margin if c else ""

    fieldnames = ["arm", "file", "seed", "status", "peptide_rmsd_to_3kyt",
                  "peptide_contact_residues", "h12_state", "h12_margin"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})
    print(f"wrote {args.out} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
