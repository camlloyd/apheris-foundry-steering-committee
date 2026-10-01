#!/usr/bin/env python3
"""The rubric's accuracy metrics, computed locally.

    python3 local_metrics.py --refs <kit>/refs_rorgamma --runs <dir> \
        --ref-state H12-out --out metrics.csv

The challenge rubric scores prediction accuracy (60 %) as the H12-domain RMSD
from the reference plus the mean of **GDT-HA**, **LDDT-PLI** and **ligand
BiSyRMSD**. This computes all four.

Use `run_metrics.sh` first — if Foundry's own `apheris-data compute-metrics`
produces these, prefer its numbers and cite the platform. This script exists so
the metrics exist either way, and as a cross-check on whichever is used.

Definitions implemented (stated because small choices change the number):

* **GDT-HA** — superimpose on CA atoms, then mean over the four high-accuracy
  cutoffs {0.5, 1, 2, 4 Å} of the fraction of CA within that distance.
  Superposition is refined iteratively from the best-fitting subset rather than
  a single global fit, which is what the GDT family specifies.
* **lDDT** — the standard superposition-free local score: for every pair of
  atoms within a 15 Å inclusion radius in the reference, the fraction of
  distances preserved within {0.5, 1, 2, 4 Å}, averaged.
* **LDDT-PLI** — the same, restricted to protein-atom-to-ligand-atom pairs, so
  it measures the interface rather than the fold.
* **Ligand RMSD** — symmetry-corrected via RDKit substructure matching when
  RDKit is present (the "BiSy" part: the best RMSD over all automorphisms of
  the ligand graph). Falls back to a plain atom-order RMSD, and the column
  `ligand_rmsd_symmetry_corrected` records which was used. **Note:** full
  BiSyRMSD additionally symmetrises over receptor chain assignments; with a
  single receptor chain that reduces to the ligand-symmetry case implemented
  here. Say so rather than claiming the published metric exactly.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import pathlib
import sys

import numpy as np

GDT_HA_CUTOFFS = (0.5, 1.0, 2.0, 4.0)
LDDT_THRESHOLDS = (0.5, 1.0, 2.0, 4.0)
LDDT_INCLUSION = 15.0

try:
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    HAVE_RDKIT = True
except ImportError:
    HAVE_RDKIT = False


# ----------------------------------------------------------------- GDT-HA
def gdt_ha(pred_ca: np.ndarray, ref_ca: np.ndarray, superpose_on, rmsd) -> float:
    """Mean fraction of CA within {0.5,1,2,4} A after iterative superposition."""
    fracs = []
    for cut in GDT_HA_CUTOFFS:
        mask = np.ones(len(ref_ca), dtype=bool)
        best = 0.0
        for _ in range(20):
            if mask.sum() < 3:
                break
            fitted = superpose_on(pred_ca, ref_ca, mask)
            d = np.linalg.norm(fitted - ref_ca, axis=1)
            new_mask = d <= cut
            best = max(best, new_mask.mean())
            if new_mask.sum() < 3 or np.array_equal(new_mask, mask):
                break
            mask = new_mask
        fracs.append(best)
    return 100.0 * float(np.mean(fracs))


# ------------------------------------------------------------------- lDDT
def _lddt_from_pairs(ref_d: np.ndarray, pred_d: np.ndarray) -> float:
    if ref_d.size == 0:
        return float("nan")
    diff = np.abs(ref_d - pred_d)
    return 100.0 * float(np.mean([(diff < t).mean() for t in LDDT_THRESHOLDS]))


def lddt_backbone(pred_ca: np.ndarray, ref_ca: np.ndarray) -> float:
    """Superposition-free lDDT over CA atoms."""
    rd = np.linalg.norm(ref_ca[:, None, :] - ref_ca[None, :, :], axis=-1)
    pd = np.linalg.norm(pred_ca[:, None, :] - pred_ca[None, :, :], axis=-1)
    iu = np.triu_indices(len(ref_ca), k=1)
    sel = rd[iu] < LDDT_INCLUSION
    return _lddt_from_pairs(rd[iu][sel], pd[iu][sel])


def lddt_pli(pred_prot: np.ndarray, pred_lig: np.ndarray,
             ref_prot: np.ndarray, ref_lig: np.ndarray) -> float:
    """lDDT restricted to protein-ligand atom pairs — the interface score."""
    if len(pred_lig) == 0 or len(ref_lig) == 0:
        return float("nan")
    n = min(len(pred_lig), len(ref_lig))
    rd = np.linalg.norm(ref_prot[:, None, :] - ref_lig[None, :n, :], axis=-1)
    pd = np.linalg.norm(pred_prot[:, None, :] - pred_lig[None, :n, :], axis=-1)
    sel = rd < LDDT_INCLUSION
    return _lddt_from_pairs(rd[sel], pd[sel])


# ----------------------------------------------------------- ligand RMSD
def _automorphisms(mol):
    try:
        matches = mol.GetSubstructMatches(mol, uniquify=False, useChirality=False,
                                          maxMatches=20000)
        return list(matches)
    except Exception:
        return []


def ligand_rmsd(pred_xyz: np.ndarray, ref_xyz: np.ndarray, smiles: str | None):
    """Symmetry-corrected ligand RMSD when RDKit can build the graph."""
    n = min(len(pred_xyz), len(ref_xyz))
    if n == 0:
        return float("nan"), False
    p, r = pred_xyz[:n], ref_xyz[:n]
    plain = float(np.sqrt(((p - r) ** 2).sum(1).mean()))
    if not (HAVE_RDKIT and smiles):
        return plain, False
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() != n:
        return plain, False
    best = plain
    for perm in _automorphisms(mol):
        if len(perm) != n:
            continue
        v = float(np.sqrt(((p[list(perm)] - r) ** 2).sum(1).mean()))
        best = min(best, v)
    return best, True


# ------------------------------------------------------------------ main
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refs", required=True, help="refs_rorgamma directory")
    ap.add_argument("--runs", required=True, help="directory of predictions (searched recursively)")
    ap.add_argument("--ref-state", choices=["H12-in", "H12-out"], default="H12-out",
                    help="which crystal to score against (the direction's target state)")
    ap.add_argument("--rgkit", default=None)
    ap.add_argument("--smiles", default=None,
                    help="ligand SMILES, for symmetry-corrected RMSD")
    ap.add_argument("--out", default="metrics.csv")
    args = ap.parse_args()

    for p in filter(None, [args.rgkit, str(pathlib.Path(__file__).resolve().parent)]):
        sys.path.insert(0, p)
    from rgkit import pocket_sets as P
    from rgkit import structio as S

    refs = pathlib.Path(args.refs)
    ref_file = "3KYT_clean.cif" if args.ref_state == "H12-in" else "4ZJW_clean.cif"
    ref = S.load_chain(refs / ref_file, source=args.ref_state)

    try:
        rlp = P.pick_ligand(ref)
        ref_lig = np.array([[a.pos.x, a.pos.y, a.pos.z]
                            for ch in ref.structure[0] if ch.name == rlp.chain_id
                            for res in ch if res.name == rlp.name
                            for a in res if a.element.name != "H"])
    except Exception:
        rlp, ref_lig = None, np.zeros((0, 3))

    rows = []
    cifs = sorted(p for p in pathlib.Path(args.runs).rglob("*.cif")
                  if "__MACOSX" not in str(p))
    print(f"scoring {len(cifs)} predictions against {ref_file} "
          f"({'RDKit present' if HAVE_RDKIT else 'no RDKit — plain ligand RMSD'})\n")

    for cif in cifs:
        row = {"file": str(cif.relative_to(args.runs)), "ref_state": args.ref_state}
        try:
            pred = S.load_chain(cif, source=cif.name)
            pairs = S.align_pair(ref, pred)
            ia = np.array([p[0] for p in pairs])
            ib = np.array([p[1] for p in pairs])
            rca, pca = ref.ca[ia], pred.ca[ib]

            row["n_aligned"] = len(pairs)
            row["gdt_ha"] = round(gdt_ha(pca, rca, S.superpose_on, S.rmsd), 2)
            row["lddt"] = round(lddt_backbone(pca, rca), 2)

            try:
                plp = P.pick_ligand(pred)
                pred_lig = np.array([[a.pos.x, a.pos.y, a.pos.z]
                                     for ch in pred.structure[0] if ch.name == plp.chain_id
                                     for res in ch if res.name == plp.name
                                     for a in res if a.element.name != "H"])
            except Exception:
                pred_lig = np.zeros((0, 3))

            if len(pred_lig) and len(ref_lig):
                # Fit on protein CA only, then carry the ligand through the same
                # transform. The padding rows are placeholders excluded by the mask.
                mob = np.vstack([pca, pred_lig])
                pad = np.repeat(rca.mean(0)[None, :], len(pred_lig), axis=0)
                tgt = np.vstack([rca, pad])
                fit_mask = np.zeros(len(mob), dtype=bool)
                fit_mask[:len(pca)] = True
                plig_fit = S.superpose_on(mob, tgt, fit_mask)[len(pca):]
                row["lddt_pli"] = round(lddt_pli(pca, pred_lig, rca, ref_lig), 2)
                lr, sym = ligand_rmsd(plig_fit, ref_lig, args.smiles)
                row["ligand_rmsd"] = round(lr, 3)
                row["ligand_rmsd_symmetry_corrected"] = int(sym)
            else:
                row["lddt_pli"] = ""
                row["ligand_rmsd"] = ""
                row["ligand_rmsd_symmetry_corrected"] = 0

            vals = [row["gdt_ha"], row["lddt_pli"], row["ligand_rmsd"]]
            row["status"] = "ok"
            print(f"  {row['file'][:52]:<54} GDT-HA {row['gdt_ha']:6.2f}  "
                  f"lDDT {row['lddt']:6.2f}  LDDT-PLI {str(row['lddt_pli']):>6}  "
                  f"ligRMSD {str(row['ligand_rmsd']):>7}")
        except Exception as exc:
            row["status"] = f"error: {exc}"[:90]
            print(f"  {row['file'][:52]:<54} ERROR {exc}")
        rows.append(row)

    fields = ["file", "ref_state", "status", "n_aligned", "gdt_ha", "lddt",
              "lddt_pli", "ligand_rmsd", "ligand_rmsd_symmetry_corrected"]
    out = pathlib.Path(args.out)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    sym = [r for r in rows if r.get("ligand_rmsd_symmetry_corrected")]
    if sym and HAVE_RDKIT and args.smiles:
        from rdkit import Chem as _C
        _m = _C.MolFromSmiles(args.smiles)
        if _m is not None:
            n_auto = len(_m.GetSubstructMatches(_m, uniquify=False, useChirality=False))
            if n_auto <= 1:
                print(f"\nnote: this ligand has {n_auto} graph automorphism, so the "
                      "symmetry correction is a no-op here — the value is a plain\n"
                      "      ligand RMSD after protein superposition. Say that rather "
                      "than claiming BiSyRMSD.")
    ok = [r for r in rows if r["status"] == "ok"]
    if ok:
        print(f"\nmean over {len(ok)} scored:  GDT-HA "
              f"{np.mean([r['gdt_ha'] for r in ok]):.2f}   "
              f"lDDT {np.mean([r['lddt'] for r in ok]):.2f}")
    print(f"-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
