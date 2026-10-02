"""
Self-test for the Apheris multi-state kit.

Builds synthetic structures with KNOWN geometry and checks that every metric
returns what it should. Run this first on Day 1 -- if it passes, the kit is
working and any weird number you see later is about the data or the model,
not about the code.

    python test_kit.py
"""
from __future__ import annotations

import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from structure_io import Residue, Structure, load_structure, write_pdb
from kinase_state import (
    classify_dfg, find_kinase_motifs, transfer_motifs,
)
from state_recovery import (
    apply_transform, contact_recovery, kabsch, ligand_centroid_distance,
    pocket_contacts, residue_correspondence, state_recovery, tm_score,
)

PASS, FAIL = "  ok  ", " FAIL "
_failures = []


def check(name, condition, detail=""):
    print(f"[{PASS if condition else FAIL}] {name}" + (f"  ({detail})" if detail else ""))
    if not condition:
        _failures.append(name)


# --------------------------------------------------------- synthetic kinase

def build_kinase(d1: float, d2: float, seed: int = 0, start_resid: int = 1,
                 with_ligand: bool = False, ligand_at=None) -> Structure:
    """
    A fake kinase domain with the motifs the classifier looks for, placed so
    that D1 and D2 come out at the values you ask for.

    Layout: 200 residues, VAIK Lys at index 70, aC Glu at index 90,
    HRD at 140, DFG at 160.
    """
    rng = np.random.default_rng(seed)
    n = 200
    seq = ["ALA"] * n
    seq[70:74] = ["VAL", "ALA", "ILE", "LYS"]      # VAIK -> b3-Lys at index 73
    seq[90] = "GLU"                                 # aC-Glu
    seq[94] = "LEU"                                 # aC-Glu+4, the D1 anchor
    seq[140:143] = ["HIS", "ARG", "ASP"]            # HRD catalytic loop
    seq[160:163] = ["ASP", "PHE", "GLY"]            # DFG

    st = Structure(path="<synthetic>")

    # A gentle helix as the backbone so the core has real 3D structure.
    for i in range(n):
        t = i * 0.35
        ca = np.array([12.0 * np.cos(t), 12.0 * np.sin(t), 1.5 * i]) \
            + rng.normal(0, 0.05, 3)
        st.residues.append(
            Residue(seqid=start_resid + i, name=seq[i], chain="A", atoms={"CA": ca})
        )

    b3_lys = st.residues[73]
    ac_glu = st.residues[90]
    ac_glu_p4 = st.residues[94]
    dfg_phe = st.residues[161]

    # Pin the three measured atoms to give the requested D1/D2 exactly.
    phe_cz = np.array([0.0, 0.0, 0.0])
    b3_lys.atoms["CA"] = np.array([d2, 0.0, 0.0])           # D2 = |Lys CA - Phe CZ|
    ac_glu_p4.atoms["CA"] = np.array([0.0, d1, 0.0])        # D1 = |Glu+4 CA - Phe CZ|
    dfg_phe.atoms["CZ"] = phe_cz
    dfg_phe.atoms["CA"] = phe_cz + np.array([0.0, 0.0, 3.0])

    # Salt bridge: Lys NZ next to Glu OE1, so autodetection finds this pair.
    b3_lys.atoms["NZ"] = b3_lys.atoms["CA"] + np.array([1.5, 2.0, 0.0])
    ac_glu.atoms["CA"] = b3_lys.atoms["NZ"] + np.array([4.0, 1.0, 0.0])
    ac_glu.atoms["OE1"] = b3_lys.atoms["NZ"] + np.array([2.6, 0.3, 0.0])
    ac_glu.atoms["CD"] = ac_glu.atoms["OE1"] + np.array([0.6, 0.4, 0.0])

    if with_ligand:
        centre = np.array(ligand_at if ligand_at is not None else [4.0, 2.0, 0.0])
        lig = Residue(seqid=900, name="LIG", chain="A", atoms={}, is_polymer=False)
        for k in range(12):
            lig.atoms[f"C{k+1}"] = centre + rng.normal(0, 1.2, 3)
        st.hetero.append(lig)

    return st


# ------------------------------------------------------------------- tests

def test_superposition():
    print("\n--- superposition ---")
    rng = np.random.default_rng(1)
    P = rng.normal(0, 10, (50, 3))

    R_true, _ = np.linalg.qr(rng.normal(0, 1, (3, 3)))
    if np.linalg.det(R_true) < 0:
        R_true[:, 0] *= -1
    t_true = np.array([5.0, -3.0, 11.0])
    Q = P @ R_true + t_true

    R, t = kabsch(Q, P)
    rmsd = np.sqrt(((apply_transform(Q, R, t) - P) ** 2).sum(axis=1).mean())
    check("Kabsch recovers a known rigid transform", rmsd < 1e-8, f"rmsd={rmsd:.2e}")

    check("TM-score of identical coords is 1.0",
          abs(tm_score(P, P) - 1.0) < 1e-9, f"tm={tm_score(P, P):.6f}")

    shifted = P + np.array([20.0, 0.0, 0.0])
    check("TM-score of badly wrong coords is low",
          tm_score(shifted, P) < 0.1, f"tm={tm_score(shifted, P):.4f}")


def test_motif_detection():
    print("\n--- motif detection ---")
    st = build_kinase(d1=8.0, d2=16.0)
    m = find_kinase_motifs(st)
    print("    detected:", m.describe(st))

    check("finds the DFG-Asp", m.dfg_asp == 161, f"got {m.dfg_asp}, want 161")
    check("finds the DFG-Phe", m.dfg_phe == 162, f"got {m.dfg_phe}, want 162")
    check("finds the b3-Lys", m.b3_lys == 74, f"got {m.b3_lys}, want 74")
    check("finds the aC-Glu", m.ac_glu == 91, f"got {m.ac_glu}, want 91")
    check("finds aC-Glu+4", m.ac_glu_p4 == 95, f"got {m.ac_glu_p4}, want 95")
    check("finds the HRD-His", m.hrd_his == 141, f"got {m.hrd_his}, want 141")


def test_dfg_classification():
    print("\n--- DFG classification (KinCore D1/D2) ---")
    cases = [
        ("DFG-in",    8.0, 16.0, "D1<11, D2>14"),
        ("DFG-out",  13.0, 10.0, "D1>11, D2<14"),
        ("DFG-inter", 9.0,  9.0, "D1<11, D2<11"),
        ("unassigned", 13.0, 16.0, "D1>11, D2>14 -> outside all boxes"),
    ]
    for want, d1, d2, why in cases:
        st = build_kinase(d1=d1, d2=d2)
        res = classify_dfg(st)
        ok = res.label == want
        check(f"D1={d1:4.1f} D2={d2:4.1f} -> {want}", ok,
              f"got '{res.label}'; {why}")
        if ok and res.d1 is not None:
            check(f"    D1/D2 measured correctly for {want}",
                  abs(res.d1 - d1) < 1e-6 and abs(res.d2 - d2) < 1e-6,
                  f"D1={res.d1:.3f} D2={res.d2:.3f}")


def test_motif_transfer_across_numbering():
    print("\n--- motif transfer across different numbering ---")
    ref = build_kinase(d1=13.0, d2=10.0, start_resid=242)   # crystal-like numbering
    pred = build_kinase(d1=13.0, d2=10.0, start_resid=1)    # prediction numbering

    ref_motifs = find_kinase_motifs(ref)
    moved = transfer_motifs(ref, ref_motifs, pred)
    print("    reference:", ref_motifs.describe(ref))
    print("    transferred:", moved.describe(pred))

    check("DFG-Phe transfers with the right offset",
          moved.dfg_phe == ref_motifs.dfg_phe - 241,
          f"{ref_motifs.dfg_phe} -> {moved.dfg_phe}")

    res = classify_dfg(pred, moved)
    check("transferred motifs still classify correctly",
          res.label == "DFG-out", f"got '{res.label}'")


def test_state_recovery():
    print("\n--- dual-reference state recovery ---")
    ref_in = build_kinase(d1=8.0, d2=16.0, seed=0)
    ref_out = build_kinase(d1=13.0, d2=10.0, seed=0)

    # The region that actually differs between the two references.
    region = {160, 161, 162, 163, 74, 95}

    # A "prediction" that is the DFG-out reference with a little noise and a
    # random rigid transform applied -- exactly what a model output looks like.
    rng = np.random.default_rng(7)
    pred = build_kinase(d1=13.0, d2=10.0, seed=0)
    R_true, _ = np.linalg.qr(rng.normal(0, 1, (3, 3)))
    if np.linalg.det(R_true) < 0:
        R_true[:, 0] *= -1
    for r in pred.residues:
        for a in list(r.atoms):
            r.atoms[a] = r.atoms[a] @ R_true + np.array([30.0, -12.0, 7.0]) \
                + rng.normal(0, 0.1, 3)

    result = state_recovery(
        pred,
        {"DFG-in": ref_in, "DFG-out": ref_out},
        region_seqids={"DFG-in": region, "DFG-out": region},
    )
    print("   ", result.summary().replace("\n", "\n    "))

    check("picks the correct state through a rigid transform",
          result.verdict == "DFG-out", f"got '{result.verdict}'")
    check("margin is decisive", result.decisive, f"margin={result.margin:.2f} A")

    # And the other way round, to be sure it isn't just always saying DFG-out.
    pred_in = build_kinase(d1=8.0, d2=16.0, seed=3)
    result_in = state_recovery(
        pred_in,
        {"DFG-in": ref_in, "DFG-out": ref_out},
        region_seqids={"DFG-in": region, "DFG-out": region},
    )
    check("picks DFG-in for a DFG-in prediction",
          result_in.verdict == "DFG-in", f"got '{result_in.verdict}'")


def test_ligand_metrics():
    print("\n--- ligand metrics ---")
    ref = build_kinase(d1=13.0, d2=10.0, seed=0, with_ligand=True,
                       ligand_at=[4.0, 2.0, 0.0])

    same = build_kinase(d1=13.0, d2=10.0, seed=0, with_ligand=True,
                        ligand_at=[4.0, 2.0, 0.0])
    d = ligand_centroid_distance(same, ref)
    check("ligand centroid distance ~0 for an identical pose",
          d is not None and d < 0.5, f"{d:.3f} A" if d is not None else "None")

    wrong = build_kinase(d1=13.0, d2=10.0, seed=0, with_ligand=True,
                         ligand_at=[4.0, 2.0, 25.0])
    d_wrong = ligand_centroid_distance(wrong, ref)
    check("ligand centroid distance is large for a wrong pose",
          d_wrong is not None and d_wrong > 20,
          f"{d_wrong:.3f} A" if d_wrong is not None else "None")

    contacts = pocket_contacts(ref)
    check("pocket contacts are found", len(contacts) > 0, f"{len(contacts)} residues")

    rec = contact_recovery(same, ref)
    check("contact recovery is 1.0 for an identical pose",
          rec["recall"] == 1.0, f"recall={rec['recall']}")

    rec_wrong = contact_recovery(wrong, ref)
    check("contact recovery drops for a wrong pose",
          rec_wrong["recall"] is None or rec_wrong["recall"] < 0.5,
          f"recall={rec_wrong['recall']}")


def test_file_roundtrip():
    print("\n--- file I/O roundtrip ---")
    st = build_kinase(d1=13.0, d2=10.0, with_ligand=True)
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "synthetic.pdb")
        write_pdb(st, p)
        back = load_structure(p)

        check("polymer residues survive the roundtrip",
              len(back.residues) == len(st.residues),
              f"{len(back.residues)} vs {len(st.residues)}")
        check("ligand survives the roundtrip",
              len(back.ligands()) == 1, f"{len(back.ligands())} ligands")

        res = classify_dfg(back)
        check("classification is identical after roundtrip",
              res.label == "DFG-out", f"got '{res.label}'")

        pairs = residue_correspondence(st, back)
        check("sequence correspondence covers the whole chain",
              len(pairs) == len(st.residues), f"{len(pairs)} pairs")


if __name__ == "__main__":
    print("=" * 66)
    print("Apheris multi-state kit -- self test")
    print("=" * 66)

    test_superposition()
    test_motif_detection()
    test_dfg_classification()
    test_motif_transfer_across_numbering()
    test_state_recovery()
    test_ligand_metrics()
    test_file_roundtrip()

    print("\n" + "=" * 66)
    if _failures:
        print(f"{len(_failures)} FAILED: " + ", ".join(_failures))
        sys.exit(1)
    print("all checks passed -- the kit is working")
    print("=" * 66)
