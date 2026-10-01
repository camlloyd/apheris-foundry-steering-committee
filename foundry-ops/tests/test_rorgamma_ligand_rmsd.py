import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

import rorgamma_ligand_rmsd as rlr


def _embedded_mol(smiles, seed=42):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    return mol


def test_bisy_rmsd_is_zero_for_an_identical_conformer():
    mol = _embedded_mol("CCO")  # ethanol, no heavy-atom symmetry needed here
    score = rlr.ligand_bisy_rmsd(mol, Chem.Mol(mol))
    assert score == pytest.approx(0.0, abs=1e-3)


def test_bisy_rmsd_is_zero_after_a_symmetric_atom_swap():
    # CF4: all four F atoms are chemically equivalent (tetrahedral symmetry).
    # Swapping two F atoms' 3D positions produces the identical point cloud,
    # just under a different atom-index assignment -- BiSyRMSD must recognize
    # this as the same pose (score ~0), unlike a naive index-order RMSD.
    mol = _embedded_mol("FC(F)(F)F")
    pred = Chem.Mol(mol)
    conf = pred.GetConformer()
    f_indices = [a.GetIdx() for a in pred.GetAtoms() if a.GetSymbol() == "F"]
    i, j = f_indices[0], f_indices[1]
    pos_i, pos_j = conf.GetAtomPosition(i), conf.GetAtomPosition(j)
    conf.SetAtomPosition(i, pos_j)
    conf.SetAtomPosition(j, pos_i)

    score = rlr.ligand_bisy_rmsd(mol, pred)
    assert score == pytest.approx(0.0, abs=1e-3)


def test_naive_index_rmsd_would_have_been_nonzero_for_the_same_swap():
    # Sanity check that the swap in the test above is a real test of symmetry
    # handling, not a no-op -- a plain index-matched RMSD on the same swapped
    # coordinates must be clearly nonzero. GetConformerRMS compares two
    # conformers of the SAME mol, so add the swapped positions as a second
    # conformer rather than a separate Mol.
    mol = _embedded_mol("FC(F)(F)F")
    swapped_conf = Chem.Conformer(mol.GetConformer())
    f_indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "F"]
    i, j = f_indices[0], f_indices[1]
    pos_i, pos_j = swapped_conf.GetAtomPosition(i), swapped_conf.GetAtomPosition(j)
    swapped_conf.SetAtomPosition(i, pos_j)
    swapped_conf.SetAtomPosition(j, pos_i)
    swapped_conf_id = mol.AddConformer(swapped_conf, assignId=True)

    naive = AllChem.GetConformerRMS(mol, 0, swapped_conf_id, prealigned=True)
    assert naive > 0.5


def test_bisy_rmsd_is_nonzero_for_a_genuinely_different_conformation():
    mol1 = _embedded_mol("CCCCCC", seed=1)
    mol2 = _embedded_mol("CCCCCC", seed=999)
    score = rlr.ligand_bisy_rmsd(mol1, mol2)
    assert score > 0.1
