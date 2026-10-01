"""Ligand BiSyRMSD -- symmetry-corrected ligand pose RMSD, the
AlphaFold3/PoseBusters convention. See rorgamma_metrics.py for the other
missing judging-rubric metric (lDDT-PLI).

Uses RDKit's `rdMolAlign.GetBestRMS`, which enumerates the ligand's graph
automorphisms (e.g. a CF4's four fluorines are interchangeable) and returns
the best-matching RMSD after alignment, rather than a naive index-order RMSD
that would be fooled by chemically-equivalent atoms sitting at swapped
positions.
"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import rdMolAlign


def ligand_bisy_rmsd(ref_mol: Chem.Mol, pred_mol: Chem.Mol) -> float:
    """Symmetry-corrected RMSD between a reference and predicted ligand
    conformer. Both must be RDKit Mols with the same connectivity (same
    molecule) and an embedded 3D conformer."""
    return rdMolAlign.GetBestRMS(pred_mol, ref_mol)
