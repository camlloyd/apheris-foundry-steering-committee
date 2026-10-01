"""Metrics the judging rubric asks for that the apheris-data compute-metrics
Foundry task doesn't expose under these exact names: lDDT-PLI (protein-ligand
interface lDDT) and Ligand BiSyRMSD (symmetry-corrected ligand pose RMSD,
see rorgamma_ligand_rmsd.py for the latter -- this module is lDDT-PLI only).
"""
from __future__ import annotations

import numpy as np

DEFAULT_THRESHOLDS = (0.5, 1.0, 2.0, 4.0)


def lddt_pli(
    ref_protein: np.ndarray,
    ref_ligand: np.ndarray,
    pred_protein: np.ndarray,
    pred_ligand: np.ndarray,
    cutoff: float = 6.0,
    thresholds: tuple[float, ...] = DEFAULT_THRESHOLDS,
) -> float | None:
    """Protein-ligand interface lDDT.

    Considers every (ligand atom, protein atom) pair whose distance in the
    *reference* is within `cutoff` -- not just index-matched pairs, since a
    ligand atom can be close to many protein atoms and vice versa. For each
    such pair, compares its reference distance to its distance in the
    prediction and scores the fraction preserved within each threshold in
    `thresholds`, then averages across thresholds and pairs.

    Protein/ligand atom arrays must be in corresponding order between the
    reference and prediction (same molecule, same atom ordering, different
    coordinates) -- this does not do atom matching itself.

    Returns None if no interface pairs exist within `cutoff` (nothing to
    score).
    """
    ref_protein = np.asarray(ref_protein, dtype=float)
    ref_ligand = np.asarray(ref_ligand, dtype=float)
    pred_protein = np.asarray(pred_protein, dtype=float)
    pred_ligand = np.asarray(pred_ligand, dtype=float)

    ref_dist = np.linalg.norm(ref_ligand[:, None, :] - ref_protein[None, :, :], axis=-1)
    pred_dist = np.linalg.norm(pred_ligand[:, None, :] - pred_protein[None, :, :], axis=-1)

    mask = ref_dist <= cutoff
    if not mask.any():
        return None

    errors = np.abs(pred_dist[mask] - ref_dist[mask])
    preserved_fractions = [float(np.mean(errors < t)) for t in thresholds]
    return float(np.mean(preserved_fractions))
