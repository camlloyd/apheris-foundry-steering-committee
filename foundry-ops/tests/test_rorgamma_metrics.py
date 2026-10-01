import numpy as np
import pytest

import rorgamma_metrics as rm


def test_lddt_pli_perfect_prediction_scores_one():
    ref_protein = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    ref_ligand = np.array([[1.0, 0.0, 0.0]])
    # prediction is an exact rigid-body copy (translated) -> all interface
    # distances preserved exactly
    shift = np.array([5.0, 5.0, 5.0])
    pred_protein = ref_protein + shift
    pred_ligand = ref_ligand + shift

    score = rm.lddt_pli(ref_protein, ref_ligand, pred_protein, pred_ligand)
    assert score == pytest.approx(1.0)


def test_lddt_pli_only_considers_pairs_within_cutoff():
    # protein atom 0 is within cutoff of the ligand (dist=1), atom 1 is not
    # (dist=10, cutoff=6) -> only atom 0 contributes, and distorting atom 1's
    # position in the prediction must not affect the score.
    ref_protein = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    ref_ligand = np.array([[1.0, 0.0, 0.0]])
    pred_protein = np.array([[0.0, 0.0, 0.0], [999.0, 0.0, 0.0]])  # atom 1 moved far
    pred_ligand = np.array([[1.0, 0.0, 0.0]])  # unchanged, still perfect for atom 0

    score = rm.lddt_pli(ref_protein, ref_ligand, pred_protein, pred_ligand, cutoff=6.0)
    assert score == pytest.approx(1.0)


def test_lddt_pli_distorted_prediction_scores_less_than_one():
    ref_protein = np.array([[0.0, 0.0, 0.0]])
    ref_ligand = np.array([[1.0, 0.0, 0.0]])  # ref distance = 1.0
    pred_protein = np.array([[0.0, 0.0, 0.0]])
    pred_ligand = np.array([[10.0, 0.0, 0.0]])  # pred distance = 10.0, off by 9

    score = rm.lddt_pli(ref_protein, ref_ligand, pred_protein, pred_ligand)
    # error of 9 exceeds every threshold (0.5,1,2,4) -> 0/4 preserved
    assert score == pytest.approx(0.0)


def test_lddt_pli_partial_preservation_across_thresholds():
    # single interface pair, error of exactly 1.5 -> preserved at thresholds
    # 2 and 4, not at 0.5 and 1 -> score = 2/4 = 0.5
    ref_protein = np.array([[0.0, 0.0, 0.0]])
    ref_ligand = np.array([[1.0, 0.0, 0.0]])  # ref distance = 1.0
    pred_protein = np.array([[0.0, 0.0, 0.0]])
    pred_ligand = np.array([[2.5, 0.0, 0.0]])  # pred distance = 2.5, error = 1.5

    score = rm.lddt_pli(ref_protein, ref_ligand, pred_protein, pred_ligand)
    assert score == pytest.approx(0.5)


def test_lddt_pli_no_interface_pairs_returns_none():
    ref_protein = np.array([[0.0, 0.0, 0.0]])
    ref_ligand = np.array([[100.0, 0.0, 0.0]])  # far outside any cutoff
    pred_protein = np.array([[0.0, 0.0, 0.0]])
    pred_ligand = np.array([[100.0, 0.0, 0.0]])

    assert rm.lddt_pli(ref_protein, ref_ligand, pred_protein, pred_ligand, cutoff=6.0) is None


def test_lddt_pli_averages_over_multiple_pairs():
    # two independent interface pairs: one perfect, one fully wrong (error 9)
    # -> average of 1.0 and 0.0 = 0.5
    ref_protein = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    ref_ligand = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    pred_protein = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    pred_ligand = np.array([[1.0, 0.0, 0.0], [10.0, 0.0, 0.0]])

    # NB: this treats the two ligand "atoms" independently paired with their
    # matching protein atom by index -- see lddt_pli's docstring for the
    # pairing convention (all ligand atoms x all protein atoms, not just
    # index-matched), so compute the real expected value from that same
    # all-pairs convention instead of assuming index-pairing here.
    score = rm.lddt_pli(ref_protein, ref_ligand, pred_protein, pred_ligand)
    # all-pairs: 2 protein atoms x 2 ligand atoms = 4 pairs, all within cutoff
    # (every ref distance is 1.0). Pred distances: ligand[0] vs both protein
    # atoms = 1.0 (perfect, preserved at all 4 thresholds); ligand[1] vs both
    # protein atoms = 10.0 (error 9, preserved at 0 thresholds).
    # -> (4 + 4 + 0 + 0) / (4 pairs * 4 thresholds) = 8/16 = 0.5
    assert score == pytest.approx(0.5)
