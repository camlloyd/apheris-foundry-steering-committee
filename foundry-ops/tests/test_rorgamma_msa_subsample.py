import pytest

import rorgamma_msa_subsample as ms

QUERY_ENTRY = ">101\nASLTEIEHLVQ\n"
HIT_1 = ">hit1\tscore1\nASLTEIEHLVQ\n"
HIT_2 = ">hit2\tscore2\nAS-TEIEHLVQ\n"
HIT_3 = ">hit3\tscore3\nASLTEI-HLVQ\n"
A3M = QUERY_ENTRY + HIT_1 + HIT_2 + HIT_3


def test_subsample_keeps_query_entry_first():
    out = ms.subsample_a3m(A3M, depth=1, seed=0)
    assert out.startswith(QUERY_ENTRY)


def test_subsample_to_depth_1_keeps_only_the_query():
    out = ms.subsample_a3m(A3M, depth=1, seed=0)
    entries = out.strip().split("\n>")
    assert len(entries) == 1


def test_subsample_to_depth_n_keeps_query_plus_n_minus_1_hits():
    out = ms.subsample_a3m(A3M, depth=2, seed=0)
    entries = [e for e in out.split(">") if e]
    assert len(entries) == 2
    assert out.startswith(QUERY_ENTRY)


def test_subsample_depth_greater_than_available_keeps_everything():
    out = ms.subsample_a3m(A3M, depth=100, seed=0)
    assert out == A3M


def test_subsample_is_deterministic_for_a_given_seed():
    out1 = ms.subsample_a3m(A3M, depth=2, seed=7)
    out2 = ms.subsample_a3m(A3M, depth=2, seed=7)
    assert out1 == out2


def test_subsample_raises_on_depth_less_than_one():
    with pytest.raises(ValueError):
        ms.subsample_a3m(A3M, depth=0, seed=0)
