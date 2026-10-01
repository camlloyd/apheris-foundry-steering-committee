"""Subsample an .a3m alignment to a shallower depth.

Tests whether partial MSA perturbation (not full removal, which just
destroys fold quality -- see the "Very Low" pLDDT MSA-free arms) unlocks
ligand-induced state steering the way it's been shown to unlock alternative
apo-state sampling (AF-Cluster and friends) -- but never tested with a
ligand in the picture, per apheris_kit_rorgamma/README.md's own framing of
the literature gap.
"""
from __future__ import annotations

import random


def subsample_a3m(text: str, depth: int, seed: int) -> str:
    """Keep the query entry (first record) plus a random sample of `depth -
    1` hit entries. Deterministic for a given seed. If fewer hits exist than
    requested, keeps all of them."""
    if depth < 1:
        raise ValueError(f"depth must be >= 1, got {depth}")

    # Each record starts with ">" and runs to the next ">" (or end of string).
    records = [">" + chunk for chunk in text.split(">")[1:]]

    query, hits = records[0], records[1:]
    n_hits_wanted = depth - 1
    if n_hits_wanted >= len(hits):
        return text

    rng = random.Random(seed)
    chosen = rng.sample(hits, n_hits_wanted)
    return query + "".join(chosen)
