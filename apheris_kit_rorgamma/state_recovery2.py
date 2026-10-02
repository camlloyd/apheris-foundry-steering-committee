"""
Target-agnostic two-state scoring: state recovery when your target is NOT
necessarily a kinase.

The core kit (state_recovery.py) needs you to say which residues define the
state -- fine for kinases (DFG window), useless if the pharma data drop turns
out to be GPCRs, nuclear receptors, transporters, or anything else. This
module derives the state-defining region directly from the two reference
structures:

    superpose the two references on ALL paired CA atoms
    -> per-residue displacement
    -> the state-defining region = the residues that actually moved

That is the only honest way to do it: the region comes from the data, not
from a motif annotation that may not exist for your target class.

Everything else (dual-reference verdict, margin, decisiveness) is inherited
from state_recovery.py -- this module just fills in the region automatically.

Usage:
    from state_recovery2 import derive_state_region, state_recovery_generic
    reg = derive_state_region(ref_a, ref_b)
    result = state_recovery_generic(pred, {"stateA": ref_a, "stateB": ref_b})
"""
from __future__ import annotations

import numpy as np

from structure_io import Structure
from state_recovery import (
    StateRecoveryResult, apply_transform, kabsch, residue_correspondence,
    state_recovery,
)


def displacement_profile(
    ref_a: Structure,
    ref_b: Structure,
    a_chain: str | None = None,
    b_chain: str | None = None,
) -> list:
    """
    Per-residue CA displacement between two references after whole-structure
    superposition. Returns [(seqid_in_a, displacement_A), ...], sorted by
    seqid. This is the raw signal everything below is built on, and it is
    worth eyeballing (or plotting) before trusting the derived region.
    """
    pairs = residue_correspondence(ref_a, ref_b, a_chain, b_chain)
    if len(pairs) < 10:
        return []

    a_xyz = np.array([ra.atom("CA") for ra, _ in pairs])
    b_xyz = np.array([rb.atom("CA") for _, rb in pairs])
    R, t = kabsch(b_xyz, a_xyz)
    b_fit = apply_transform(b_xyz, R, t)
    disp = np.linalg.norm(b_fit - a_xyz, axis=1)

    return sorted(
        ((ra.seqid, float(d)) for (ra, _), d in zip(pairs, disp)),
        key=lambda x: x[0],
    )


def _segments(seqids: list, max_gap: int = 3) -> list:
    """Group sorted seqids into contiguous runs, bridging small gaps."""
    if not seqids:
        return []
    segs = [[seqids[0]]]
    for s in seqids[1:]:
        if s - segs[-1][-1] <= max_gap + 1:
            segs[-1].append(s)
        else:
            segs.append([s])
    return segs


def derive_state_region(
    ref_a: Structure,
    ref_b: Structure,
    min_disp: float = 2.0,
    min_region: int = 6,
    min_segment: int = 3,
    a_chain: str | None = None,
    b_chain: str | None = None,
) -> dict:
    """
    The state-defining region between two reference structures, in ref_a
    numbering.

    A residue counts as "moved" if its CA displacement after global
    superposition exceeds max(min_disp, median displacement). Moved residues
    are merged into segments (small gaps bridged, tiny segments dropped as
    noise). If the segments are too small overall, fall back to the top
    min_region residues by displacement -- a small region is less sensitive
    but still better than a global RMSD call.

    Returns a dict with the region set, the threshold used, the full
    displacement profile, and warnings. READ the warnings: if the two
    references barely differ, there is no state to recover and every
    downstream number is noise.
    """
    profile = displacement_profile(ref_a, ref_b, a_chain, b_chain)
    out = {"region": set(), "threshold": None, "profile": profile,
           "segments": [], "warnings": []}
    if not profile:
        out["warnings"].append(
            "fewer than 10 aligned CA pairs between the two references -- "
            "are these the same protein and construct?"
        )
        return out

    disps = np.array([d for _, d in profile])
    thr = max(min_disp, float(np.median(disps)))
    out["threshold"] = thr

    moved = sorted(s for s, d in profile if d > thr)
    segs = [sg for sg in _segments(moved) if len(sg) >= min_segment]
    region = {s for sg in segs for s in sg}

    if len(region) < min_region:
        top = sorted(profile, key=lambda x: -x[1])[:min_region]
        region = {s for s, _ in top}
        segs = _segments(sorted(region))
        out["warnings"].append(
            f"only {len(moved)} residues moved above {thr:.1f} A; fell back to "
            f"the top {min_region} most-displaced residues. The two references "
            f"may be the same state -- check the profile before scoring."
        )

    if float(np.max(disps)) < min_disp:
        out["warnings"].append(
            f"maximum CA displacement between the references is "
            f"{float(np.max(disps)):.2f} A -- these look like the SAME state. "
            f"State recovery on this pair is meaningless."
        )

    out["region"] = region
    out["segments"] = segs
    return out


def map_region(region: set, ref_a: Structure, ref_b: Structure,
               a_chain=None, b_chain=None) -> set:
    """Translate a region in ref_a numbering into ref_b numbering."""
    pairs = residue_correspondence(ref_a, ref_b, a_chain, b_chain)
    a_to_b = {ra.seqid: rb.seqid for ra, rb in pairs}
    return {a_to_b[s] for s in region if s in a_to_b}


def state_recovery_generic(
    prediction: Structure,
    references: dict,
    region: set | None = None,
    pred_chain: str | None = None,
    ref_chains: dict | None = None,
) -> StateRecoveryResult:
    """
    Drop-in replacement for state_recovery() that derives the state-defining
    region itself when you don't have one.

    references maps {state_name: Structure}. With two references the region
    comes from their mutual displacement; with more, it is the union over all
    pairs. The region is mapped into each reference's own numbering before
    scoring, exactly as state_recovery() expects.

    Pass an explicit region (in the FIRST reference's numbering) if you have
    prior knowledge -- e.g. you know helix 12 is the whole story for your
    nuclear receptor and don't want the activation-loop wobble diluting it.
    """
    names = list(references)
    if len(names) < 2:
        raise ValueError("need at least two reference states")

    notes = []
    if region is None:
        union = None
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                der = derive_state_region(references[names[i]], references[names[j]])
                notes.extend(f"[{names[i]} vs {names[j]}] {w}"
                             for w in der["warnings"])
                union = der["region"] if union is None else union | der["region"]
        region = union
        notes.append(
            f"auto-derived state region: {len(region)} residues in "
            f"{names[0]} numbering"
        )

    # map the region into every reference's own numbering
    first = references[names[0]]
    region_seqids = {names[0]: set(region)}
    for name in names[1:]:
        region_seqids[name] = map_region(region, first, references[name])

    result = state_recovery(
        prediction, references,
        region_seqids=region_seqids,
        pred_chain=pred_chain, ref_chains=ref_chains,
    )
    result.notes = notes + result.notes
    return result
