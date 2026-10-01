"""state_recovery2 - CA displacement profiling, state-region derivation and the
generic two-state recovery scorer.

The state region is *derived from the data*, not asserted.  Given an active and
an inactive reference, positions whose CA moves more than
``max(abs_floor, median_displacement)`` are collected, merged into contiguous
segments, and short segments are dropped.  For RORgamma this is expected to
recover helix 12 / AF-2; if it instead recovers a scattered set of loops, that
is a signal the chosen reference pair does not encode a clean two-state switch
and the fallback in the plan's section 7 applies.

The region is always expressed in the numbering of the ACTIVE reference, which
acts as the canonical frame.  Every other structure reaches it through a
sequence alignment, so no numbering assumption leaks between entries.
"""

from __future__ import annotations

import dataclasses
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .structio import ChainModel, align_pair, rmsd, superpose_on

ABS_FLOOR_A = 2.0       # the plan's max(2 A, median) threshold
MERGE_GAP = 2           # bridge gaps of <= this many positions when merging
MIN_SEGMENT = 4         # drop segments shorter than this (noise, not a helix)
MAX_ITERS = 6


@dataclasses.dataclass
class Profile:
    """Per-position CA displacement between two aligned chains."""

    numbers_a: List[int]
    numbers_b: List[int]
    disp: np.ndarray            # (N,) Angstrom, parallel to numbers_a
    core_mask: np.ndarray       # (N,) bool - positions used for the final fit
    core_rmsd: float
    overall_rmsd: float
    n_aligned: int
    iterations: int


@dataclasses.dataclass
class Region:
    """A derived state region, in ACTIVE-reference numbering."""

    numbers: List[int]
    segments: List[Tuple[int, int]]
    threshold: float
    median_disp: float
    max_disp: float
    mean_disp_in_region: float

    def as_dict(self) -> dict:
        return {
            "numbers": self.numbers,
            "segments": [list(s) for s in self.segments],
            "threshold_A": round(self.threshold, 3),
            "median_disp_A": round(self.median_disp, 3),
            "max_disp_A": round(self.max_disp, 3),
            "mean_disp_in_region_A": round(self.mean_disp_in_region, 3),
            "n_positions": len(self.numbers),
        }


def displacement_profile(
    a: ChainModel,
    b: ChainModel,
    abs_floor: float = ABS_FLOOR_A,
    max_iters: int = MAX_ITERS,
) -> Profile:
    """CA displacement of `b` relative to `a` after core-only superposition.

    The fit set is refined iteratively: an all-position fit gives a first
    displacement profile, the high-displacement positions are excluded, and the
    fit is repeated.  This stops a large rigid-body swing of one helix from
    rotating the whole domain and hiding itself.
    """
    pairs = align_pair(a, b)
    if len(pairs) < 20:
        raise ValueError(
            f"only {len(pairs)} aligned positions between {a.source} and "
            f"{b.source}; refusing to profile"
        )
    ia = np.array([p[0] for p in pairs])
    ib = np.array([p[1] for p in pairs])
    ref = a.ca[ia]
    mob = b.ca[ib]

    core = np.ones(len(pairs), dtype=bool)
    disp = None
    iters = 0
    for iters in range(1, max_iters + 1):
        fitted = superpose_on(mob, ref, core)
        disp = np.linalg.norm(fitted - ref, axis=1)
        thr = max(abs_floor, float(np.median(disp)))
        new_core = disp <= thr
        if new_core.sum() < max(10, int(0.3 * len(pairs))):
            # Refinement is eating the structure - keep the previous fit set.
            break
        if np.array_equal(new_core, core):
            core = new_core
            break
        core = new_core

    fitted = superpose_on(mob, ref, core)
    disp = np.linalg.norm(fitted - ref, axis=1)
    return Profile(
        numbers_a=[a.numbers[i] for i in ia],
        numbers_b=[b.numbers[i] for i in ib],
        disp=disp,
        core_mask=core,
        core_rmsd=rmsd(fitted[core], ref[core]),
        overall_rmsd=rmsd(fitted, ref),
        n_aligned=len(pairs),
        iterations=iters,
    )


def _segments(numbers: Sequence[int], merge_gap: int = MERGE_GAP) -> List[Tuple[int, int]]:
    """Merge a sorted list of residue numbers into (start, end) segments."""
    if not numbers:
        return []
    nums = sorted(numbers)
    segs = [[nums[0], nums[0]]]
    for n in nums[1:]:
        if n - segs[-1][1] <= merge_gap + 1:
            segs[-1][1] = n
        else:
            segs.append([n, n])
    return [(s[0], s[1]) for s in segs]


def derive_region(
    profile: Profile,
    abs_floor: float = ABS_FLOOR_A,
    merge_gap: int = MERGE_GAP,
    min_segment: int = MIN_SEGMENT,
) -> Region:
    """Auto-derive the state region from a displacement profile."""
    disp = profile.disp
    median = float(np.median(disp))
    thr = max(abs_floor, median)
    hot = [n for n, d in zip(profile.numbers_a, disp) if d > thr]
    segs = _segments(hot, merge_gap=merge_gap)

    kept_segs, kept_nums = [], []
    present = set(profile.numbers_a)
    for lo, hi in segs:
        span = [n for n in range(lo, hi + 1) if n in present]
        if len(span) >= min_segment:
            kept_segs.append((lo, hi))
            kept_nums.extend(span)

    if kept_nums:
        idx = {n: i for i, n in enumerate(profile.numbers_a)}
        in_region = np.array([disp[idx[n]] for n in kept_nums])
        mean_in = float(in_region.mean())
    else:
        mean_in = 0.0

    return Region(
        numbers=sorted(kept_nums),
        segments=kept_segs,
        threshold=thr,
        median_disp=median,
        max_disp=float(disp.max()),
        mean_disp_in_region=mean_in,
    )


def _mapped_indices(
    canonical: ChainModel,
    other: ChainModel,
    region_numbers: Sequence[int],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Indices into (canonical, other) for region positions and core positions."""
    pairs = align_pair(canonical, other)
    want = set(region_numbers)
    reg_c, reg_o, core_c, core_o = [], [], [], []
    for i, j in pairs:
        if canonical.numbers[i] in want:
            reg_c.append(i)
            reg_o.append(j)
        else:
            core_c.append(i)
            core_o.append(j)
    return (
        np.array(reg_c, dtype=int),
        np.array(reg_o, dtype=int),
        np.array(core_c, dtype=int),
        np.array(core_o, dtype=int),
    )


def region_rmsd(
    pred: ChainModel,
    ref: ChainModel,
    region_numbers_in_canonical: Sequence[int],
    canonical: ChainModel,
) -> Dict[str, float]:
    """RMSD of the state region between `pred` and `ref`, after core-only fit.

    All three structures are tied together through `canonical` (the active
    reference), so the region means the same set of residues everywhere even
    when author numbering differs.
    """
    # Region numbers expressed in the reference's own numbering.
    _, reg_ref_idx, _, _ = _mapped_indices(canonical, ref, region_numbers_in_canonical)
    ref_region_numbers = {ref.numbers[j] for j in reg_ref_idx}

    pairs = align_pair(ref, pred)
    reg_r, reg_p, core_r, core_p = [], [], [], []
    for i, j in pairs:
        if ref.numbers[i] in ref_region_numbers:
            reg_r.append(i)
            reg_p.append(j)
        else:
            core_r.append(i)
            core_p.append(j)

    if len(core_r) < 10 or len(reg_r) < 3:
        return {
            "region_rmsd": float("nan"),
            "core_rmsd": float("nan"),
            "n_region": len(reg_r),
            "n_core": len(core_r),
        }

    ref_all = ref.ca[np.array(core_r + reg_r)]
    pred_all = pred.ca[np.array(core_p + reg_p)]
    fit_mask = np.zeros(len(ref_all), dtype=bool)
    fit_mask[: len(core_r)] = True
    fitted = superpose_on(pred_all, ref_all, fit_mask)

    return {
        "region_rmsd": rmsd(fitted[len(core_r):], ref_all[len(core_r):]),
        "core_rmsd": rmsd(fitted[: len(core_r)], ref_all[: len(core_r)]),
        "n_region": len(reg_r),
        "n_core": len(core_r),
    }


def state_recovery_generic(
    pred: ChainModel,
    refs: Dict[str, ChainModel],
    region_numbers: Sequence[int],
    canonical_key: str = "active",
    decisive_margin: float = 1.0,
) -> Dict[str, object]:
    """Score one prediction against an active/inactive reference pair.

    Returns the region RMSD to each reference, the signed margin
    (d_active - d_inactive; negative means the prediction looks active) and a
    categorical call that is withheld when the margin is below
    `decisive_margin`.
    """
    canonical = refs[canonical_key]
    out: Dict[str, object] = {}
    dists: Dict[str, float] = {}
    for name, ref in refs.items():
        res = region_rmsd(pred, ref, region_numbers, canonical)
        dists[name] = res["region_rmsd"]
        out[f"region_rmsd_{name}"] = res["region_rmsd"]
        out[f"core_rmsd_{name}"] = res["core_rmsd"]
        out[f"n_region_{name}"] = res["n_region"]

    d_a = dists.get("active", float("nan"))
    d_i = dists.get("inactive", float("nan"))
    margin = d_a - d_i
    out["margin"] = margin
    if np.isnan(margin):
        out["call"] = "unscored"
    elif abs(margin) < decisive_margin:
        out["call"] = "non-decisive"
    else:
        out["call"] = "active" if margin < 0 else "inactive"
    out["decisive_margin"] = decisive_margin
    return out
