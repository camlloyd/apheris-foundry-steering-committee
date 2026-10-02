# Task 1 — decision-boundary recalibration: results

**Pre-registered:** `../preregistrations/H5_TASK1_BOUNDARY_PREREGISTRATION.md`,
rule (a) committed before this ran.
**Scope:** rescore of the existing 258-compound pool (`~/h3/official_state_calls_all.csv`),
no new predictions.

## Result

`D_ref` (mutual H12-crop backbone RMSD between the two official references,
5VB7 and 6T4I) = **0.353 Å**, symmetric both directions (0.35302 vs 0.35302,
agreeing to 5 decimal places — floating-point noise only).

As predicted analytically in the pre-registration (a symmetric distance
metric scored against exactly two reference points places them
symmetrically around `delta = 0` by construction, for any `D_ref`), rule
(a) is **bit-for-bit identical to the current rule**:

| rule | balanced accuracy | pooled | agonist | antagonist |
|---|---|---|---|---|
| current (delta > 0 → agonist) | 69.1% | 69.0% | 77.3% | 60.8% |
| (a) recentered on D_ref midpoint | 69.1% | 69.0% | 77.3% | 60.8% |

**Per the pre-registered decision rule: does not improve balanced accuracy
(it cannot — it is the same rule), so the original numbers stand
unchanged.** The 77.3%/60.8% split is not an artifact of the boundary
being mis-centered relative to the two reference crystals.

## Post-hoc boundary sweep (sensitivity analysis, not a result)

Swept the raw decision variable `delta = rmsd_to_antagonist - rmsd_to_agonist`
over [-6, +6] Å in steps of 0.25, full CSV in `task1_boundary_sweep.csv`.

- Balanced accuracy peaks at **69.06%, at `delta = 0.0`** — exactly where
  the current rule already sits. A secondary near-peak exists at
  `delta ≈ -0.75` (68.8%, agonist 96.9%/antagonist 40.8%), but it is lower,
  not a different regime worth adopting.
- The curve is not a broad plateau: balanced accuracy falls from 69% to
  ~52% within 0.25 Å of either side of the optimum, and degrades further,
  monotonically, out to the extremes (all-agonist floor ~50%/all-antagonist
  floor ~50%, as expected from class balance).
- **Conclusion:** the existing boundary already sits at (or within noise
  of) the empirical optimum for this pool. There is no "free" accuracy
  being left on the table by a mis-placed threshold — moving it in either
  direction trades agonist accuracy for antagonist accuracy roughly
  one-for-one, it does not improve both. This is consistent with (2) from
  the pre-registration's framing being false and (1) being the right
  read: **antagonist really is the harder class for this model on real
  ligands, not a scoring-threshold artifact.**

## Provenance

`task1_boundary_recalibration.py`, run 2026-10-02, reading only
`~/h3/official_state_calls_all.csv` (already on disk) plus one cheap
non-GPU `apheris-data` call to score the two 24-residue reference crops
against each other for `D_ref`.
