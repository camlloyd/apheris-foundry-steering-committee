# Task 1 pre-registration — decision-boundary recalibration

**Timestamped:** 2026-10-02, before any rescore under the new rule ran.
**Scope:** rescore of the existing 258-compound ChEMBL pool (120 curated +
138 scored-of-150 random) already on disk. No new GPU jobs, no new
predictions — this only changes how an existing `(rmsd_to_agonist,
rmsd_to_antagonist)` pair per structure is turned into a call.

## Question

The current agonist/antagonist split (77.3% / 60.8%) could mean two very
different things: (1) antagonist structures are a genuinely harder class for
this model to get right, or (2) the decision boundary itself is
mis-centered, so some genuine antagonist predictions are being scored as
agonist by an arbitrary cutoff, not because the model failed to represent
the state. Only (2) is fixable without touching the model. This task tests
(2) by recalibrating the boundary and nothing else.

## Candidates considered, one chosen in advance

- **(a) midpoint of the two reference crystals' mutual H12 displacement.**
  Score the agonist reference (5VB7) against the antagonist reference
  (6T4I) through the same `h12_state.call_states` path used for every
  other structure, giving one number, `D_ref` = their mutual H12-crop
  backbone RMSD. Recenter the decision variable
  `delta = rmsd_to_antagonist - rmsd_to_agonist` on whatever `D_ref`
  implies the references' own midpoint to be, rather than assuming it
  sits at `delta = 0`.
- **(b) per-distance z-scoring against each reference's own within-state
  spread across a reference set**, then compare the two z-scores.

**Chosen: (a).** It uses only the two reference structures the existing
pipeline already depends on — no second reference set has to be assembled
or its membership justified, and no spread/variance estimate has to be
fit. (b) requires deciding what counts as "the reference set" for each
state (how many structures, which ones), which is itself a free parameter
this project's pre-registration discipline exists to avoid introducing
on judgment night. (a) has none.

## Mechanical prediction, stated before running it

Because `delta` is computed from a backbone RMSD that is symmetric
(`RMSD(X,Y) == RMSD(Y,X)` under this metric), scoring the two references
against each other gives `rmsd_to_agonist(6T4I) == rmsd_to_antagonist(5VB7)
== D_ref` and `rmsd_to_agonist(5VB7) == rmsd_to_antagonist(6T4I) == 0`. That
places the two references at `delta = +D_ref` (agonist) and
`delta = -D_ref` (antagonist) — already symmetric around `delta = 0` by
construction, for any two-reference pair, regardless of what `D_ref` turns
out to be. **Rule (a)'s boundary is therefore algebraically identical to
the current rule** (`delta > 0 -> agonist`, i.e. whichever reference is
closer). This is predicted analytically, before measuring `D_ref` at all;
measuring it is still worth doing as a mechanical check that the pipeline
behaves as the algebra says, and because the sweep in the second half of
this task needs `D_ref` as a reference point on the x-axis regardless.

## Decision rule, decided in advance

If the empirical rescore matches the analytic prediction (no change),
report that explicitly as the finding — "the 77.3/60.8 split is not a
boundary-placement artifact of this type" — and keep the original numbers
on the slide. Do not try candidate (b) afterward and present whichever
wins; (b) was ruled out above for reasons independent of the outcome here.
