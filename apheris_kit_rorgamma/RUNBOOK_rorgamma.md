# RORgamma Runbook — Day 1 handoff (Lily -> Foundry ops, 17:15–17:30)

## Status at handoff (all gates closed)

| Gate | Deadline | Status |
|------|----------|--------|
| Scorer validation on crystals | 11:30 | PASSED 5/5 — 3KYT->H12-in, 4ZJW->H12-out, 5IXK->H12-out, 5NTK->H12-out, 4ZJR->H12-in (decoupling decoy correctly called, margin 2.33 A). `python test_scorer_gate.py` re-runs it in ~10 s. |
| Foundry schema | 12:15 | OPEN — probe on the laptop with Hub access (below). Emitted jobs use best-guess param names until then. |
| Reference pair | — | 3KYT (H12-in) vs 4ZJW (H12-out), chosen empirically; 4ZJR is the documented decoupling decoy. See `refs_rorgamma/ref_pair_analysis.json`, `displacement_profiles.png`. |

## Pre-registered success criterion (do not edit after results are known)

PRIMARY, inverse-agonist direction (`rorgamma_invago`): a steering arm (A1
pocket-conditioned, A3 target-template) recovers H12-out in >= 60% of seeds
with decisive margins (>= 1 A), AND A4 (wrong template) does not beat A0
(baseline) by more than +10 pp, AND A2 (shared pocket) sits within 10 pp of A0.
Steering means arms beat BOTH the baseline and the controls.

If A0 already recovers H12-out >= 60%: the honest headline is "the model
defaults to the therapeutic state" — still a result, not a steering claim.
If nothing recovers H12-out: negative result, report it as such.

SECONDARY: ligand centroid distance and pocket-contact recall are reported
NEXT TO state recovery, never merged (pose != state).

Threats to validity (named up front): contact-set leakage (pocket sets derive
from the same crystals as refs/templates — mitigated by A2/A4), template
memory (3KYT/4ZJW are pre-2025 — mitigated by A4/A5 and any post-cutoff
pharma-drop data), unvalidated H12 region (mitigated by the crystal gate +
sensitivity alternatives below).

## Core arms and seeds (agreed with ops review)

- Core arms BOTH directions: A0, A1, A2, A3, A4, A5, A8 (A2 = shared-pocket
  control — kept in the core: it is the only clean test that pocket
  conditioning per se does not flip the state).
- Seeds: 5 on every core arm in the `rorgamma_invago` direction (the recovery
  RATE is the headline; n=2 cannot support a rate), 2 on `rorgamma_25hc`.
- If the params schema supports samples-per-job, use it (one job, N samples);
  else submit one job per seed.
- Single model first: Boltz-2 arms only until the story is complete; B0–B3
  (OpenFold3 mirrors) only after.
- Compute: 7x5 + 7x2 = 49 predictions, ~1–3 min each on the Hub GPUs — well
  inside Build I–II. Quota-tight floor: 5 seeds only on A0/A1/A4 (invago) +
  2 seeds on A0/A1/A2 (25hc) = 27.

## Submission (Foundry ops)

    # 1. probe the live schema (laptop with apheris-foundry installed + authed)
    apheris-foundry workflows --help
    python foundry_emit.py probe            # re-emit if the schema differs
    python foundry_emit.py emit --runs runs/rorgamma_invago --out foundry_jobs_invago
    python foundry_emit.py emit --runs runs/rorgamma_25hc  --out foundry_jobs_25hc

    # 2. smoke arm FIRST (gate 2, by 12:15): one job end-to-end
    apheris-foundry workflows run --workflow predict \
      --input foundry_jobs_invago/A0_baseline.request.json \
      --model-params @foundry_jobs_invago/A0_baseline.params.json
    apheris-foundry jobs logs <job_id>
    apheris-foundry jobs download <job_id>

    # 3. then the whole grid as ONE batch
    bash foundry_jobs_invago/submit_all.sh
    bash foundry_jobs_25hc/submit_all.sh

If the schema rejects the pocket block: re-emit with the probed names; if the
Hub exposes NO pocket constraints, that arm list ("unsupported" in
foundry_manifest.json) is itself product feedback for Apheris — say it in the
pitch. If Foundry is entirely unavailable, fall back to the sandbox HPC
Boltz-2 and state so honestly.

## Collect + score (as results land)

    python foundry_emit.py collect --jobs foundry_jobs_invago/ --runs runs/rorgamma_invago
    python foundry_emit.py collect --jobs foundry_jobs_25hc/  --runs runs/rorgamma_25hc

    python score_run.py score --config targets_rorgamma.json \
      --target rorgamma_invago --runs runs/rorgamma_invago --out results_invago.csv
    python score_run.py score --config targets_rorgamma.json \
      --target rorgamma_25hc --runs runs/rorgamma_25hc --out results_25hc.csv

Predictions must land as `runs/<target>/<arm_id>/<model files>.cif`.

## Figure

    python figure.py results_invago.csv --out-prefix figure_day1_invago \
      --target-state H12-out --md summary_invago.md
    python figure.py results_25hc.csv --out-prefix figure_day1_25hc \
      --target-state H12-in --md summary_25hc.md

## Sensitivity checks (Build II, after the main table stands)

1. Region alternatives (edit `state_region_seqids` in targets_rorgamma.json,
   re-score only): auto region `[286,287,288,290,479..486]`; full H12 window
   `range(470, 508)`.
2. Alternative inactive reference: swap `ref_other_state`/`ref_target_state`
   to `refs_rorgamma/5IXK.cif` (chain B) or `5NTK.cif` (chain A), re-score.
3. A9_peptide_cofold (25hc grid): if the SRC2-2 peptide lands in the groove
   AND H12 packs, that is orthogonal evidence the agonist state assembled —
   the slide's coactivator, made computational.

## Raw prediction outputs are not tracked

`runs/**/*_model.cif` (the raw OpenFold3 structures themselves) are
git-ignored — they regenerate from this runbook's own commands, run
against the pinned image digests in `collect_compute_stats.sh`/
`h2_compute_usage.md`, so tracking them a second time is pure
duplication. `RAW_PREDICTIONS_MANIFEST.csv` lists every one that existed
at commit time (relative path, size, sha256), so a missing or
corrupted re-run is detectable without having to carry the bytes.

## Pitch line

"State-specific pocket conditioning steers co-folding models between RORgamma
agonist and inverse-agonist H12 states — and the controls show naive
conditioning can't." Supporting facts: the agonist's state-specific contacts
include the agonist lock itself (H479, Y502); the inverse agonist's reach
toward H12 (480) and the back pocket; the scorer calls the 4ZJR decoupling
decoy honestly; 4ZJR/9VZQ show ligand labels lie about state — we derived the
pair from displacement profiles instead.
