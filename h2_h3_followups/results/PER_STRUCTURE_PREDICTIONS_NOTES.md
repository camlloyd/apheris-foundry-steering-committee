# per_structure_predictions.csv — schema notes and a bug this surfaced

> **⚠ Correction (2026-10-02, no re-run):** the H12-local pLDDT values
> reported throughout `H3_RESULTS_REPORT.md` (87.50 / 77.52 / 83.19 / 77.00,
> and the C4–C6 extension) are computed over construct residues
> **480–487**, not the 479–486 stated there — an off-by-one bug, detailed
> below. The values are correct as measured; the documented window was off
> by one. **The arm ordering and the p = 0.0079 separation have not been
> re-verified under the corrected window — this is stated as an open
> question, not predicted to hold.**

**No submission schema exists to match.** Checked `rorgt_candidate_kit/README.md`,
`HACKATHON.md`, and everything under `dataset/` — the README says outright
that "the set submissions are scored on is not in [this kit]," and nothing
specifies a format. This CSV is our own, built for judges to recompute a
ROC curve or move our decision threshold themselves rather than trust our
own aggregate numbers.

## Coverage and reconciliation

1,943 rows: every H1 (51), H2 (200), H2-post-hoc (350), H3 core+extended
(350), and ChEMBL-pool (942, including 12 rows for compounds attempted but
never successfully predicted) structure, plus H4's D1 arm (50) — the one
round that postdates the official rescore and had to be scored fresh here
(`lib.h12_state.call_states`, same official reference, no new GPU job).
D2 is **not** duplicated as a separate H4 round: it's the same 50
structures already counted under `h2_posthoc`'s `coactivator` arm, scored
once.

Reconciled against every previously published aggregate, computed
independently from this CSV alone:

| check | from this CSV | published |
|---|---|---|
| rows in the historical "1,191 excl. pool batches" grouping (h1+h2+h2_posthoc+h3+pool120) | 1,191 | 1,191 |
| decisive (\|margin\|≥1Å) antagonist, any fold quality, in that grouping | 209 | 209 |
| decisive antagonist restricted to well-folded (pLDDT≥70) | 170 | 170 |
| well-folded structures in that grouping | 990 | 990 |
| pooled ChEMBL per-compound accuracy (258 compounds actually scored) | 178/258 | 178/258 |

All five match exactly. One historical quirk worth knowing if you go
looking for it: the original "1,191 excl. pool batches" label is slightly
inaccurate — it already included pool120 (batch A, 240 structures) inside
its own "H3 (apo/C1–C6)" bucket, and only truly excluded rand150 (batch B)
and D1 (which didn't exist yet). Reproduced here as-is, not re-labeled,
since the numbers it produces are the ones already on the record.

## A bug this surfaced: h12_state_classifier.py's H12-local pLDDT is off by one residue

Building `h12_local_plddt` for this CSV required reproducing the
token-index arithmetic `h12_state_classifier.py` already uses for H3.
Checked it directly against ground truth (`gemmi`'s own residue/atom
layout for a known structure) rather than trusting it, and it doesn't
hold up: that script computes `h12_toks = [r - offset for r in REGION]`
and indexes the per-residue pLDDT dict with those values directly, but
`atom_tokens` in `*_scores.json` is **0-based** (token 0 = predicted
residue 1), while `r - offset` produces a **1-based** residue number. The
result is indexed one token too high throughout — every "H12-local
pLDDT" number anywhere in this project's H3 reporting (87.50 / 77.52 /
83.19 / 77.00 and the C4–C6 extension) is actually the mean over
construct residues **480–487**, not the claimed **479–486**.

This CSV's own `h12_local_plddt` column does **not** have this bug — it's
computed fresh here with the corrected `token = residue - offset - 1`,
verified against `gemmi`'s atom/residue order directly (see
`build_per_structure_predictions.py`'s `h12_local_plddt()` docstring).

**Not fixed in `h12_state_classifier.py` itself in this pass** — that
would silently change every number in `H3_RESULTS_REPORT.md` and the
classifier's own calibration table without a decision to re-report them.
Per instruction: a correction note, not a re-run. The arm ordering
(apo > C2 > C1 ≈ C3) and the p = 0.0079 separation have **not** been
re-verified under the corrected (480–487) window — this is an open
question, not something predicted to hold either way.

## Provenance

`build_per_structure_predictions.py`, run 2026-10-02. Reads
`~/h3/official_state_calls_all.csv` (already on disk) and each
structure's own `*_scores.json` for pLDDT; no new predictions. D1's
official-frame rescore is the one piece of new scoring this script runs,
using the already-pinned `rorgt_candidate_kit` reference pair.
