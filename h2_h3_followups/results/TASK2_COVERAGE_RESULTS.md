# Task 2 — accuracy/coverage curve: results

**Scope:** rescore of the existing 258-compound pool, no new predictions.
Full curve: `task2_accuracy_coverage_curve.csv` (margin threshold 0.0–4.0 Å,
step 0.1).

## Does the evaluation harness permit abstention?

**Checked first, as instructed — the answer is: we don't know, and the
materials provided don't say.** `rorgt_candidate_kit/README.md` states
explicitly that "the set submissions are scored on is not in [this kit]";
nothing in `HACKATHON.md`, the kit's `README.md`, or any file under
`rorgt_candidate_kit/` documents a submission format or an abstention
policy for the held-out evaluation. `lib/h12_state.py`'s own `UNSCORED`
state exists only for a structural reason (an empty H12 crop — wrong
offset/chain), not as a designed "below-confidence, withhold the call"
option. We are not assuming abstention is free; we are stating plainly
that this is unknown rather than picking the assumption that favors the
better number.

Given that, **the headline number reported anywhere in this project
remains the full-coverage figure (69.0% pooled)** — the curve below
characterizes a predictor that *could* abstain, it does not substitute a
better number for one that can't.

## The curve

| margin ≥ | coverage | n scored | accuracy on scored |
|---|---|---|---|
| 0.0 (full coverage, current headline) | 100.0% | 258 | 69.0% |
| 0.8 | 76.7% | 198 | 74.8% |
| 1.0 (this project's existing "decisive" cut) | 67.4% | 174 | 79.3% |
| 2.0 | 44.2% | 114 | 81.6% |
| 2.9 | 37.2% | 96 | 80.2% |

(174 scored at margin ≥ 1.0, not 258 − 12 = 246; the "12 excluded" in
earlier reports is a separate, documented RDKit featurizer failure during
compound selection, upstream of this margin-based withholding — the two
exclusion mechanisms are not the same thing and shouldn't be conflated on
a slide.)

Accuracy rises close to monotonically as coverage drops — this is a
genuinely calibrated confidence signal (withheld calls are disproportionately
wrong ones), not noise. It plateaus around 80–82% in the 2.0–2.9 Å range
rather than continuing to climb toward 100%, which is itself informative:
even the most decisive-margin calls this model makes top out around
80% accuracy, they don't approach certainty.

## How to use this

This is a characterization of a calibrated predictor operating at a
chosen operating point, not a free upgrade from 69.0% to 79.3% or 81.6%.
Quoting a higher number without the coverage it costs would misrepresent
what was measured. If the (unknown) held-out harness does allow
abstention, this curve is the right thing to hand judges to justify an
operating point; if it scores abstentions as wrong, the right operating
point is margin ≥ 0 (the full-coverage number) and this curve is
characterization only, not a quoted accuracy.

## Provenance

`task2_coverage_curve.py`, run 2026-10-02, reading only
`~/h3/official_state_calls_all.csv` (already on disk, no new predictions
or GPU jobs).
