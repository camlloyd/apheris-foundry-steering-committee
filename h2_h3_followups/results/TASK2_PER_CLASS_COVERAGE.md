# Per-class accuracy/coverage curves — agonist vs. antagonist

**Scope:** existing 258-compound pool (128 agonist, 130 antagonist), no new
predictions, no threshold tuning. Full curve: `task2_per_class_coverage.csv`
(margin threshold 0.0–4.0 Å, step 0.1).

## Where each class crosses 80% accuracy

- **Agonist crosses 80% almost immediately** — at margin ≥ 0.1, with
  **100% coverage retained**. Agonist accuracy is already close to 80% at
  full coverage (77.3%) and clears it with almost no calls withheld.
- **Antagonist never reaches 80% accuracy anywhere in the swept range
  (0.0–4.0 Å).** Its best point in the whole sweep is **74.3% at
  margin ≥ 1.3** (57% coverage, n=74) — close, but it does not clear the
  bar, and accuracy does not climb steadily toward it.

## This is weaker than the pooled curve in `TASK2_COVERAGE_RESULTS.md` implied

That report showed a clean, close-to-monotonic pooled accuracy/coverage
curve and read it as "a genuinely calibrated confidence signal." Split by
class, that's not the full story:

| margin ≥ | agonist coverage | agonist acc | antagonist coverage | antagonist acc |
|---|---|---|---|---|
| 0.0 | 100.0% | 77.3% | 100.0% | 60.8% |
| 0.5 | 89.1% | 93.9% | 86.2% | 50.0% |
| 1.0 | 70.3% | 91.1% | 64.6% | 66.7% |
| 1.5 | 57.0% | 91.8% | 52.3% | 73.5% |
| 2.0 | 52.3% | 92.5% | 36.2% | 65.0% |
| 2.5 | 52.3% | 92.5% | 30.0% | 61.5% |
| 3.0 | 48.4% | 93.6% | 23.1% | 50.0% |
| 3.5 | 48.4% | 93.6% | 20.0% | 42.3% |

**Agonist accuracy climbs cleanly and almost monotonically with margin, from
77.3% to the low-90s.** **Antagonist accuracy does not** — it rises from
60.8% to a peak of ~74% around margin 1.3–1.5, then *degrades* as margin
tightens further (down to 42.3% by margin ≥ 3.5, though n shrinks to the
low dozens there and the noise grows accordingly). The pooled curve's
apparent monotonic improvement in `TASK2_COVERAGE_RESULTS.md` is driven
almost entirely by the agonist half of the pool; it was not showing a
property that holds for both classes, and presenting it unqualified would
overstate what abstention buys you on the harder class.

## Does abstention rescue antagonists more than agonists?

**No — the opposite.** Agonist reaches a high-accuracy, high-confidence
regime quickly and at low cost (80%+ accuracy at 100% coverage already);
antagonist needs to give up roughly 40%+ of its coverage just to approach
(not reach) 80%, and then degrades again if pushed further. Quantified:
to hit the closest antagonist gets to 80% (74.3%), the pool gives up 43%
of antagonist coverage; agonist needs to give up nothing. **Abstention
is a tool for improving agonist confidence, not for rescuing antagonist
calls** — on this pool, withholding low-margin calls does not solve the
antagonist-accuracy problem, it just shrinks the antagonist denominator
without reliably raising the numeral past where it already was.

## Honest summary

This sharpens, and partly undercuts, the Task 2/Task 3 framing: margin is
a real signal, but mostly for agonist calls. For antagonists — the harder
class this whole project keeps finding — margin does not deliver a clean
path to high confidence; it has a soft ceiling around 65–74% that more
selective thresholds don't clearly beat, within the noise this sample size
allows. Don't present the pooled coverage curve on a slide without this
caveat.

## Provenance

`task2_per_class_coverage.py`, run 2026-10-02, reading only
`~/h3/official_state_calls_all.csv` (already on disk, no new predictions).
