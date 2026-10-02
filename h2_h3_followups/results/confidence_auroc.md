# Confidence-signal AUROC — margin vs. pLDDT

**Scope:** existing 258-compound pool, no new predictions, no classifier fit.

n = 258 compounds (0 dropped: missing margin or pLDDT).

## AUROC for predicting call correctness

| signal | AUROC | 95% CI (bootstrap, n=2000) |
|---|---|---|
| decisive margin | 0.681 | [0.606, 0.756] |
| mean pLDDT | 0.631 | [0.547, 0.709] |

Margin's point estimate is higher, but **the CIs overlap substantially**
(margin's lower bound, 0.606, sits well inside pLDDT's interval). This
does not support a claim that margin is a significantly better correctness
signal than pLDDT — only that its point estimate is higher. Both clear
0.5 (better than chance) with CIs that exclude 0.5.

## Margin vs. pLDDT correlation

Pearson r = **0.446** across the 258 compounds — a moderate positive
correlation, not near-independent. The two signals share real information
(plausibly both reflecting how well the H12 region folded at all) rather
than contributing fully separate evidence. Reporting both on a slide as
if they were independent confirmations of the same conclusion would
overstate the evidence; r = 0.45 means they partially double-count.

## Reading this against Task 3

Task 3's quartile table (non-monotonic for pLDDT, clean for margin) and
this AUROC comparison are answering different questions and don't have to
agree — and here they only partially agree. Margin's AUROC point estimate
is higher (0.681 vs 0.631), consistent with Task 3's framing, but the gap
is **not statistically clean** given the CI overlap — this is weaker
support for "margin is the better signal" than Task 3's summary alone
implied. Both signals are informative; margin is not decisively better
than pLDDT by this test.
