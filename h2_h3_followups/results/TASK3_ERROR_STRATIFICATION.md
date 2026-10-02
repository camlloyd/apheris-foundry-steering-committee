# Task 3 — error stratification: results

**Scope:** rescore/re-analysis of the existing 258-compound pool. No new
classifier fit, no touch to `h12_state_classifier.py`, no new predictions.
Per-compound data: `task3_per_compound.csv`.

Overall: **80/258 (31.0%) misclassified** — 51 antagonist, 29 agonist (of
130 and 128 respectively: 39.2% vs 22.7% error rate), matching the
77.3%/60.8% accuracy split already reported.

## Binding mode (orthosteric vs. allosteric) — not possible with this data

**Checked first, as instructed.** `chembl_subset_120.json` and
`chembl_subset_random150.json` (the only annotation that travels with the
pool) carry `molecule_chembl_id`, `smiles`, `activity_label`, `potency_nm`,
`scaffold` — no orthosteric/allosteric field, and nothing elsewhere in the
kit maps ChEMBL IDs to a binding-mode label. The brief calls this split
"half the problem," and this project still has not addressed it — not
because it was skipped, but because the data to do it with real labels
doesn't exist in what we have. Falling back to scaffold class, as the
task instructions anticipate.

## By scaffold

Scaffolds with ≥3 compounds in the pool (most scaffolds have 1–2; 207 of
217 distinct scaffolds don't meet this bar, so read all of this as
suggestive, not powered):

| scaffold (truncated) | n | errors | error rate |
|---|---|---|---|
| `O=S(=O)(c1ccccc1)N(Cc1ccc(N2CCNCC2)cc1)C1CCC1` | 3 | 3 | **100%** |
| `O=S(=O)(c1ccccc1)N1CCOc2ccccc21` | 5 | 3 | 60% |
| `O=S(=O)(c1ccccc1)N1C[C@H](CC2CCOCC2)Oc2ccc(-c3ccccc3)cc21` | 3 | 1 | 33% |
| `O=C(NCc1ccccc1)c1ccc2[nH]c3ccccc3c2c1` | 3 | 1 | 33% |
| `O=S(=O)(NCc1ccccc1)N1CCc2cc(N3CCNCC3)ccc2C1` | 4 | 1 | 25% |
| `O=C(Cc1ccccc1)Nc1ccc(-c2ccccc2)cc1` | 5 | 0 | 0% |
| `O=C(Cc1ccccc1)Nc1ccc(CCN2CCCc3ccccc32)cc1` | 6 | 0 | 0% |
| `O=C(Cc1ccccc1)Nc1ccc(CNCC2CCCCC2)cc1` | 3 | 0 | 0% |
| `O=C1C=CC2c3nc(-c4ccnc5ccccc45)nc(-c4ccccc4)c3CC[C@@H]2C1` | 4 | 0 | 0% |

**Pattern worth a slide, flagged with its n:** every 100%/60%/33%-error
scaffold above is a **benzenesulfonamide** core (`O=S(=O)(c1ccccc1)N...`).
Every 0%-error scaffold above is a **phenylacetamide** core
(`O=C(Cc1ccccc1)Nc1ccc...`) or one other rigid polycyclic scaffold. That's
a real, chemically coherent split, not noise scattered across unrelated
cores — but it rests on 3–6 compounds per scaffold, so report it as "the
errors cluster by chemotype, sulfonamides look like the harder class in
this pool" rather than a validated rule.

## By mean pLDDT (quartile)

| quartile | pLDDT range | n | errors | error rate |
|---|---|---|---|---|
| Q1 (lowest) | 88.4–90.4 | 64 | 32 | 50.0% |
| Q2 | 90.4–91.9 | 65 | 22 | 33.8% |
| Q3 | 92.0–93.4 | 64 | 9 | **14.1%** |
| Q4 (highest) | 93.4–95.0 | 65 | 17 | 26.2% |

Mostly the expected direction (lower confidence → more errors), but
**not monotonic** — Q4 (highest confidence) has a higher error rate than
Q3. Reported as measured, not smoothed into a clean story it isn't: pLDDT
alone is not a reliable, monotonic error predictor across the whole range
in this pool (consistent with H3's own finding that pLDDT tracks
pharmacological class, not a clean confidence-correctness relationship in
general).

## By decisive margin

| margin band | n | errors | error rate |
|---|---|---|---|
| [0, 0.5) | 79 | 36 | 45.6% |
| [0.5, 1.0) | 60 | 24 | 40.0% |
| [1.0, 2.0) | 73 | 13 | 17.8% |
| [2.0, ∞) | 46 | 7 | **15.2%** |

This one **is** clean and monotonic — consistent with Task 2's coverage
curve (same signal, re-cut by band instead of swept threshold). Margin is
a better-behaved error predictor than raw pLDDT in this data.

## Headline for this task

The orthosteric/allosteric split the brief calls out explicitly could not
be tested — no ground-truth label for it exists in what we have, and that
gap is worth stating on a slide rather than silently working around.
Where we *can* stratify: antagonist compounds error at nearly 2× the rate
of agonist ones (39.2% vs 22.7%), the hardest-performing scaffolds in the
pool share a benzenesulfonamide core (small n, flagged as such), margin is
a clean monotonic error predictor, and raw pLDDT is not.

## Provenance

`task3_error_stratification.py`, run 2026-10-02, reading
`~/h3/official_state_calls_all.csv`, `~/h3/chembl_subset_120.json`,
`~/h3/chembl_subset_random150.json`, and each structure's already-written
`*_scores.json` for per-atom pLDDT. No new GPU jobs.
