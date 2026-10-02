# H4 — Results Report

**Hypothesis registered:** `../preregistrations/H4_PREREGISTRATION.md`, timestamped
2026-10-02, overnight, before D1 launched and before any scoring ran.
**This report written:** 2026-10-02, after both arms completed and the
pre-registered scoring script was run exactly as specified (no changes made
to the script between registration and running it).

---

## 1. Background

H1/H2 found the model never reaches H12-out on this target (geometric
criterion). H3 found the opposite ligand-dependence in *confidence*
(pLDDT) than the naive story predicts. H4 asks a third, independent
question: when the model is given the LXXLL coactivator peptide (3KYT
chain C, `KHKILHRLLQDS`) to co-fold alongside the receptor, does *where it
places that peptide* depend on which ligand is bound — agonist (D1, HC2)
vs inverse agonist (D2, 4P1, reused from an earlier H2 post-hoc run,
identical protocol)?

## 2. Pre-registered design (unchanged from registration)

5 seeds (42–46) × 10 diffusion samples = 50 structures/arm, 100 total.
Primary readout: peptide heavy-atom RMSD to its position in 3KYT chain C,
after superposing on the LBD (chain A) only. Decision rule: D1 ≈ D2 (within
D1's own seed spread) → H4 supported; D1 clearly ≠ D2 beyond seed noise →
H4 rejected, with the pre-specified expectation that rejection would look
like "peptide displaced/absent specifically when the inverse agonist is
bound."

## 3. Results

### 3.1 Primary: peptide RMSD to 3KYT (receptor-only superposition)

| arm | seed means (42,43,44,45,46) | mean ± sd (seed-level) | all-structure range |
|---|---|---|---|
| D1 (agonist, HC2) | 34.78, 35.87, 41.49, 37.53, 35.66 | **37.07 ± 2.67 Å** | 19.2–55.4 Å |
| D2 (inverse agonist, 4P1) | 4.22, 4.36, 4.45, 4.23, 4.26 | **4.30 ± 0.10 Å** | 3.5–5.4 Å |

Exact permutation test (5-vs-5 seed-level means): observed diff
D1 − D2 = **+32.76 Å**, **p = 0.0079** (complete separation, the tightest
this design can produce — same convention as H2/H3), Cohen's d = 17.4.

### 3.2 Secondary: peptide–LBD interface contact count

| arm | seed means | mean |
|---|---|---|
| D1 | 8.8, 8.3, 9.0, 9.1, 8.5 | 8.74 |
| D2 | 10.3, 10.3, 10.2, 10.4, 10.3 | 10.30 |

p = 0.0079 (complete separation), diff D1 − D2 = −1.56. D2 structures have
*more* peptide-receptor contacts, not fewer.

### 3.3 Secondary: official H12 state call (organizers' 5VB7/6T4I reference, window 484–507)

| arm | agonist | antagonist |
|---|---|---|
| D1 | 48/50 | 2/50 |
| D2 | 50/50 | 0/50 |

Both arms call agonist H12 state almost unanimously — consistent with
H1/H2's finding that this model defaults to H12-in regardless of ligand,
independent of the coactivator-peptide question this round asks.

## 4. Applying the pre-registered decision rule

D1 completed cleanly (50/50) and the D1/D2 contrast is unambiguous on all
three readouts, by a wide, non-overlapping margin (p = 0.0079 on both
peptide RMSD and contact count). Per the rule as written, this is
**H4 rejected**: the model does *not* place the coactivator peptide the
same way regardless of ligand.

**But the direction is the opposite of the one the pre-registration
anticipated.** The decision rule's rejection branch was written expecting
"peptide displaced/absent specifically when the inverse agonist is
bound." What the data show is the reverse: in **D2 (4P1, inverse
agonist)**, the peptide lands within ~4.3 Å of its exact crystallographic
position in 3KYT — itself an *agonist*-bound structure — essentially every
seed, with almost no spread (sd = 0.10 Å across seed means). In **D1
(HC2, agonist)**, the same peptide is folded into a completely different
location relative to the receptor, 19–55 Å from the 3KYT placement, with
large seed-to-seed spread (sd = 2.67 Å on seed means, individual
structures ranging over 36 Å) — and still makes *more* receptor contacts
in D2 than D1 (10.3 vs 8.7), so this is not a case of the peptide simply
floating free in D1 versus bound in D2: both arms bind the peptide
somewhere, only D2's binding mode coincides with the crystal pose.

Stated plainly: giving the model the inverse-agonist ligand reproduced the
coactivator's exact crystallographic pose far better than giving it the
agonist did — the opposite of the pharmacologically expected direction,
and opposite of what the pre-registration's "H4 rejected" scenario
described.

## 5. Interpretation (flagged as interpretation, not a registered claim)

No template conditioning was used in either arm (`template_paths: []` in
both configs) and both share identical MSA/weights/seed settings, so this
is not template leakage. The most defensible reading is that **3KYT's own
coordinates are close to an attractor in the model's co-folding output for
this receptor+peptide combination**, and which ligand is present in the
query changes how strongly that attractor is expressed — not that the
model is reasoning about agonist/inverse-agonist pharmacology correctly.
Given H1/H2's finding that the model's H12 geometry itself is
ligand-insensitive (always H12-in) and H3's finding that its *confidence*
is ordered opposite the naive story, a third ligand-dependent effect
pointing the "wrong" way is consistent with a pattern across this whole
project: real, measurable ligand-conditioned behavior exists, but it does
not track pharmacological class in the direction a mechanistic model would
predict.

This interpretation is not pre-registered and should be labeled as such
wherever it is quoted.

## 6. Limitations

- n = 1 target, 2 ligands, no replication of D1 against an independent
  agonist to confirm the direction holds beyond HC2 specifically.
- D2 was not re-run for this round; it is the same 50 structures used as
  an earlier H2 post-hoc coactivator probe, under the identical protocol —
  reused per the pre-registration, not re-generated.
- The large D1 RMSD values (19–55 Å) indicate the peptide is folded far
  from the 3KYT frame, not merely "worse-fit" within it; a visual/contact
  check (§3.2) rules out "peptide unbound/floating" as the explanation.
- As with every other report in this project: pLDDT and RMSD are the
  model's own outputs, not ground truth about real receptor-coactivator
  biology.

## 7. Provenance

- Scoring run: `h2_h3_followups/scripts/h4_coactivator_score.py`, run
  2026-10-02 per the command in the pre-registration, unmodified from the
  version written and smoke-tested the night before (smoke test used
  copies of D2 structures only, never D1, and was deleted before D1
  completed — see prior session log).
- Raw scores: `h2_h3_followups/results/h4_results.csv` (100 rows).
- Permutation test: exact 5-vs-5, same implementation convention as
  `per_residue_discriminator.py` (H2/H3).
