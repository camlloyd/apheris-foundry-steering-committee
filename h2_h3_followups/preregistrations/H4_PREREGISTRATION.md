# H4 pre-registration — coactivator placement, agonist vs inverse agonist

**Timestamped before launch:** 2026-10-02, overnight run, prior to any D1/D2 scoring.

**H4.** The model has learned H12 conformation and coactivator placement
independently rather than as coupled. Prediction: peptide placement in D2
does not differ from D1 by more than seed noise.

**Design.** Full RORγ LBD construct (265–507) + SRC2-2 LXXLL coactivator
peptide (3KYT chain C, residues 686–697, `KHKILHRLLQDS`) + ligand, 5 seeds
(42–46) × 10 diffusion samples, full MSA (`h2/msa_full`, receptor chain
only — same setup already proven to run for D2).

| arm | ligand | purpose | status |
|---|---|---|---|
| D1 | HC2 (25-HC, agonist) | positive control — does the model place the peptide at all? | **launching tonight** |
| D2 | 4P1 (inverse agonist) | the question | **already exists**: `h2_posthoc/out/coactivator/` (50 structures, run earlier as a post-hoc H2 probe, same protocol — 5 seeds × 10 samples, same MSA, same peptide sequence, only the ligand differs) |

**Primary readout (pre-specified):** peptide heavy-atom RMSD to its position
in 3KYT chain C, after superposing on the LBD (chain A) only — not a joint
fit, which would let the peptide's placement inflate the receptor-fit score
or vice versa. Compare D1 vs D2 means, seed-level, against the D1 seed
spread as the noise scale (same exact-permutation-test convention as H2/H3).

**Secondary readouts:** peptide–LBD interface contact recall (heavy atoms
within 4.5 Å, same cutoff as `per_residue_discriminator.py`), and the H12
state call from the official `lib/h12_state.py` (5VB7/6T4I, window 484–507)
on each structure's receptor chain.

**Decision rule, decided in advance:**
- D1 ≈ D2 (RMSD difference within D1's own seed spread) → **H4 supported**:
  the model places the coactivator the same way regardless of H12 state,
  i.e. it has not learned the H12–coactivator coupling as a joint
  constraint — a real, quotable mechanistic finding, consistent with
  everything H1–H3 already showed (the model defaults to one placement
  strategy regardless of ligand).
- D1 clearly ≠ D2 (peptide is displaced/absent specifically when the
  inverse agonist is bound, beyond seed noise) → **H4 rejected**: the model
  *has* captured some version of the coupling, at least as a co-folding
  correlation — the more interesting and more publishable outcome, but not
  the one either prior result (H1/H2's geometric negative, or the official
  rescore's 17% antagonist recovery) predicts.
- Ambiguous, partial, or either arm fails to complete: **reported as such**,
  held out of the final slide, kept as "what we'd check next" per the
  session's own rule (a half-finished D1/D2 contrast is worse on a slide
  than an honest "we started this and ran out of time").

**Rules for how this gets used tonight:** the deck is built as if this
result doesn't exist. No placeholder slide for it. Only if D1 completes
cleanly and the D1/D2 contrast is unambiguous does it earn a slide at the
end, added in the morning — not before.

## Provenance

- Peptide sequence: 3KYT chain C, residues 686–697 (`KHKILHRLLQDS`), the
  exact LXXLL motif this project has referenced throughout (see
  `H3_RESULTS_REPORT.md` appendix).
- D1 query: `h2_posthoc/in_coact_d1/query.json` (new, this session).
- D2 query (already run): `h2_posthoc/in_coact/query.json` (identical
  structure, ligand swapped to 4P1 — this is the file `run_coactivator.sh`
  already used).
- Scoring (to run in the morning, not tonight): `h2_h3_followups/scripts/
  h4_coactivator_score.py` — peptide RMSD, interface contacts, and official
  H12 state call per structure, D1 vs D2.
