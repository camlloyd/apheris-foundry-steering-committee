# RORγ Two-State Steering — Day 1 Findings Handoff

Foundry Ops → Science Lead, end of Day 1 (2026-10-01). Covers everything
generated and scored since `HANDOFF.md`/`RUNBOOK_rorgamma.md`. Source data:
`apheris_kit_rorgamma/runs/`, scored with `score_run.py` (Lily's own,
unmodified).

## Setup, and the one pivot that colors every result below

All predictions in this document ran on **OpenFold3**
(`of3-ob-2025-06-30-174k.pt`), not Boltz-2 — the Hub/`apheris-foundry` CLI
was never reachable from the compute host (see `AGENTS.md`'s "Resolution
(2026-10-01)" section for the full investigation). Every finding below is
therefore a statement about this OpenFold3 checkpoint, not about co-folding
models in general. That distinction matters for how strongly any of this
gets stated on a slide.

Core arms seeded per the pre-registered plan: invago A0/A1/A4 at 5 seeds,
25hc A0/A1/A2 at 2 seeds; everything else at 1 seed (quota-tight floor).

## Primary result: 0% H12-out recovery, invago direction — now comprehensively tested

Every arm, every seed, every lever tried landed H12-in. That now includes,
after today's additional runs:

| Arm | Lever | Seeds | State recovery | Ligand dist (Å) | Contact recall |
|---|---|---|---|---|---|
| A0_baseline | none | 5 | 0% | 6.62 | 0.41 |
| A1_pocket_specific | pocket | 5 | 0% | 2.65 | 0.65 |
| A2_pocket_shared | shared pocket (control) | 1 | 0% | 1.08 | 0.85 |
| A3_template_target | correct template | 1 | 0% | 12.80 | 0.10 |
| A4_template_wrong | wrong template (control) | 5 | 0% | 1.30 | 0.68 |
| A5_msa_free | no MSA | 1 | 0% | 3.78 | 0.25 |
| A6_pocket_plus_template | pocket + correct template | 1 | 0% | 3.47 | 0.60 |
| A7_pocket_msa_free | pocket, no MSA | 1 | 0% | 8.60 | 0.15 |
| A8_conflict | pocket + wrong template | 1 | 0% | 3.56 | 0.65 |
| **A10_msa_depth32** (new) | partial MSA (depth 32), no conditioning | 5 | 0% | 1.30 | 0.78 |
| **A11_msa_depth32_pocket** (new) | partial MSA (depth 32) + pocket | 5 | 0% | 1.62 | 0.75 |
| B0-B3 (OpenFold3 "mirrors") | — | 1 each | 0% | 1.5-12.8 | 0.10-0.75 |

Not "steering fails to beat baseline" — nothing in the design space, run
alone or combined, ever reaches the target state. Per the pre-registered
criterion, this is reported as a negative result, not a failure.

## What was ruled out, and how

1. **Template conditioning was genuinely dispatched, not dropped.** A3
   (handed the correct H12-out template, 4ZJW) is statistically
   indistinguishable from A4, the deliberately-wrong-template control
   (0.807 Å apart, same verdict, 0.08 Å margin difference). Verified
   directly against the actual Nextflow task config: `template_paths` was
   correctly set. The signal is being overridden, not silently ignored.
2. **MSA-free isn't a clean "prior removed" test.** A5/A7/B2 have mean
   pLDDT 43.5-43.7 — AlphaFold's own "Very Low confidence" band (<50).
   Full MSA removal breaks the model rather than cleanly ablating the
   evolutionary prior.
3. **Partial MSA depth reduction, alone, doesn't move the state either.**
   A10 (random subsample to depth 32 of the 2020-sequence ColabFold
   alignment, no conditioning): 0/5.
4. **Partial MSA + pocket conditioning together, the interaction-effect
   hypothesis, is also ruled out.** A11 (depth 32 + the same pocket
   contacts as A1): 0/5. A weakened prior didn't let the structural lever
   finally compete.

The lever space actually explored is now large: pocket, template,
wrong-template, pocket+template conflict, full MSA, zero MSA, partial MSA,
and partial MSA + pocket. All eleven land H12-in. That is a comprehensive,
not a preliminary, negative result on this checkpoint.

## Secondary finding: pose accuracy and state recovery decouple cleanly

A10/A11 (depth-32 MSA) have ligand centroid distances of 1.30/1.62 Å vs.
A0's (full MSA) 6.62 Å — dramatically better pose — while the state call is
identical (H12-in, every seed) across all three. MSA depth measurably
changes *where* the ligand sits without changing *which* conformational
state the protein adopts. Same lesson KinConfBench established on kinases
(pose ≠ state), now independently observed on a nuclear receptor via a
different lever (MSA depth, not just templating).

## 25-HC direction: the reported "100%" needs a quality floor before it's quotable

`score_run.py`'s margin test is purely relative (no absolute quality
floor), so three broken MSA-free folds in that direction (A5/A7/B2, same
pLDDT 43.5-43.7 problem as above) were counted as valid "H12-in, decisive"
successes. **Fix before this goes on a slide:** require pLDDT ≥ 70
(AlphaFold's "Confident" band) and a region-RMSD floor before accepting a
state call. Doesn't touch the invago 0%, but inflates the apparent 25hc
success rate as currently computed.

## Bonus finding: the challenge's own framing, demonstrated on their target

A0's ligand placement: 0.34-0.37 Å in the agonist (well-represented)
direction vs. 0.93-12.91 Å in the inverse-agonist (under-represented)
direction. The model places the common ligand precisely and the rare one
unreliably — the hackathon's own "under-represented conformations" problem,
shown on RORγ itself, no extra experiment required.

## Judging-rubric metrics (organizer-provided)

Rubric: state-call correctness (H12 RMSD) + mean(GDT-HA, lDDT-PLI, Ligand
BiSyRMSD). `apheris-data`'s `compute-metrics` Foundry task provides GDT-HA
directly (verified live) but no field under either of the other two names.
Implemented both from scratch, TDD'd, 100% test coverage:
- `foundry-ops/rorgamma_metrics.py::lddt_pli` — pure distance-geometry
  (protein-ligand interface lDDT), no chemistry library needed.
- `foundry-ops/rorgamma_ligand_rmsd.py::ligand_bisy_rmsd` — wraps RDKit's
  symmetry-aware `GetBestRMS` rather than reimplementing automorphism
  matching.

Not yet wired into an end-to-end scoring script against our actual
predictions — the functions are built and tested, integration is
follow-up work.

## Suggested pitch framing

"Pocket and template conditioning cannot move this co-folding model's
inverse-agonist state call — and neither can removing or thinning the MSA,
alone or in combination with structural conditioning. But pose accuracy and
state recovery decouple cleanly (MSA depth fixes ligand placement without
touching the state), and the model is specifically unreliable on the
under-represented ligand class — both demonstrated on a real therapeutic
target, both honest, neither requiring a steering claim we can't support."

## What's left, if there's time

- Wire the two new metrics into an actual scoring pass against the full
  prediction set (judging-rubric numbers, not just Lily's internal scorer).
- Add the pLDDT/RMSD quality floor to the 25hc numbers before they're
  quoted anywhere.
- Genuinely different levers not yet tried: sequence-similarity MSA
  *clustering* (not random subsampling) in case a conformation-linked
  sub-family exists in this alignment; actual Boltz-2, if Hub access ever
  resolves.
