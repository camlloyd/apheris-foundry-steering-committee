# H3 — pre-registration

**Written 2026-10-01 14:56 UTC. No H3 prediction has been run or seen at the
time of writing.** H1 and H2 are retired and reported on their own evidence
(see `H2_PREREGISTRATION.md`, `titration_summary.md`, `~/h2_posthoc/`). H3
below is a new hypothesis for new data.

---

## Motivation

H1 and H2 together show: whenever the model has enough MSA depth to fold the
RORγ LBD at all (depths 128/512, the two post-hoc fold-stability sweeps, every
template condition including a genuinely alternative-pose template), it folds
H12 into the agonist-locked pose. 0/100+ structures across every conditioning
lever tried reached H12-out by the RMSD state-call criterion. The geometric
state call is a strict pass/fail gate (pLDDT ≥ 70 AND decisive RMSD margin) and
may be discarding softer, continuous signal that a binary gate can't see.

H3 asks a different, more tractable question: is the model's **local
confidence at the state-defining residues** (H12, 479–486) sensitive to ligand
identity at all — independent of whether it ever crosses the geometric
threshold into a different RMSD basin.

## Primary claim (registered; this is what H3 can actually establish)

**H12 local pLDDT (mean per-residue pLDDT over region 479–486) differs between
ligand conditions.** Specifically: the inverse-agonist arm (C2) and the
agonist-class arms (apo, C1) show significantly different H12-local pLDDT
distributions, assessed seed-level (mean over diffusion samples per seed, 5
seeds per arm), by exact permutation test on the 5-vs-5 seed means, reporting
the exact p-value and an effect size (standardized mean difference /
Cliff's delta). Complete separation of 5 vs 5 seed-level means corresponds to
p ≈ 0.008 under the permutation null (2/C(10,5)) — strict, and treated as
such.

**This is deliberately narrower than "H12 confidence is ordered opposite to
crystallographic stability."** That stronger claim is an *interpretation*,
not the registered primary claim, for the reason given below.

### Why the stronger claim is demoted to interpretation, not registered

25-hydroxycholesterol (the agonist in 3KYT, arm C1) is a small, flexible
sterol. 4P1 (the inverse agonist, arm C2) is a larger, rigid synthetic ligand
that fills the pocket more completely. H12 packs directly over the pocket, so
the model's H12 confidence may simply track **pocket occupancy / ligand pose
certainty**, not pharmacological state. That confound predicts the *same*
C2 > C1 ordering as a genuine state story, so this design alone cannot
distinguish them. It additionally predicts a **monotonic** apo < C1 < C2
pattern (more/stiffer ligand → more confidence, regardless of direction),
which the state story does not predict. AF3-class confidence metrics are
separately documented as unreliable specifically in ligand-binding contexts,
which is independent prior reason to suspect the pocket-confound reading over
a state-encoding one.

**Decided now, before any data exists:**
- If the result is **apo ≈ C1 < C2** (inverse agonist lower, agonist-class
  arms similar to each other): read this as the state-dependent pattern —
  interesting, but not fully diagnostic without arm C3 (below).
- If the result is **monotonic apo < C1 < C2**: read this as the
  pocket-filledness confound, explicitly, in the writeup — not as evidence for
  state encoding.
- **Apo ambiguity, decided now:** apo RORγt crystallizes with H12 ordered
  (active-like) despite the empty pocket. A drop in apo H12-local pLDDT is
  therefore compatible with *both* stories (no ligand to stabilize the pocket,
  or no ligand to signal "agonist") and will be reported as ambiguous on its
  own, not used alone to argue for either reading.

## Arms

| arm | ligand | pharmacology | rigidity | why |
|---|---|---|---|---|
| apo | none | — | — | missing baseline H1/H2 never ran cleanly (A7 became a pocket+MSA arm, confounding ligand presence with conditioning) |
| C1 | HC2 (20/25-hydroxycholesterol) | agonist | flexible, small sterol | the natural agonist ligand, already in 3KYT |
| C2 | 4P1 | inverse agonist | rigid, pocket-filling synthetic | the inverse agonist used throughout H1/H2 |
| C3 (if run) | 6F1 "BIO592" (PDB 5IZ0/5IXK context; SMILES from RCSB CCD: `CCN1c2ccc(cc2OCC1=O)N(CC(F)(F)F)S(=O)(=O)c3ccc(c(c3)C)C`) | **agonist** (confirmed by FRET coactivator-recruitment assay in Fauber et al.) | rigid synthetic, non-sterol | breaks the sterol-flexibility confound: if H12 pLDDT tracks pharmacology (agonist), C3 should pattern with C1 despite being rigid like C2; if it tracks rigidity/pocket-filling, C3 should pattern with C2 despite being an agonist |

Construct, numbering and MSA: identical full LBD (3KYT numbering 265–507) used
throughout H1/H2. **Full MSA depth** (the ~2020-sequence a3m already fetched
at `~/h2/msa_full/`, not a truncated ladder) for every arm — this hypothesis is
about ligand identity, not MSA depth, so every arm uses the condition that
folds cleanly. No template conditioning (template lever already tested
separately in H2's post-hoc template×depth arm). Same 5 seeds (42–46), same 10
diffusion samples per seed, as every other arm in this project — 50 structures
per arm, 150 (or 200 with C3).

## Secondary readouts (reported beside the primary, never merged into it)

1. **Ligand/pocket-residue mean pLDDT per arm** — tests the pocket-confound
   directly. If pocket-lining residues (not H12) show the same C2 > C1
   ordering as H12 itself, that's evidence the effect is about pocket
   occupancy generally, not specifically about the state-defining residues.
2. **H12 region (479–486) coordinate RMSD distributions to both 3KYT and
   4ZJW**, per arm — reported as distributions, not collapsed into the binary
   state call. If the inverse-agonist arm shows a consistent 2–3 Å shift
   without crossing the pre-registered 1 Å decisive margin, that is coordinate
   signal the geometric H2 criterion was built to discard, and it will be
   reported as such — a secondary, continuous readout, not folded into "H2
   SUPPORTED/NOT SUPPORTED."

## Unit of analysis (decided now)

**Seed-level mean over the 10 diffusion samples**, not pooled structures. 5
seeds per arm → 5 numbers per arm → permutation test on those. Pooling all 50
structures per arm as if independent would inflate n and understate the real
degrees of freedom (samples within a seed share the same MSA/seed-level
stochasticity); this is stated here so it cannot be chosen after seeing which
version gives a better p-value.

## Threats to validity, named now

1. **Pocket-occupancy confound** (above) — the central threat; C3 is the
   mitigation, not a guarantee. If C3 cannot be run, this is reported as an
   unresolved limitation, not elided.
2. **Apo ambiguity** (above) — decided in advance how it will be read.
3. **n=1 target, n≤4 ligands.** Nothing here generalizes beyond RORγ LBD
   without a cross-target follow-up (explicitly out of scope for this round —
   see HANDOFF_H2.md's novelty note).
4. **pLDDT is a model confidence estimate, not a ground-truth disorder
   measurement.** A pLDDT difference is evidence the model's *estimate of its
   own certainty* is ligand-dependent; it is not, by itself, evidence that the
   biological state is ligand-dependent. This is the whole reason the primary
   claim is framed as "H12 local confidence is ligand-dependent" rather than
   "ligand controls H12 state."

## Declared negative

If no significant separation is found between any arms at the seed level
(permutation p not below a pre-specified α, suggested α = 0.05, exact
threshold a function of n once C3's inclusion is finalized), the conclusion
is: H12 local confidence in this model is not detectably ligand-dependent for
these ligands at this MSA depth — a third negative in the same direction as
H1 and H2, now at the level of the model's own confidence estimate rather than
its geometry.

## Sign-off

Criterion fixed at 2026-10-01 14:56 UTC, before any H3 job was submitted.
Any later change to this file must be recorded as a dated amendment below,
with its reason, and the original criterion left legible.

### Amendments

(none)
