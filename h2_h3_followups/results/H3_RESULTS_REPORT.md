# H3 — Results Report

> **⚠ See `OFFICIAL_RESCORE_FINDINGS.md` in this same folder.** This report's
> pLDDT window (479–486) is this project's own auto-derived choice, not the
> organizers' official H12 window (484–507). The pLDDT separation reported
> below replicates (and is in fact cleaner) on the official window — see
> that file §5 — but report numbers from here with that context, not as
> the organizers' own definition.

**Hypothesis registered:** `H3_PREREGISTRATION.md`, timestamped 2026-10-01 14:56 UTC, before any H3 job ran.
**This report written:** 2026-10-01, after all 4 arms completed and scored.

---

## 1. Background

H1 (pocket/template steering of RORγ LBD toward the H12-out inverse-agonist
state) and H2 (MSA-depth titration of the same question) were both retired as
negatives: across 122+ structures and nine conditioning levers, **0** ever
reached H12-out by the pre-registered geometric criterion (mean pLDDT ≥ 70
AND region-RMSD margin ≥ 1 Å to the correct reference). Whenever the model
folds the domain at all, it defaults to the agonist-locked H12-in pose.

H3 asks a narrower, more tractable question: is the model's **local
confidence at the state-defining residues** (H12, residues 479–486, 3KYT
numbering) sensitive to ligand identity at all — independent of whether it
ever crosses the geometric threshold into a different conformation.

## 2. Pre-registered design

| arm | ligand | pharmacology | rigidity |
|---|---|---|---|
| apo | none | — | — |
| C1 | HC2 (20/25-hydroxycholesterol) | agonist | flexible, small sterol |
| C2 | 4P1 | inverse agonist | rigid, pocket-filling synthetic |
| C3 | 6F1 "BIO592" (PDB 5IZ0/5IXK; SMILES from RCSB CCD: `CCN1c2ccc(cc2OCC1=O)N(CC(F)(F)F)S(=O)(=O)c3ccc(c(c3)C)C`) | **agonist** (confirmed by FRET coactivator-recruitment assay, Fauber et al.) | rigid synthetic, non-sterol |

- Construct: full LBD, 3KYT numbering 265–507, identical to every H1/H2 arm.
- MSA: full depth (~2020 sequences), not a truncated ladder — H3 is about ligand identity, not MSA depth.
- No template conditioning.
- Seeds 42–46, 10 diffusion samples/seed → 50 structures/arm, 200 total.
- **Primary claim (registered):** H12-local pLDDT (mean per-residue pLDDT over 479–486) differs between ligand conditions, tested as an exact permutation test on 5-vs-5 **seed-level means** (not pooled structures — decided in advance).
- **Explicitly demoted to interpretation, not the registered claim:** "H12 confidence is ordered opposite to crystallographic stability." Reason: 25-HC (C1) is a small flexible sterol, 4P1 (C2) is a rigid pocket-filling synthetic — a pocket-occupancy confound predicts the same C2 > C1 ordering as a genuine state story, so this comparison alone can't distinguish them. C3 was added specifically to break this confound: it is rigid/pocket-filling like C2 but pharmacologically an agonist like C1.
- **Decided in advance:** apo < C1 < C2 (monotonic) → read as pocket-filledness confound. apo ≈ C1 < C2 → read as state-dependent. Apo alone is acknowledged as ambiguous either way (apo crystallizes H12-ordered despite an empty pocket).
- Secondary, non-primary readouts: pocket-residue pLDDT (confound check) and H12 coordinate RMSD distributions (continuous signal the binary H2 criterion discards).

## 3. Results

### 3.1 Primary: H12-local pLDDT, seed-level means

| arm | seed means (42,43,44,45,46) | mean |
|---|---|---|
| apo | 86.9, 87.3, 88.3, 87.9, 87.1 | **87.50** |
| C1 (agonist, flexible) | 75.9, 79.0, 77.5, 78.8, 76.3 | **77.52** |
| C2 (inverse agonist, rigid) | 85.2, 81.2, 81.9, 84.2, 83.4 | **83.19** |
| C3 (agonist, rigid) | 76.8, 78.5, 76.5, 78.6, 74.6 | **77.00** |

| comparison | diff | exact permutation p (5v5) | Cohen's d | Cliff's δ |
|---|---|---|---|---|
| C2 vs C1 | +5.67 | **0.0079** (complete separation) | +3.70 | +1.00 |
| C2 vs apo | −4.31 | **0.0079** | −3.51 | −1.00 |
| **C3 vs C1** | −0.52 | **0.60** (not significant) | −0.34 | −0.20 |
| C3 vs C2 | −6.20 | **0.0079** | −3.78 | −1.00 |

p = 0.0079 is the tightest attainable value for a 5-vs-5 exact permutation test (2/C(10,5)) — i.e. complete, non-overlapping separation of the two seed sets.

**Primary claim: confirmed.** H12-local pLDDT is ligand-dependent, at the maximum statistical confidence this design can produce.

**But the pattern matches neither pre-registered interpretation.** Actual ordering: **apo (87.5) > C2 (83.2) > C1 ≈ C3 (77.0–77.5)**.
- Not the state-dependent pattern (predicted apo ≈ C1 < C2) — apo is the highest, not tied with C1.
- Not the monotonic pocket-confound pattern (predicted apo < C1 < C2) — the order is reversed from that.

**What *is* resolved cleanly: the rigidity/pocket-filling confound is rejected.** C3 (rigid, pocket-filling, like C2) patterns statistically with C1 (flexible, like itself pharmacologically an agonist) — p = 0.60, indistinguishable — not with C2. Since C3 shares C2's physical property (rigid, fills the pocket) but behaves like C1, H12-local confidence tracks **pharmacological class, not ligand rigidity**. That is exactly the question the 4th arm was added to answer, and it answers it.

### 3.2 Secondary: pocket-residue pLDDT (back-pocket residues 362, 379, 380, 396, 401, 480)

| arm | mean pocket pLDDT |
|---|---|
| apo | 92.32 |
| C1 | 92.73 |
| C2 | 93.85 |
| C3 | 93.10 |

| comparison | diff | p |
|---|---|---|
| C2 vs C1 | +1.13 | 0.0079 |
| C2 vs apo | +1.54 | 0.0079 |
| C3 vs C2 | −0.75 | 0.024 |
| C3 vs C1 | +0.37 | 0.15 (ns) |

Same direction as H12 (C2 highest among ligand-bound arms) but **5–10× smaller magnitude** (≤1.5 points vs. 5–10 points at H12). The ligand-dependence is concentrated at the state-defining residues, not a generic whole-pocket confidence shift.

### 3.3 Secondary: H12 coordinate RMSD (continuous signal, not collapsed into the binary H2 state call)

| arm | RMSD→3KYT (H12-in) | RMSD→4ZJW (H12-out) | margin (out−in) |
|---|---|---|---|
| apo | 0.54 Å | 7.71 Å | +7.17 |
| C1 | 0.60 Å | 7.70 Å | +7.10 |
| C2 | 1.52 Å | 5.97 Å | +4.44 |
| C3 | 1.05 Å | 6.90 Å | +5.85 |

All four arms still call decisively **H12-in** (nowhere near the 1 Å crossover). But C2 shows a real, consistent ~1–3 Å geometric loosening relative to apo/C1 — smaller margin, closer to 4ZJW — with C3 intermediate. This is sub-criterion signal the binary H2 gate would have discarded.

## 4. Interpretation (explicitly separated from the registered statistical result)

The primary statistical claim is as strong as this design can make it (p = 0.0079 across the key comparisons) and the confound question is cleanly resolved. The direction is the opposite of naive expectation going in: the model is **most** confident about H12 with no ligand at all, **least** confident with either agonist (sterol or synthetic), and intermediate with the inverse agonist — not "inverse agonist destabilizes H12 confidence, agonist stabilizes it."

One plausible (not registered, offered as interpretation only) explanation: both agonists make intimate contact with H12-adjacent pocket residues — residue 480 is reported in the kit's own pocket mapping as sitting *inside* H12 itself — and that direct contact, not pharmacological label, may be what perturbs local confidence. 4P1 (the inverse agonist) is bulkier and may engage the pocket more peripherally relative to H12, leaving the model free to default confidently to the agonist-locked pose regardless of the ligand's actual pharmacology, consistent with everything H1 and H2 already showed.

## 5. Project-level synthesis (H1 + H2 + H3)

Across three independent approaches — geometric state-calling (H1), a fine-grained MSA-depth/fold-stability sweep plus template and coactivator probes (H2 and its post-hoc follow-ups), and now model-confidence analysis (H3) — ligand pharmacology never moves this model's RORγ LBD prediction into the H12-out state. H3 adds a new, real, statistically airtight finding: the model's *confidence*, unlike its geometry, **is** ligand-sensitive — just not in the direction, or for the reason, that was hypothesized going in.

## 6. Known limitations (carried over, still apply)

1. n = 1 target (RORγ LBD), n = 4 ligand conditions. No cross-target claim is made.
2. pLDDT is the model's own confidence estimate, not a ground-truth disorder measurement. A pLDDT difference is evidence the model's *estimate of its own certainty* is ligand-dependent — not, by itself, evidence that the underlying biological state is ligand-dependent.
3. The explanatory account in §4 is interpretation, explicitly flagged as such, not a registered or tested claim.

## Appendix: data provenance

- Structures: `~/h3/out/{apo,c1,c2,c3}/` — 50 CIFs + scores.json per arm, apheris-openfold3 v0.15.1, full MSA (apheris-msa v0.9.0 fetch), weights `of3-ob-2025-06-30-174k.pt`.
- Per-residue pLDDT extraction and permutation tests: `~/h3/h3_perresidue.pkl`, computed from each structure's `atom_predicted_local_distance_difference_test` / `atom_tokens` arrays in its sibling `_scores.json`.
- Whole-structure pLDDT / RMSD-to-reference scoring: `~/h3/titration_results.csv`, via the project's `score_titration.py` against the kit's `rgkit` scorer (`state_recovery2`, `structio`) and `apheris_kit_rorgamma/refs_rorgamma/{3KYT,4ZJW}_clean.cif`.
- Pocket residue set (362, 379, 380, 396, 401, 480) and H12 region (479–486): `apheris_kit_rorgamma/targets_rorgamma.json`.
- BIO592/6F1 SMILES: RCSB Chemical Component Dictionary (`https://files.rcsb.org/ligands/view/6F1.cif`), cross-checked against PDB 5IZ0.

## 7. Post-report follow-up (same day, before final submission): contact control, window check, extended pool

Three checks run against the existing n=3-ligand dataset (apo/C1/C2/C3) plus a
newly launched 3-ligand extension, to see whether the §3/§4 finding survives
scrutiny and generalizes.

### 7.1 Contact control (`per_residue_discriminator.py`), re-run and confirmed

Re-ran the pre-registered arbitration test on the deposited C1/C2/C3 pool.
Result reproduces exactly what §4/§5 already concluded, now shown with the
actual contact footprints and numbers:

```
c1 vs c2:  H12 diff=-5.67 p=0.0079 d=-3.70   non-H12 diff=-3.37 p=0.0079 d=-5.75  -> CONTACT-PERTURBATION
c2 vs c3:  H12 diff=+6.20 p=0.0079 d=+3.78   non-H12 diff=+2.42 p=0.0079 d=+3.03  -> CONTACT-PERTURBATION
Overall: CONTACT-PERTURBATION model supported across all comparisons.
```

The non-H12 ligand-contact-residue effect size is equal to or larger than the
H12 effect size in both comparisons. **This control does not invalidate the
statistical finding (ligand identity still predicts H12-local pLDDT at
p=0.0079) — it confirms the §4 interpretation that the effect is a general
ligand-contact-region confidence suppression, not something unique to the
state-defining residues.** Report it that way; do not claim H12-specificity.

### 7.2 Window check: 479–486 vs. the deck's 501–507

The deck references an H12 window through residue 507; this project's
primary region (479–486) was auto-derived and is the only part of H12 both
reference crystals resolve (`INTEGRITY_NOTES.md`). 4ZJW does not resolve
502–507 at all, so an RMSD-to-4ZJW state call cannot be computed on that
window — but the **pLDDT** signal can be, since it only needs the predicted
structure, not a reference. Re-extracted per-residue pLDDT at 501–507 from
the same structures, same seed-level permutation test:

| window | apo | C1 (agonist) | C2 (inv. agonist) | C3 (agonist) | c1 vs c3 |
|---|---|---|---|---|---|
| 479–486 (primary) | 87.50 | 77.52 | 83.19 | 77.00 | p=0.60 (ns) |
| 501–507 (deck's window) | 79.59 | 56.31 | 67.17 | 58.40 | p=0.33 (ns) |

Same ordering (apo > C2 > C1 ≈ C3), same pairwise significances (C1 vs C2 and
C2 vs C3 both p=0.0079; C1 vs C3 not significant both times — the
rigidity-confound rejection replicates), substantially larger absolute
pLDDT gap at 501–507. **The separation holds on the challenge's own window —
the numbering objection does not change the finding.** The RMSD-based
geometric state call (H1/H2) still cannot be computed on 501–507 because
4ZJW lacks coordinates there; that limitation is unchanged and is about the
reference crystal, not this project's choice of window.

### 7.3 Extended ligand pool (n=3 → n=6): three more inverse agonists

The original H3 set had 2 agonists (HC2, BIO592) vs. 1 inverse agonist (4P1)
— a real design gap (n=3, 2-vs-1). Three more pharmacologically-annotated
inverse agonists already referenced in `refs_rorgamma/` were run overnight,
no conditioning, 5 seeds × 10 samples each (same protocol as C1–C3):

| arm | ligand | source PDB | pharmacology | note |
|---|---|---|---|---|
| C4 | 6EW / BIO399 | 5IXK | inverse agonist | co-crystallized H12-out |
| C5 | 99N | 5NTK | inverse agonist | co-crystallized H12-out (swung-out H12, extreme case) |
| C6 | 4P3 | 4ZJR | inverse agonist | **co-crystallized H12-in** — the kit's own "state-function decoupling decoy": pharmacologically inverse-agonist but the reference crystal keeps H12 packed. This is the compound that discriminates a pharmacology-tracking signal from a conformation-tracking one. |

SMILES pulled from RCSB CCD (`files.rcsb.org/ligands/view/{6EW,99N,4P3}.cif`).
Results for C4–C6, and classifier accuracy pooled across all 6 ligands
against known pharmacology labels, are in `h3_extended_pool_results.md`
(written once the jobs finished).
