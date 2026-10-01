# H3 — Results Report

**Hypothesis registered:** `H3_PREREGISTRATION.md`, timestamped 2026-10-01 14:56 UTC, before any H3 job ran.
**This report written:** 2026-10-01, after all 4 arms completed and scored.

---

## 1. Background

H1 (pocket/template steering of RORγ LBD toward the H12-out inverse-agonist
state, 41 predictions, 0/25 in the inverse-agonist direction) and H2
(MSA-depth titration plus post-hoc follow-ups — fine-grained fold-stability
sweep, an alternative-pose template arm, a coactivator co-fold; 550
structures) were both retired as negatives: across **591 structures total**
and nine-plus conditioning levers, **0** ever reached H12-out by the
pre-registered geometric criterion (mean pLDDT ≥ 70 AND region-RMSD margin
≥ 1 Å to the correct reference). Whenever the model folds the domain at all,
it defaults to the agonist-locked H12-in pose.

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
