# Report additions — RORγ H12 entry (Day 2)

Three drop-in sections for the final report/slides. Section 1 resolves a numbering
question with primary PDB data and MUST be reflected anywhere "H12" is named.

## 1. Numbering and reference-coverage note (integrity)

The challenge slide defines H12 as residues **501–507** (UniProt RORC numbering). The
analysis kit's auto-derived "H12 region" is **479–486** — the H11–H12 junction. Direct
inspection of the deposited PDB files (1 Oct 2026) establishes:

| Reference | State | Ligand | 479–486 | 501–507 |
|---|---|---|---|---|
| 3KYT | agonist ("in") | HC2 (25-hydroxycholesterol) — the C1 arm ligand | 8/8 modeled | 7/7 modeled |
| 4ZJW | inverse agonist ("out") | 4P1 — the C2 arm ligand | 8/8 modeled (displaced conformation) | **0/7 modeled (disordered)** |

3KYT also contains the coactivator peptide as chain C (residues 686–697, sequence
KHKILHRLLQDS — canonical LXXLL motif), which is the reference for the coactivator
co-fold comparison.

Consequences:

1. H1–H3 "H12-local" metrics were computed over the **junction segment (479–486)**.
   This segment is state-informative — it sits ~6–7.7 Å apart between the two
   references, and the H3 margins (+7.1 agonist vs +4.4 inverse agonist) quantify
   exactly that displacement — but it is not the segment the challenge scores.
2. For the true H12 (501–507), the "out" reference is **absence of electron density**,
   so RMSD-to-4ZJW is undefined by construction over that window. State scoring over
   501–507 must be confidence/geometry-based (pLDDT, RMSD-to-3KYT, helix geometry),
   not RMSD-to-4ZJW.
3. All H3 per-residue analyses have been re-extracted over 501–507 from existing
   outputs (zero GPU) and both windows are reported side by side. Figures label the
   window explicitly.

## 2. Corrected H3 headline (replaces "primary claim holds")

> H3's pre-registered primary contrast **passed**: local pLDDT on the state-defining
> region is higher with the inverse agonist than with either agonist
> (C2 − C1 = +5.67 points, complete seed-level separation, p = 0.0079 — the tightest
> achievable at 5-vs-5 — Cliff's δ = 1.0). The C3 arm (rigid synthetic agonist)
> **rejects the ligand-rigidity/pocket-filledness confound**: C3 is
> indistinguishable from the flexible sterol C1 (p = 0.60) and fully separated from
> the equally rigid inverse agonist C2 (p = 0.0079). The effect partitions by
> pharmacological class, not ligand physics.
>
> The pre-registered **ordering prediction failed**: apo is highest (87.5), not tied
> with the agonist arms. We therefore do **not** claim anti-correlation with
> crystallographic order. The supported claim is narrower: local confidence on the
> state-defining region is ligand-dependent and partitions by pharmacological class,
> with agonists *suppressing* confidence.
>
> We note a confidence–accuracy decoupling: C1 shows the best junction geometry
> (0.60 Å to the 3KYT reference) and the *lowest* confidence — on this target, local
> pLDDT does not rank local accuracy at the state-defining residues.

## 3. Classifier methods (for the 60% prediction-accuracy component)

> Pool state calls use a conservative two-feature voting rule
> (`h12_state_classifier.py`). Features per compound (median across seeds/samples):
> junction-window pLDDT (479–486), junction RMSD margin (RMSD→4ZJW − RMSD→3KYT), and —
> after re-extraction — H12-window pLDDT (501–507). A state call requires at least two
> agreeing votes; otherwise the compound is reported as "uncertain" rather than forced.
> Structures below the global pLDDT floor (70) are recorded as "unscorable" and listed,
> never silently dropped (project convention from H2). Thresholds were calibrated from
> the four H3 arms (known pharmacology: C1/C3 agonist, C2 inverse agonist) by midpoint
> grid search (`--calibrate`); defaults before calibration are the H3 midpoints
> (junction pLDDT 80.0; margin 5.8 Å). Where deposited-set labels are available, we
> report overall and per-class accuracy with the confusion matrix; where blind, we
> submit calls with calibrated confidence and say so.
>
> Limitations, stated once: thresholds are calibrated on three ligands of known
> pharmacology at a single target; the classifier reads a *confidence signal*, not a
> conformation — in 322 structures the model has never once adopted H12-out
> coordinates, so coordinate-level antagonist calls remain beyond the model on this
> target.
