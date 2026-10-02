# H3 extended pool — three more ligands, contact control, window check

> **⚠ See `OFFICIAL_RESCORE_FINDINGS.md` in this same folder** for the
> official-reference re-score, which upgrades §3's C6 finding from
> "the pLDDT signal direction suggests pharmacology over conformation" to a
> direct coordinate-level confirmation (C6 calls decisive agonist 19/50,
> 0/50 decisive antagonist against the real 6T4I reference) — and supersedes
> this file's §4 classifier-accuracy numbers with a real-label ChEMBL pool
> test (75.8–91.6%, not a self-referential n=3-calibrated classifier).

Run same day as `H3_RESULTS_REPORT.md`, after that report was written, as a
set of controls meant to either strengthen or break the H3 finding before
submission. Full protocol match to C1–C3: full LBD construct, full MSA
(~2020 seqs), no template/pocket conditioning, seeds 42–46, 10 diffusion
samples/seed, 50 structures/arm. SMILES for the three new ligands pulled
directly from the RCSB Chemical Component Dictionary
(`files.rcsb.org/ligands/view/{6EW,99N,4P3}.cif`).

| arm | ligand | source PDB | pharmacology (literature) | role |
|---|---|---|---|---|
| C4 | 6EW / BIO399 | 5IXK | inverse agonist | co-crystallized H12-out |
| C5 | 99N | 5NTK | inverse agonist | co-crystallized H12-out, extreme swung-out H12 |
| C6 | 4P3 | 4ZJR | inverse agonist | co-crystallized **H12-in** — the kit's own pre-registered "state-function decoupling decoy" (`targets_rorgamma.json`): pharmacologically inverse-agonist, but the reference crystal keeps H12 packed in the agonist position. |

## 1. Geometric state call — unchanged negative, now over more structures

Every one of the 150 new structures (C4+C5+C6) still calls decisively
**H12-in** by the pre-registered RMSD criterion; none approaches the ≥1 Å
margin. One structure in C5 (99N) came within 0.01 Å of the threshold — the
closest any structure in this whole project has come — but still did not
cross.

| arm | mean margin (out−in) | min margin | max margin | structures crossing ≥1 Å |
|---|---|---|---|---|
| apo | −7.17 | −7.41 | −6.98 | 0/50 |
| C1 | −7.10 | −7.43 | −3.24 | 0/50 |
| C2 | −4.44 | −7.34 | −2.52 | 0/50 |
| C3 | −5.85 | −7.30 | −2.36 | 0/50 |
| C4 | −6.15 | −7.32 | −2.97 | 0/50 |
| C5 | −4.12 | −6.07 | **−0.01** | 0/50 |
| C6 | −6.02 | −7.39 | −3.03 | 0/50 |

Combined with H1+H2 (591 structures, 0 crossings, see `README.md`), the
project-wide geometric negative is now **0/941 structures across four
independent approaches and 7 ligand conditions** ever reached H12-out by the
pre-registered criterion.

## 2. H12-local pLDDT — the n=3 "class" pattern does NOT generalize cleanly to n=6

This is the headline finding of this follow-up, and it cuts against the
clean story in `H3_RESULTS_REPORT.md` §3–4. Re-running
`per_residue_discriminator.py` over all six ligand arms (full output:
`h3_contact_control_extended.txt`):

| arm | pharmacology | H12-local pLDDT (mean of 5 seed means) |
|---|---|---|
| apo | — | 87.50 |
| C1 (HC2) | agonist | 77.52 |
| C3 (BIO592) | agonist | 77.00 |
| C2 (4P1) | inverse agonist | 83.19 |
| **C4 (6EW/BIO399)** | inverse agonist | **77.56** |
| **C5 (99N)** | inverse agonist | **63.34** |
| **C6 (4P3)** | inverse agonist | **82.22** |

C2 and C6 pattern together (83.2, 82.2 — both clearly above the agonists).
**C4 is statistically indistinguishable from both agonists** (vs C1:
diff=−0.04, p=0.97; vs C3: diff=−0.56, p=0.59) — it does not separate at
all. **C5 is the most extreme value in the entire dataset, lower than
either agonist** (63.3 vs ~77), the opposite direction from what a
"inverse-agonist ⇒ higher H12 pLDDT" class rule predicts.

Pairwise permutation verdicts (contact-perturbation arbitration test,
non-H12 ligand-contact residues vs. H12):

```
c1 vs c2: CONTACT-PERTURBATION   c1 vs c4: AMBIGUOUS (no separation at all)
c1 vs c5: H12-SPECIFIC            c1 vs c6: CONTACT-PERTURBATION
c2 vs c3: CONTACT-PERTURBATION   c3 vs c4: AMBIGUOUS
c3 vs c5: AMBIGUOUS               c3 vs c6: AMBIGUOUS
```

**Honest reading:** the clean, statistically airtight 2-agonist-vs-1-
inverse-agonist separation reported in §3.1 of `H3_RESULTS_REPORT.md` does
not hold up as a general "pharmacological class" rule once more ligands of
the same nominal class are added. Of the three new inverse agonists, one
(C6/4P3) replicates the original pattern, one (C4/6EW) shows no detectable
effect at all, and one (C5/99N) goes in the opposite direction more
strongly than any agonist. **This is exactly the n=3 → n=6 stress test the
pre-submission plan called for, and it does not survive it as a class
claim.** The correct scope for the §3/§4 finding is now: *"H12-local pLDDT
is ligand-identity-dependent and the effect is concentrated at
ligand-contact residues generally, not specifically pharmacology-encoded —
demonstrated on 4P1/HC2/BIO592 (n=3) and not replicated as a clean class
effect on three additional inverse agonists (n=3 more)."* Do not state a
general "pharmacology, not rigidity" rule on the slide; state the n=3
result and flag that the n=6 extension did not confirm it generalizes.

## 3. The 4P3/C6 discriminator result (the one genuinely clean finding here)

C6 (4P3) is pharmacologically an inverse agonist but is co-crystallized with
H12 still packed in the agonist position (4ZJR, the kit's own
"state-function decoupling decoy"). The question posed in the pre-
submission plan: does the model's confidence signal track the ligand's
*pharmacology* (read inverse-agonist-like) or the *reference crystal's
conformation* (read agonist-like)?

- H12-local pLDDT for C6: 82.22 — closer to C2 (83.19, the calibration
  inverse agonist) than to either agonist (77.0–77.5).
- The per-compound classifier (`h12_state_classifier.py`, threshold fit on
  n=3) calls C6 **"inverse-agonist-like" in 48/50 structures (96%)**.

**Result: the signal tracks pharmacology, not the reference crystal's
resolved conformation**, for this one compound. This is the cleanest single
result in this follow-up, exactly because it was designed in advance to be
decisive either way — but it is one compound, and section 2 above shows the
same signal fails to separate at all for a different inverse agonist (C4).
Report both facts together; do not lead with C6 alone.

## 4. Pooled classifier accuracy against known pharmacology (n=6 ligands, 300 ligand-bound structures)

`h12_state_classifier.py --work ~/h3 --labels h3_labels_extended.csv`,
full output in `h3/state_calls_pool_c1-c6.csv`. Apo (50 structures) excluded
per the classifier's documented scope (it does not attempt to detect apo).

```
agonist-like:         100/100 = 100.0%
inverse-agonist-like:  92/200 =  46.0%
OVERALL:              192/300 =  64.0%
```

Per-arm breakdown (not in the script's own output, computed from the same
CSV):

| arm | pharmacology | called 'agonist-like' | called 'inverse-agonist-like' |
|---|---|---|---|
| C1 | agonist | 49/50 | 1/50 |
| C3 | agonist | 44/50 | 6/50 |
| C2 (calibration ligand) | inverse agonist | 18/50 | 32/50 (64%) |
| C4 | inverse agonist | 38/50 | 12/50 (24%) |
| C5 | inverse agonist | 50/50 | 0/50 (0%) |
| C6 | inverse agonist | 2/50 | 48/50 (96%) |

**64% overall accuracy is real but almost entirely a reflection of the
agonist class being easy (100%) and the inverse-agonist class being a coin
flip at best (46%), with huge, inconsistent per-ligand variance (0% to
96%)** — not a stable property of "inverse agonists" as a class. The
classifier generalizes acceptably to new agonists (C4's and C5's own
agonist-like calls happen to be numerically "wrong" by label since both are
inverse agonists, but note C1/C3 — the actual agonists — are called
correctly essentially every time). It does **not** generalize reliably to
new inverse agonists: performance ranges from worse-than-chance (C5, 0%) to
excellent (C6, 96%) depending on the specific compound. The documented
calibration caveat in `h12_state_classifier.py` (threshold fit on n=3, not
validated outside that set) is fully borne out by this test — report the
64% number only with that caveat attached, and report the per-arm spread
alongside it, not instead of it.

## 5. What this changes about the project's claims

- **Does not invalidate:** the n=3, pre-registered, p=0.0079 H3 result on
  HC2/4P1/BIO592, the contact-perturbation interpretation of it (§1 above
  reproduces it exactly), or the project-wide geometric negative (now
  0/941, strengthened by 150 more structures).
- **Does invalidate:** any framing of the H3 finding as a general
  "pharmacological class" signal, or of the classifier as reliable on
  unseen inverse agonists. Both must be scoped to n=3/n=6 explicitly, with
  the C4/C5/C6 spread shown, not summarized away.
- **Adds one clean result:** C6/4P3 — the model's H12 confidence tracks
  this compound's pharmacological label, not its reference crystal's
  resolved (agonist-locked) conformation. Reported as a single-compound
  finding, not generalized.

## Appendix: data provenance

- Structures: `~/h3/out/{c4,c5,c6}/` — 50 CIFs + scores.json per arm, same
  pipeline/weights/MSA source as C1–C3 (`run_h3.sh` pattern, extended to
  three more query JSONs under `~/h3/in_c{4,5,6}/query.json`).
- Contact control: `h2_h3_followups/scripts/per_residue_discriminator.py`,
  raw output in `h3_contact_control.txt` (n=3) and
  `h3_contact_control_extended.txt` (n=6).
- Window check: `h2_h3_followups/scripts/h3_window_check.py`, raw output in
  `h3_window_check.txt`.
- Pooled classifier + per-arm accuracy: `h12_state_classifier.py`, labels
  file `~/h3/h3_labels_extended.csv`, scored CSV
  `~/h3/state_calls_pool_c1-c6.csv`.
- SMILES: RCSB CCD, `files.rcsb.org/ligands/view/{6EW,99N,4P3}.cif`.
