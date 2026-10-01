# H2 — pre-registration

**Written 2026-10-01 13:11 UTC. No MSA-titration prediction had been run or
seen at the time of writing.** H1's results were complete and are summarised
below; nothing in this document may be edited once titration results exist.

---

## Status of H1 (retired, not revised)

**H1 (pre-registered, tested):** state-specific pocket conditioning steers
co-folding models between the RORγ H12-in and H12-out states, while naive /
shared-pocket conditioning does not.

**Result: refuted.** 0/25 predictions in the inverse-agonist direction recovered
H12-out, across five conditioning levers (pocket-specific, pocket-shared,
target template, wrong template, MSA-free) and their combinations. The baseline
recovered H12-out in 0/5 seeds, so the ≥60 % criterion sat in its strong regime
and nothing cleared it.

H1 is retired on its own evidence. H2 below is a **new hypothesis for new data**,
not a rewrite of H1. The H1 criterion and its negative result stand unchanged in
the final report.

### What H1 established, which motivates H2

| observation | value |
|---|---|
| H12 separation between the two crystal states (479–486) | 7.78 Å |
| All 22 structurally sound predictions, distance to H12-in | 1.55 ± 0.44 Å |
| …distance to H12-out | 5.91 ± 0.60 Å |
| Closest any prediction came to H12-out | 5.53 Å |
| Predicted ensemble spread vs inter-state distance | **17× tighter** |
| Full-MSA arms: mean pLDDT | 85–93 (fold good, H12 locked in) |
| Zero-MSA arms (A5/A7/B2): mean pLDDT | 43.5–43.7 (fold destroyed) |
| Zero-MSA arms: radius of gyration | 20.5–20.6 Å vs 18.3 Å crystal |

Template conditioning changed structures only at noise level and in no
consistent direction: no-template vs correct-state template = 0.625 Å;
no-template vs wrong template = 0.736 Å; correct vs wrong template = 0.807 Å.

**The two MSA endpoints bracket an untested interval.** Full MSA folds the
domain and locks H12-in. Zero MSA releases everything but destroys the fold.
Nothing between has been tried.

---

## H2

**OpenFold3's H12-in bias is governed by MSA depth rather than by pocket or
template conditioning. There exists an intermediate MSA depth at which the LBD
still folds and the H12-out state becomes accessible.**

Prior art: MSA subsampling to bias AlphaFold2 toward user-defined conformational
states (Sala, Hildebrand & Meiler) — cited on the challenge deck's own slide 12.

### Primary criterion (fixed before launch)

At some tested depth N, **≥ 2 of 5 seeds** produce a structure that satisfies
**both**:

1. **Fold-quality floor** — mean pLDDT ≥ 70, and
2. **State call** — region RMSD (479–486) to 4ZJW below that to 3KYT, with
   margin ≥ 1.0 Å.

Both conditions must hold in the same structure. A structure failing the floor
is recorded as `unscorable` and counts toward neither numerator nor denominator
of the fold-quality rate, but **does** count as a failure for the state-recovery
rate at that depth.

### Depths tested

2048, 512, 128, 32, 8 sequences, plus the two existing endpoints (full MSA from
H1, and zero MSA from H1's A5/A7/B2). Five seeds per depth: 42, 43, 44, 45, 46.
`--num-diffusion-samples 10`, so 50 structures per depth.

Depth is measured as **number of sequences retained in the .a3m after the query**,
taken from the top of the file as returned by `apheris-msa fetch-msa` — i.e. no
re-ranking, no diversity selection. Stated here because the subsampling policy
changes what the experiment means and must not be chosen after seeing results.

### Secondary readouts (reported beside the primary, never merged into it)

- mean pLDDT per depth — the fold-stability curve
- ligand centroid distance and pocket-contact recall — pose, reported separately
  from state (pose ≠ state)
- radius of gyration — independent fold-integrity check

### Declared negative

If no depth satisfies the primary criterion, the conclusion is:

> On RORγ LBD, fold stability and H12-state accessibility cannot be separated by
> MSA depth. The state is inaccessible at every depth that preserves the fold.

That is a reportable result and will be reported as the headline if it occurs.
It is a stronger claim than H1's negative, because it closes the mechanism that
the literature suggests should work.

### Threats to validity, named now

1. **Depth is confounded with fold quality.** Shallower MSAs both release
   constraints and degrade folding. The pLDDT floor separates these, but a
   result at the floor's edge is ambiguous and will be reported as such.
2. **Subsampling policy.** Taking the top N sequences retains the most similar
   homologues. A diversity-based policy might behave differently; not tested,
   and that limitation will be stated.
3. **One target, one ligand.** Nothing here generalises beyond RORγ LBD with 4P1
   without further work.
4. **Single architecture.** All H1 predictions and all H2 predictions are
   OpenFold3. The only Boltz-2 data point is the H1 smoke run (A0, H12-in,
   ligand 0.71 Å). No cross-model claim will be made beyond that single
   comparison, which is n=1 per model.

---

## Sign-off

Criterion fixed at 2026-10-01 13:11 UTC, before any titration job was submitted.
Any later change to this file must be recorded as a dated amendment below, with
its reason, and the original criterion left legible.

### Amendments

(none)
