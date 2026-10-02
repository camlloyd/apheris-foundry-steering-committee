# RORγ H12 steering

Can you tell a nuclear receptor to switch off, just by showing a co-folding
model the right ligand?

RORγt is the master transcription factor for Th17 cell differentiation — a
validated autoimmune/inflammatory drug target (psoriasis, psoriatic
arthritis) where several inverse-agonist drug programs have stalled on
selectivity. Switching it off is not a binding-site occupancy question, it's
a *conformational* one: Helix 12 (H12), the receptor's C-terminal activation
helix, either packs down to close the coactivator groove (agonist, "on") or
gets displaced so a corepressor binds instead (inverse agonist, "off").
Whether a ligand flips that switch is exactly what this hackathon asked a
structure-prediction model to get right, across two distinct mechanisms
(orthosteric steric clash and allosteric long-range destabilization) that
look nothing alike in the pocket.

This repo is the full record of trying to answer that with Foundry
(`apheris-openfold3`), four pre-registered hypotheses deep, including the
place where we found our own mistake the night before judging and re-ran
everything rather than keep the more flattering number.

## What we found, in order

| | question | n | result |
|---|---|---|---|
| **H1** | does pocket/template conditioning steer the model from H12-in to H12-out? | 51 | **no** — 0/35 in the inverse-agonist direction, including the pre-registered "money arm"; 16/16 confirm the agonist default works |
| **H2** | does MSA depth unlock the state once the fold is stable? | 550 | **no** — 0/550 at any depth that folds the domain at all |
| **H3** | is the model's *confidence*, not its geometry, ligand-dependent? | 200 | **yes** — p = 0.0079, complete 5-vs-5 seed separation — but in the opposite direction from the naive story (most confident apo, least confident with *either* ligand) |
| **H4** | does coactivator-peptide placement depend on ligand identity? | 100 | **yes, and backwards** — the *inverse agonist* query reproduces the real crystal peptide pose almost exactly (4.3 Å); the *agonist* query doesn't (37.1 Å) |

Then, scored against the organizers' own reference pair and ChEMBL labels
instead of the reference structures we'd chosen ourselves:

- **The headline "0 ever reaches the antagonist state" does not survive the
  official reference.** 209/1,191 structures (17.2% of well-folded ones)
  call decisive antagonist, including 60% of H1's own unconditioned
  baseline. We found this ourselves, the night before judging, and re-ran
  the full re-score rather than present the more dramatic number.
- **H2's negative is the one result that replicates exactly** either way —
  MSA depth genuinely never separates fold stability from state access.
- Against 258 real ChEMBL agonist/antagonist labels never seen by this
  project before scoring: **69.0% pooled accuracy** (77.3% agonist / 60.8%
  antagonist) — matching the organizers' own stated asymmetry
  ("antagonist is the hard half") almost exactly, measured independently.
- We then spent a morning trying to break that number before trusting it:
  pre-registered a decision-boundary recalibration (the boundary was
  already optimal — the asymmetry is real class difficulty, not a
  miscalibration), swept accuracy against coverage by class (abstention
  helps agonist calls, not antagonist ones — a real, reported limitation,
  not smoothed over), and stratified the errors by scaffold, confidence,
  and margin.

Full numbers, every caveat, and the statistics behind each line are in
`JUDGES_ROLLUP.md` — read that first if you're judging this. The slide
version is `PITCH_SLIDES.md`.

## How this repo is organized

```
JUDGES_ROLLUP.md          Start here. The synthesized story end to end,
                           every result, every correction, what to say and
                           what not to overclaim.
PITCH_SLIDES.md           The deck, as marp-flavored markdown.
HANDOFF.md                The Day-1 plan this project worked from: target,
                           hypothesis, experimental grid, task split.

apheris_kit_rorgamma/     H1 — the pocket/template/MSA-free steering grid
                          (15 arms, both ligand directions) and the scoring
                          library it's built on. Its RUNBOOK_rorgamma.md
                          has every command; raw predicted structures
                          aren't tracked (regenerable, see below) but every
                          reference crystal and per-structure confidence
                          score is.

h2_h3_followups/          H2 (MSA-depth titration), H3 (confidence vs.
                          ligand identity), H4 (coactivator placement), the
                          official re-score against the organizers' own
                          reference and ChEMBL pool, and the morning
                          stress-tests of that pool's accuracy numbers.
                          README.md there explains the full re-run recipe;
                          every round has its own timestamped
                          pre-registration in preregistrations/ and its own
                          results/*_REPORT.md.

foundry-ops/              The three judging-rubric metrics as tested,
                          reusable code: GDT-HA (via apheris-data),
                          lDDT-PLI, and ligand BiSyRMSD (RDKit,
                          symmetry-aware). uv-managed, Docker-reproducible,
                          CI-tested.
```

## What's tracked, what isn't

Every script, config, pre-registration, result, report, reference crystal
structure, and per-structure confidence score (`*_scores.json`) in this
repo is real and re-derivable from what's here. The one thing deliberately
*not* tracked is the raw predicted structures themselves
(`apheris_kit_rorgamma/runs/**/*_model.cif`) — they regenerate from the
commands already in `apheris_kit_rorgamma/RUNBOOK_rorgamma.md`, so
committing them a second time would just be carrying weight.
`apheris_kit_rorgamma/RAW_PREDICTIONS_MANIFEST.csv` records every one that
existed at commit time — path, size, sha256 — so a bad regeneration is
detectable without carrying the bytes.

## Running any of it

Every round is independently reproducible from its own directory — see
`apheris_kit_rorgamma/RUNBOOK_rorgamma.md` for H1, and
`h2_h3_followups/README.md` for H2 through H4 and the pool re-score. Both
assume pre-pulled Foundry module images
(`apheris-openfold3`, `apheris-msa`, `apheris-data`) run directly via
`docker run ... module-run ...` — no custom training, no synthetic-data
fine-tuning (a documented dead end per the organizers' own notes, never
attempted here).

## What this project is honest about

- n = 1 target (RORγ LBD), a handful of ligands directly, 258 compounds in
  the ChEMBL accuracy test — no cross-target claim is made anywhere.
- pLDDT and RMSD are the model's own outputs, not ground truth about real
  receptor biology.
- The orthosteric/allosteric mechanism split the brief calls "half the
  problem" could not be tested against real labels — the data to do it
  with doesn't exist in what we were given, and we say so rather than
  quietly work around it.
- Every number above that turned out to be wrong, overstated, or weaker
  than it first looked is reported that way in `JUDGES_ROLLUP.md`, not
  corrected quietly in place.
