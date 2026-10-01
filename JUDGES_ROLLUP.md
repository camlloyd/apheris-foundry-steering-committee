# RORγ H12 Steering — Rollup for Judges

Synthesizes H1 (`apheris_kit_rorgamma/`) + H2/H3 (`h2_h3_followups/`) into the
single artifact `HANDOFF.md` §8 asks for. Supersedes the pitch line in
`HANDOFF.md` §9, which was written before any result and is now stale —
`HANDOFF.md` itself is left unedited as the pre-registered record; this file
is the honest update.

## Updated pitch line

~~"State-specific pocket conditioning steers co-folding models between RORγ
agonist and inverse-agonist H12 states."~~ (pre-registered, not supported)

**"Across three independent approaches and 801 scored structures (601 testing
state recovery directly, plus 200 testing model confidence), ligand
pharmacology and pocket/template conditioning never move this model's RORγ
LBD prediction into the H12-out state — a rigorous, pre-registered negative
result — but the model's own *confidence* at the state-defining residues
(H12, 479–486) is ligand-dependent with the tightest possible statistical
significance (p = 0.0079, complete 5-vs-5 seed separation), in a direction
that rules out a naive rigidity confound."**

## What actually happened, by round

| round | question | design | n | result |
|---|---|---|---|---|
| H1 | does state-specific pocket/template conditioning steer Boltz-2/OpenFold3 between H12-in/H12-out? | 15 arms (A0–A11, B0–B3), real Foundry predictions | 35 (invago dir.) + 16 (25-HC dir.) = 51 scored | **0/35 reach H12-out in the invago direction** (incl. A4, the pre-registered "money arm"), **16/16 reach H12-in in the 25-HC direction** (apo/agonist default is confirmed, matches H1 hypothesis) |
| H2 | does MSA depth separate fold-stability from H12-state accessibility? | depths 8/32/48/64/96/128/512, 5 seeds | 200 core + 350 post-hoc = 550 | **0/550 reach H12-out at any depth that also folds the domain** (pLDDT≥70) — the state is inaccessible everywhere the fold is trustworthy |
| H3 | is the model's *confidence*, not geometry, ligand-dependent? | apo/C1(agonist)/C2(inverse agonist)/C3(agonist, added to break a confound), 5 seeds×10 samples | 200 | **Yes** — H12-local pLDDT differs by ligand (p=0.0079); ordering is apo > C2 > C1≈C3 — confidence tracks pharmacological class, not rigidity, but direction is opposite the naive hypothesis |

**Total structures scored: 51 (H1) + 550 (H2) + 200 (H3) = 801.** Zero, across
H1+H2 (601), ever reached H12-out by the pre-registered criterion (mean
pLDDT ≥ 70 AND region-RMSD margin ≥ 1 Å).

> **Open discrepancy, not resolved here:** `h2_h3_followups/README.md` states
> H1 as "41 predictions scored, 0/25 in the inverse-agonist direction."
> Re-scoring the actual prediction files on disk today
> (`apheris_kit_rorgamma/runs/rorgamma_invago/`,
> `.../rorgamma_25hc/`, via `score_run.py score`) gives 35 invago + 16 25-HC
> = 51 total, and 0/35 (not 0/25) in the invago direction — every subset I
> tried (all arms, A-only, core-7-only: 19/10) failed to reproduce 41/25
> exactly. The *qualitative* result is identical either way (0 recoveries),
> so nothing above depends on which count is right, but whoever owns the
> pitch deck should reconcile this before quoting an exact H1 n on a slide.
> Re-run: `cd apheris_kit_rorgamma && python3 score_run.py score --config
> targets_rorgamma.json --target rorgamma_invago --runs
> runs/rorgamma_invago --out results_invago.csv` (ditto `rorgamma_25hc`).

## Deliverables status (HANDOFF.md §8)

| deliverable | status |
|---|---|
| `figure_day1.png/.svg` + `summary_table.md` | **done this session**, split by direction (the merge HANDOFF.md implies doesn't make sense here — H1 ran two directions with different target states): `apheris_kit_rorgamma/figure_day1_invago.{png,svg}` + `summary_table_invago.md`, `apheris_kit_rorgamma/figure_day1_25hc.{png,svg}` + `summary_table_25hc.md` |
| `results.csv` | **done this session**: `apheris_kit_rorgamma/results_invago.csv`, `results_25hc.csv` (real re-score of the on-disk H1 predictions); H2/H3 equivalents already existed at `h2_h3_followups/results/{h2,h3}_titration_results.csv` and `h2_posthoc_results.csv` |
| `targets_rorgamma.json` | done (pre-existing) |
| `rorgamma_addon.zip` | done (pre-existing, root of repo + `h2_h3_followups/H1_handoff_package_rorgamma_addon.zip` backup) |
| evening runbook | done (pre-existing): `RUNBOOK_rorgamma.md`, `h2_h3_followups/README.md` |
| pitch outline | **this file** — the first pass at a standalone one |

## Judging-rubric metrics (state-call correctness + mean of GDT-HA, lDDT-PLI, Ligand BiSyRMSD)

All three implemented and pytest-covered in `foundry-ops/`:
`rorgamma_metrics.py` (lDDT-PLI), `rorgamma_ligand_rmsd.py` (BiSyRMSD via
RDKit), GDT-HA covered directly by apheris-data's `compute-metrics` task.
State-call correctness is `apheris_kit_rorgamma/state_recovery2.py`
(crystal-validated 5/5, see `RUNBOOK_rorgamma.md`). Nothing outstanding here.

## What to actually say to judges

1. A real, controlled, pre-registered negative result (H1+H2, 601
   structures, nine conditioning levers) — the null is reported as a
   finding per the project's own pre-registration, not hidden.
2. A genuine positive finding growing out of that negative (H3): the model's
   confidence, not its geometry, is ligand-sensitive, at the tightest
   possible statistical confidence this design permits — and a follow-up
   arm (C3) was added specifically to rule out the obvious confound, which
   it did.
3. A working, tested, judging-rubric-compliant metrics pipeline.
4. A tangible secondary artifact: `h12_state_classifier.py`, a per-compound
   confidence-based classifier derived from the H3 signal (88% self-test
   accuracy, 64% on the inverse-agonist class specifically — report with
   that caveat, n=3 ligands).

## Known limitations (unchanged, carried through every report)

n=1 target (RORγ LBD); small ligand set; pLDDT is the model's own
confidence estimate, not ground-truth disorder; see
`h2_h3_followups/INTEGRITY_NOTES.md` before quoting any H12 RMSD number
against a specific residue range — the reference pair resolves 479–486 but
not the deck's broader 501–507 window.
