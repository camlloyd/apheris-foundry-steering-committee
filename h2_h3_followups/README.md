# RORγ H12 steering — reproducibility package

> **⚠ Read `results/OFFICIAL_RESCORE_FINDINGS.md` first.** Everything below
> was scored against this project's own chosen reference structures
> (3KYT/4ZJW, window 479–486), not the organizers' official pair
> (5VB7/6T4I, window 484–507). Re-scored against the official kit, the
> headline "0/941 ever reached H12-out" does not hold — 209/1,191
> structures (17.2% of well-folded ones) call decisive antagonist. H2's
> negative is the one result that replicates unchanged. Numbers below are
> kept as-is (scored against our own reference, honestly labeled); the
> corrected numbers are reported side by side in that file, not substituted
> in place of these.

Single folder: what ran, in what order, with what settings, on what compute,
and how to re-run any of it. Everything here was produced on this VM with
pre-pulled Apheris Foundry module images (`apheris-openfold3-v0.15.1`,
`apheris-msa-v0.9.0`, `apheris-data-v0.26.0`) run directly via `docker run ...
module-run ...` — no custom training, no synthetic-structure fine-tuning (the
documented dead end was never attempted).

## Story in one paragraph

H1 (pocket/template steering, nine levers, 41 predictions scored, 0/25 in the
inverse-agonist direction recovered H12-out) and H2 (MSA-depth titration,
fine grid at 8/32/48/64/96/128/512 sequences, 200 core + 350 post-hoc
structures — the fold-stability cliff bisected at 48/64/96, a template arm
re-run against a genuinely alternative-pose reference structure, and a
coactivator co-fold) together account for **591 scored structures, 0 of
which ever reached H12-out** by the pre-registered geometric criterion (mean
pLDDT ≥ 70 AND region-RMSD margin ≥ 1 Å). (Earlier drafts of this project's
own documents used a "322" figure for this total — that number did not trace
back to a real count and has been corrected here; 591 is independently
verified against the structures on disk, see the counts in each round's
`results/` CSV.)
H3 then asked a narrower question — is the model's *confidence*, not its
geometry, ligand-dependent? Yes, with the tightest possible statistical
confidence (p = 0.0079, complete 5-vs-5 seed separation) — but a follow-up
discriminating analysis found the effect is a **general ligand-contact-region
confidence suppression**, not something specific to the state-defining
residues. `results/H3_RESULTS_REPORT.md` has the full writeup; `h12_state_classifier.py`
turns that signal into a usable (if imperfect — 88% self-test accuracy,
64% specifically on the inverse-agonist class) per-compound classifier.

H4 then asked whether coactivator-peptide *placement* (not state, not
confidence) tracks ligand identity. It does, with the same statistical
confidence as H3 (p = 0.0079, complete 5-vs-5 seed separation) — but
backwards from the pharmacologically expected direction: the inverse
agonist (D2, 4P1) query reproduces the peptide's real crystallographic
pose almost exactly (4.3 Å), while the agonist (D1, HC2) query folds it
into a different location entirely (37.1 Å), with more peptide-receptor
contacts in D2 than D1. `results/H4_RESULTS_REPORT.md` has the full
writeup, including why this isn't template leakage and what's registered
vs. interpreted.

## Folder contents

```
preregistrations/   H2, H3, and H4 pre-registrations, each timestamped
                     before any job in that round ran. None were edited
                     afterward; H3's amendment log is empty.
scripts/             the runnable pipeline:
  run_msa_titration.sh       MSA fetch + depth-ladder + predict-openfold3,
                             real fixed docker invocations (see inline
                             comments for what was wrong in the first draft
                             and why — missing --weights, wrong --seeds
                             syntax, wrong query JSON shape, etc.)
  run_metrics.sh              apheris-data compute-metrics, probe + run,
                             manifest-based (the only input form that scores
                             predictions independently rather than as one
                             ensemble)
  score_titration.py          pLDDT floor + H12 RMSD state call, against
                             the pre-registered criterion
  local_metrics.py            independent local reimplementation of
                             GDT-HA/lDDT/LDDT-PLI/ligand RMSD, cross-checked
                             against Foundry's own compute-metrics numbers
  h12_state_classifier.py     the H3-derived confidence classifier --
                             READ ITS DOCSTRING, the calibration caveat is
                             load-bearing
  per_residue_discriminator.py the contact-perturbation-vs-H12-specific
                             arbitration test, reusable on new ligand sets
  collect_compute_stats.sh    hardware, pinned image digests, wall time,
                             structure counts
  rgkit/                      the shared scoring library (structure I/O,
                             state-recovery RMSD) all of the above import
results/             every CSV/pickle/report the scripts above produced
logs/                driver logs for each round (H2, H2 post-hoc, H3) --
                     real docker stdout, including the errors that got fixed
                     along the way, left in rather than cleaned up
INTEGRITY_NOTES.md    numbering verification + what the H12-out reference
                     structure actually resolves (read before quoting any
                     H12 RMSD number against a specific residue range)
raw_runs/            the actual working directories behind every round
                     (configs, logs, query inputs, MSA files, result
                     metadata) -- see raw_runs/README.md for what's
                     tracked here vs. regenerated from a manifest
```

## Re-running any of it

```bash
export OF3=quay.io/apheris/foundry-hackathon:apheris-openfold3-v0.15.1
export MSA=quay.io/apheris/foundry-hackathon:apheris-msa-v0.9.0
export DATA=quay.io/apheris/foundry-hackathon:apheris-data-v0.26.0
KIT=/path/to/apheris_kit_rorgamma   # refs_rorgamma/, targets_rorgamma.json

cd scripts
bash run_msa_titration.sh                                   # H2 grid
python3 score_titration.py --work ~/h2 --kit "$KIT" --rgkit .
bash run_metrics.sh probe && bash run_metrics.sh run ~/h2/out
python3 local_metrics.py --refs "$KIT/refs_rorgamma" --runs ~/h2/out \
    --ref-state H12-out --smiles '<ligand SMILES>' --out metrics.csv
bash collect_compute_stats.sh > compute_usage.md

# on any new pool of predictions (per-compound), once structures exist:
python3 h12_state_classifier.py --work <pool_predictions_dir> --kit "$KIT" \
    --rgkit . --out state_calls.csv --labels <optional_ground_truth.csv>
python3 per_residue_discriminator.py --work <dir> \
    --arms <arm1>:<label1> <arm2>:<label2> ...
```

Every stage is idempotent (skips completed depths/arms), so an interrupted
run resumes.

## Known limitations (carried through from every report in `results/`)

- n = 1 target (RORγ LBD), small ligand set (HC2, 4P1, BIO592). No
  cross-target claim is made anywhere in this project.
- The classifier's threshold is fit on that same n = 3 ligands — report its
  accuracy with that caveat attached every time, not as a general number.
- pLDDT differences are evidence about the model's confidence estimate, not
  direct evidence about the underlying biological state.
- See `INTEGRITY_NOTES.md` for the H12-window numbering clarification before
  quoting any RMSD figure against a specific residue range on a slide.
