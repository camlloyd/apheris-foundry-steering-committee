# Day-of command cheatsheet (extension modules)

Everything assumes you're in the kit directory and the pharma data drop is at
`$DROP` (a directory of structures, or structures + metadata CSV).

## 10:30 — after the opening talks

```bash
# 1. What does the Foundry Hub actually expose? (2 min, decides the grid)
apheris-foundry workflows list
python foundry_emit.py probe --workflow predict
#    -> prints which levers (template / pocket / msa / seeds) the schema
#       supports. Unsupported levers become the product-feedback slide.
```

## 11:00 — triage the data drop (45 min)

```bash
python triage.py $DROP --metadata $DROP/meta.csv --cutoff-date 2023-06-01 \
    --out triage_report.csv --write-targets pharma_targets.json --top 3
#    -> ranked target table with printed reasons. CONFIRM the top pick by
#       eye (open both reference structures). You want: two states, both
#       liganded, post-cutoff, not famous.
#    -> fill in the "smiles" field in pharma_targets.json before continuing.
```

## 11:45 — validate, generate, submit as ONE batch

```bash
python score_run.py validate --config pharma_targets.json --target <pick>
python make_inputs.py --config pharma_targets.json --target <pick> --out runs/
python foundry_emit.py emit --runs runs/<pick> --out foundry_jobs/
bash foundry_jobs/submit_all.sh        # all arms, one batch, then walk away
```

## 12:30 onwards — score as results land

```bash
python foundry_emit.py collect --jobs <download_dir> --runs runs/<pick>
python score_run.py score --config pharma_targets.json --target <pick> \
    --runs runs/<pick> --out results.csv
python figure.py results.csv --target-state <target_state> \
    --out-prefix figure_day1 --md summary_table.md
```

## Fallbacks

- **`validate` fails on a reference with "unassigned" or wrong state** -> the
  motif autodetect picked wrong anchors. Two known causes, both fixed in the
  kit but worth recognizing: (1) DFG-out references have a broken K-E salt
  bridge (alphaC-out), which defeats distance-only anchor search; (2) unmodeled
  motif residues (1KV2's DFG Gly170) break a contiguous sequence search.
  Fix by hand: add a `motif_overrides` block to the target config (see
  targets.json for worked examples) with the anchors counted from the
  structure in PyMOL, then re-run validate.
- **Hub has no pocket/template params** -> run the raw Boltz-2 YAMLs /
  OpenFold3 JSONs from `runs/<pick>/` on the Lyceum GPUs directly; the grid
  is identical, only the runner changes.
- **Pharma data is sequences+SMILES only** -> references must come from the
  PDB; use the kit's ABL1/p38a configs as templates and lean on the
  wrong-template control for the memorisation argument.
- **Non-kinase targets** -> scoring still works: state_recovery2 derives the
  state-defining region from the two references automatically. The DFG
  columns in the CSV just stay empty.
- **Nothing two-state in the drop** -> pick the target with the largest
  intra-target CA RMSD spread from triage_report.csv and split by the
  median; say on the slide that the states are data-derived, not annotated.
