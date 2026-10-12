# raw_runs — the actual working directories behind H2, H2 post-hoc, H3, and the metrics probe

These are `~/h2`, `~/h2_posthoc`, `~/h3`, `~/metrics`, and `~/probe` on the
machine this project ran on, copied in as-is. They are what the scripts in
`../scripts/` actually read and write (`WORK=${WORK:-$HOME/h2}` etc.) — the
processed CSVs and reports in `../results/` are derived from this, not the
other way around.

## What's tracked here, what isn't, and why

1,880 raw OpenFold3 structures (`*_model.cif`, 325 MB) and their matching
per-atom confidence files (`*_scores.json`, 2.1 GB — one float per atom per
structure, so they're disproportionately large relative to the structures
themselves) are **not** tracked. Both regenerate deterministically from the
fixed seeds and pinned image digests already in `../scripts/` and
`../results/h2_compute_usage.md`; carrying 2.4 GB of regenerable bytes in
git history is not worth what it buys.

In their place:

- **`RAW_MODEL_CIF_MANIFEST.csv`** / **`RAW_SCORES_JSON_MANIFEST.csv`** —
  every file that existed at commit time: relative path, size, sha256. A
  bad or partial regeneration is detectable without carrying the bytes.
- **`{h2,h2_posthoc,h3}_mean_plddt.csv`** — the one number actually used
  downstream (per-structure mean pLDDT, extracted from the full per-atom
  arrays before they were left out) kept directly, so `../scripts/` and
  `../results/` analyses that depend on it don't require regenerating
  anything to re-check.

Everything else that's small — driver logs, per-round configs
(`config.json`, `result.json`), query inputs (`query.json`), MSA files
(`.a3m`), reference/template CIFs used as scoring inputs
(`h2_posthoc/alt_templates/`, `metrics/refs_patched/`), existing summary
CSVs — is tracked directly, unchanged.

**Also excluded, deliberately, not even manifested:** `h3/official_rescore_work/`,
`h3/official_rescore_work_rand150/`, `h3/official_rescore_smoke/`, and
`metrics/out_test2/`, `metrics/out_test3/`, `metrics/in_test/` — these are
`apheris-data`'s own intermediate scratch (H12-cropped fragments,
preprocessing output) from running the official re-score, not source data
or a final result. They regenerate as a side effect of re-running
`../scripts/official_rescore.py`, the same way a `/tmp` directory would.

## Regenerating the full set

```bash
# re-run any round per ../scripts/ and ../README.md's own instructions,
# writing output to the same ~/h2 / ~/h2_posthoc / ~/h3 paths these
# manifests were taken from, then verify:
python3 -c "
import csv, hashlib, os
for manifest in ('RAW_MODEL_CIF_MANIFEST.csv', 'RAW_SCORES_JSON_MANIFEST.csv'):
    for row in csv.DictReader(open(manifest)):
        path = os.path.expanduser(os.path.join('~', row['relative_path']))
        h = hashlib.sha256()
        with open(path, 'rb') as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b''):
                h.update(chunk)
        assert h.hexdigest() == row['sha256'], row['relative_path']
"
```
