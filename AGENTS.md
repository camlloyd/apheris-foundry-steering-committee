# AGENTS.md

Agent-facing orientation for this repo. Human context and the Day-1 plan live
in [HANDOFF.md](HANDOFF.md) — read that first for the science (target,
hypothesis, experimental grid) and the task split between Science Lead and
Foundry Ops.

## What's in this repo

- `HANDOFF.md` — the Day-1 plan: target/hypothesis, the A0-A8/B3 experimental
  grid, timeline, and the Science-Lead vs Foundry-Ops task split. Treat it as
  the source of truth for what the grid arms mean and why.
- `foundry-ops/` — Foundry-side scaffold for running the grid:
  - `grid/grid_spec.json` — one entry per arm (A0-A8, B3): ligand, conditioning,
    seed count, query_id prefix.
  - `grid/pocket_sets.json`, `grid/sequences.json` — **placeholders**, marked
    `_todo`. Real values (protein sequence, ligand SMILES/CCD, pocket residue
    lists, wrong-template CIF) come from the science side's
    `targets_rorgamma.json` once it's delivered.
  - `grid/build_grid.py` — stdlib-only generator. Groups arms by their
    batch-wide `(seeds, template_paths)` signature (these are **not**
    per-query in the Nextflow predict workflow — confirmed from
    `apheris-openfold3-predict.json`'s `parameters` block), and writes one
    Nextflow-ready input dir + `request.json` per group under
    `grid/runs_input/`, plus `arm_manifest.csv` (query_id -> arm_id -> seeds ->
    conditioning) and the exact `nextflow run` command per group in
    `runs_input/commands.txt`.
  - `submit_and_collect.py` — stdlib-only. After a group's Nextflow run
    finishes, reads its `result.json`, cross-references `arm_manifest.csv`,
    and copies predictions into `preds/<arm_id>/seed_<seed>/`, appending to
    `job_tracker.csv`. Idempotent — safe to re-run as more groups land.
  - `tests/`, `pyproject.toml`, `uv.lock` — pytest suite (31 tests, gated at
    95% coverage) for the two scripts above, managed with
    [uv](https://docs.astral.sh/uv/).
  - `Dockerfile` — reproducible container for running the test suite,
    following [astral's uv Docker guide](https://docs.astral.sh/uv/guides/integration/docker/).

## Running the tests

Locally with uv (installs `.venv` from `uv.lock`, no system pip needed):

```bash
cd foundry-ops
uv run pytest
```

Or in a container, for a host with a broken/missing Python toolchain:

```bash
cd foundry-ops
docker build -t foundry-ops .
docker run --rm foundry-ops
```

CI (`.github/workflows/ci.yml`) runs the same `uv sync --locked` + `uv run
pytest` path via `astral-sh/setup-uv`.

## Where the actual workflow lives

The Nextflow predict pipeline these scripts drive is **not** in this repo —
it's in the sibling `foundry-prototypes-and-demos/foundry-predict-workflow/`
checkout (see that repo's own `skills/` for direct-Docker and prediction
contract details). `commands.txt` prints paths relative to that location.

## Conventions for agents working here

- Don't fill in `_todo` placeholders with guessed data — wait for the real
  `targets_rorgamma.json` handoff, then rerun `build_grid.py`.
- `seeds` and `template_paths` are batch-wide Nextflow model params, not
  per-query — any new arm with a different seed count or template need goes
  in its own group, not patched into an existing `request.json`.
- Keep `build_grid.py`/`submit_and_collect.py` themselves stdlib-only (no
  runtime deps) — they need to run anywhere, including hosts with a broken
  system Python. Test/dev tooling (pytest etc.) is fine to manage via uv in
  `pyproject.toml`/`uv.lock`.
- This is placeholder-driven scaffolding for a single hackathon day — don't
  generalize it beyond the RORγ grid or add speculative config.
