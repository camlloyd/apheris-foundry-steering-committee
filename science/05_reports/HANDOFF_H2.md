# Handoff — H2 MSA-depth titration (ready to paste)

H1 is finished and the answer is a clean negative. H2 is a new, pre-registered
experiment on new data — not a rewrite of H1. Everything you need is in
`h2_package/`. Read `H2_PREREGISTRATION.md` first; it was timestamped
**2026-10-01 13:11 UTC, before any titration job existed**, and must not be
edited once results land.

## What H1 settled

Scored all 41 of your predictions with the kit's own scorer.

- **0/25 inverse-agonist predictions recovered H12-out.** Every arm, every seed,
  every lever — pocket-specific, pocket-shared, target template, wrong template,
  MSA-free, and combinations.
- Baseline was 0/5, so the ≥60 % criterion sat in its strong regime and nothing
  cleared it. H1 is refuted on its own pre-registered terms.
- The two crystal states are **7.78 Å** apart across H12. All 22 structurally
  sound predictions sit **1.55 ± 0.44 Å** from H12-in and **5.91 ± 0.60 Å** from
  H12-out — an ensemble ~18× tighter than the gap, entirely at one pole. Closest
  approach to H12-out: 5.53 Å.
- **Templates did nothing systematic.** No-template vs correct-state template
  0.625 Å; no-template vs wrong template 0.736 Å; correct vs wrong 0.807 Å. All
  at seed-noise level, in no consistent direction.

Two things to know before you reuse any of this:

1. **A-arms and B-arms are the same engine.** A0 vs B0 differ by 0.013 Å over 243
   residues. Lily confirmed everything ran on OpenFold3, so the labels are
   just mislabelled, not a bug — but no slide can claim a cross-model result.
   The only Boltz-2 data point is the HPC A0 smoke run.
2. **The H1 scorer had no quality floor.** A5/A7/B2 (MSA-free) came back at
   **pLDDT 43.5–43.7**, Rg 20.6 Å vs 18.3 Å crystal, 27 Å from *both* references
   — and `score_run.py` called them "H12-in, decisive, margin 8.32" because the
   margin test is purely relative. `score_titration.py` in this package fixes
   that; it rejects anything below pLDDT 70 as `unscorable` and prints it rather
   than dropping it. Worth rescoring the 25hc direction with the same floor —
   its "100 %" currently counts three broken structures as successes.

## H2 in one line

Full MSA folds the domain and locks H12-in (pLDDT 91). Zero MSA releases
everything and destroys the fold (pLDDT 44). **Nobody has looked in between.**
H2 titrates that interval. Method precedent is on the organisers' own slide 12
(Sala/Hildebrand/Meiler, MSA subsampling to bias toward user-defined states).

## Run it

On the VM, with the three image env vars already set:

```bash
cd ~/h2_package
bash run_msa_titration.sh
```

That script: writes the query (identical to H1's `B0_of3_baseline` input, so the
titration is comparable to the H1 baseline by construction) → fetches the full
MSA once via `apheris-msa` → builds the depth ladder by truncating the `.a3m`
→ launches one predict job per depth.

**5 jobs, not 25** — `--seeds` takes an array, so all five seeds run in one job,
and `--num-diffusion-samples 10` gives 50 structures per depth. 250 structures
total. It skips any depth that already has output, so it's safe to re-run after
an interruption.

Depths: 2048, 512, 128, 32, 8. Override with `DEPTHS="..."` if the full MSA
comes back shallower than 2048 — check `~/h2/full_depth.txt` after the fetch and
drop any depth at or above the real total, since those would just duplicate the
full-MSA endpoint.

## Score it

```bash
python3 score_titration.py --work ~/h2 --kit ~/apheris_kit_rorgamma --rgkit ~/h2_package
```

Applies both pre-registered conditions — pLDDT ≥ 70 **and** H12-out with margin
≥ 1 Å — and tests the primary criterion (≥2 of 5 seeds at some depth). Writes
`titration_results.csv` and `titration_summary.md`, and prints
`H2 SUPPORTED` or `H2 NOT SUPPORTED (declared negative)`.

Tested here against the real H1 structures: it correctly passed the five sound
A0 seeds and rejected the pLDDT-43.7 MSA-free run.

## The figure

`figure_h1_collapse.png/.svg` is the H1 result and the H2 motivation in one
image — panel A is the collapse in state space, panel B is the untested MSA
window with the five depths marked. Regenerate any time with:

```bash
python3 figure_h1_collapse.py --kit ~/apheris_kit_rorgamma \
  --runs ~/runs/rorgamma_invago --rgkit ~/h2_package
```

## The 60 % — accuracy metrics (now runnable)

The rubric scores accuracy as H12 RMSD plus the mean of **GDT-HA, LDDT-PLI and
ligand BiSyRMSD**. We had none of the last three. There are now two routes, and
you should do both:

```bash
bash run_metrics.sh probe                # capture the compute-metrics schema
bash run_metrics.sh run ~/h2/out         # Foundry's own numbers
```

`run_metrics.sh probe` saves the schema and prints the one-liner that shows its
flags. The `run` block carries a provisional `--input/--output` invocation — if
it errors, the error names the real flag and that is the only place to fix it.

Independently, and **not blocked on the above**:

```bash
python3 local_metrics.py --refs ~/apheris_kit_rorgamma/refs_rorgamma \
    --runs ~/h2/out --ref-state H12-out \
    --smiles "CNC(=O)c1ccc(c(c1)c2ccc3c(c2)CCCN3C(=O)c4c(cccc4Cl)F)Cl" \
    --out metrics.csv
```

That computes GDT-HA, lDDT, LDDT-PLI and ligand RMSD locally. Already run over
the H1 predictions — mean GDT-HA 71.3, mean lDDT 88.3, and it correctly collapses
on the MSA-free arms (GDT-HA 2.9–5.7, lDDT 38.8), which is an independent
confirmation that those three structures are broken.

Prefer Foundry's numbers if the module produces them, and say in the talk that
both agree. One honest caveat already handled in the script's output: 4P1 has a
single graph automorphism, so the symmetry correction is a no-op — report it as
a ligand RMSD after protein superposition, not as BiSyRMSD.

## The 10 % — reproducibility (now a deliverable)

`WORKFLOW.md` is written: what the pipeline does, how to run it, every setting
with the reason it has that value, the validation done before any prediction was
scored, five named limitations, and how to point it at a held-out structure.

The one thing it needs is your compute numbers, which are mechanical:

```bash
bash collect_compute_stats.sh > compute_usage.md
```

That pulls GPU model, image digests (so the run is pinned, not just named), wall
time per depth from `~/h2/logs/`, and structure counts. Run it **after** the
titration finishes.

## Ask the organisers today

`ORGANISER_QUESTIONS.md` has four, with the reason each one changes what we
build. The important one is whether the held-out structure set is released
before submissions close or scored afterwards — it decides whether this evening
goes on running their data or on hardening the pipeline.

## What goes on the slide either way

> Across five conditioning levers and 25 predictions, OpenFold3 never left the
> agonist state. Its H12 ensemble is 18× tighter than the gap between the two
> crystal states and sits entirely at one pole. It places the agonist at 0.35 Å
> and the inverse agonist anywhere from 0.9 to 12.9 Å, at pLDDT 90 throughout —
> pose and state dissociate, and confidence predicts neither. We then tested
> whether MSA depth, not conditioning, controls state access.

If H2 comes back positive, the depth window is the headline and this is the
setup. If it comes back negative, the headline becomes *"fold stability and
state accessibility cannot be separated on this target"* — which closes the
mechanism the literature says should work, and is the stronger of the two
negatives.
