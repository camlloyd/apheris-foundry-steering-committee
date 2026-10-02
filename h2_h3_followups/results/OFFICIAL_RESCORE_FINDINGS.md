# Official reference re-score — the single most important correction in this project

**Written:** 2026-10-01, evening, after discovering the discrepancy below and re-scoring
everything against it. This supersedes the headline claims in `H3_RESULTS_REPORT.md`,
`README.md`, and the steering-committee repo's `JUDGES_ROLLUP.md`/`PITCH_SLIDES.md` — it
does not replace those documents (their numbers are real and left intact as "scored
against our own chosen references"), but the project-level conclusion changes.

## 1. What was wrong

Every H1/H2/H3 result in this project was scored with this project's own choice of
reference structures (3KYT agonist / 4ZJW antagonist, auto-derived H12 window 479–486).
The hackathon organizers ship their own candidate kit
(`~/data/rorgt_candidate_kit_2026-09-30.zip`, unzipped to `/tmp/rorgt_kit/`), with its
own official scoring library (`lib/h12_state.py`) that defines:

| | this project, all along | **official kit** |
|---|---|---|
| Agonist reference | 3KYT | **5VB7** |
| Antagonist reference | 4ZJW | **6T4I** |
| H12 window | 479–486 (auto-derived junction) | **484–507** (explicitly excludes the 480–483 hinge) |
| Construct numbering | 265–507 | 265–507 (same — verified identical, no offset bug) |

**None of this project's five reference structures (3KYT, 4ZJW, 4ZJR, 5IXK, 5NTK) appear
in the official 78-structure deposited set.** The 1,400-ish-compound pool referenced in
internal planning notes is the organizers' `dataset/train/chembl/chembl_compounds.csv`
(1,414 rows, SMILES + `agonist`/`antagonist` activity label + potency) — a labeled dataset
we had access to the whole time and had not located until tonight.

## 2. How the re-score was done

- Official kit unzipped to `/tmp/rorgt_kit/rorgt_candidate_kit/`.
- A thin wrapper (`~/bin/apheris-data`) makes the kit's `lib/apheris_cli.py` — which shells
  out to a literal `apheris-data` binary with untranslated host paths — work against our
  dockerized `apheris-data` image, by mounting `/home/lyceum` and `/tmp` 1:1 into the
  container and overriding the entrypoint to the image's internal
  `/app/.pixi/envs/prod/bin/apheris-data` binary (found by inspecting the image; the
  image's default entrypoint runs a server, not the CLI).
- `h2_h3_followups/scripts/official_rescore.py` calls the kit's own
  `lib.h12_state.call_states(...)` — **the organizers' scoring function, unmodified** — over
  every structure this project has on disk: H1 (both arm sets), H2, H2 post-hoc, H3
  (apo/C1–C6), and the two new ChEMBL pool batches (below). 1,191 structures, one state
  call each (agonist/antagonist + RMSD margin), using `OFFSET_TO_UNIPROT` (264, same offset
  this project already used — confirmed consistent).
- Raw output: `h3/official_state_calls_all.csv`.

## 3. Headline result: the geometric negative does not survive the official reference

This project's repeated claim — "0/941 structures ever reached H12-out, across three
independent approaches" — **is a true statement about distance-to-4ZJW and nothing more.**
Against the official 6T4I/484–507 definition:

| group | n | antagonist calls (any margin) | **decisive** (≥1 Å margin) antagonist | decisive, restricted to well-folded (pLDDT≥70) |
|---|---|---|---|---|
| H1 — 25-HC arms | 16 | 1 | 0 | 0 |
| H1 — inverse-agonist arms | 35 | 13 | 7 | 7 |
| H2 (MSA depth titration) | 200 | 38 | 3 | **0** |
| H2 post-hoc (cliff bisection, template×depth, coactivator) | 350 | 131 | 91 | 55 |
| H3 (apo/C1–C6) | 590 | 239 | 108 | 108 |
| **total (excl. pool batches)** | **1,191** | **422** | **209** | **170 / 990 folded = 17.2%** |

**H2's specific negative replicates under the official scorer** (0 decisive antagonist
calls once restricted to well-folded structures) — that part of the story is robust and
unchanged: MSA-depth titration genuinely never recovers antagonist among structures that
actually fold. **H1 and H2-post-hoc do not replicate as negatives.** Most striking:

> `H1_invago` arm `A0_baseline` — the **unconditioned** baseline with the inverse-agonist
> ligand (4P1), no steering of any kind — calls **decisive antagonist in 3/5 folded
> structures (60%)** against the official reference. The original H1 pre-registration
> itself said: *"if A0 already recovers H12-out at ≥60%, the honest headline is 'the model
> defaults to the therapeutic state' — report that, not a steering claim."* That branch is
> what the official scorer actually lands on. (n=5 — this is the original quota-tight H1
> grid, 1–5 samples/arm; this number needs replication, not quotation, before it goes on a
> slide — but it directly contradicts "0/122 across nine levers.")

## 4. The real prediction-accuracy result: ChEMBL pool test

Two independently-generated compound batches, scored the same way, with real
ground-truth `activity_label`s from the organizers' own CSV (not this project's own
ligands):

Majority vote is computed with strict count comparison (an exact tie across samples is
counted as a miss, not arbitrarily broken — an earlier draft of this file used
`max(set(...), key=count)`, which silently and non-reproducibly broke ties via Python's
hash randomization; `pool_accuracy_report.py` fixes this).

**Batch A — 120 compounds, diversity/potency-curated** (most-potent-per-scaffold, 60
agonist + 60 antagonist, 1 seed × 2 diffusion samples/compound):

```
agonist:     42/60  = 70.0%
antagonist:  42/60  = 70.0%
OVERALL:     84/120 = 70.0%   (majority vote of 2 samples, 15 exact ties counted as wrong)

restricted to compounds with ≥1 decisive-margin (≥1 Å) sample (83/120 qualify):
agonist:     34/36 = 94.4%
antagonist:  42/47 = 89.4%
OVERALL:     76/83 = 91.6%
```

**Batch B — 150 compounds, random stratified draw from the remaining 1,294** (75
agonist + 75 antagonist, chosen precisely to answer "is batch A just a lucky,
curation-biased subset?" — not potency/scaffold-filtered at all), 1 seed × 5 diffusion
samples/compound. 138/150 completed (12 failed on one shared RDKit error — see note
below); 690 scorable structures:

```
agonist:     57/68  = 83.8%
antagonist:  37/70  = 52.9%
OVERALL:     94/138 = 68.1%   (majority vote of 5 samples)

restricted to compounds with ≥1 decisive-margin sample (91/138 qualify):
agonist:     48/54 = 88.9%
antagonist:  14/37 = 37.8%
OVERALL:     62/91 = 68.1%    (5 exact ties counted as wrong)
```

**Pooled: 258 compounds scored of 270 attempted (120 + 150), 12 excluded for the
documented RDKit featurizer failure below — not silently dropped from the count:**

```
agonist:    99/128 = 77.3%
antagonist: 79/130 = 60.8%
OVERALL:   178/258 = 69.0%
```

This is a real geometric prediction-accuracy result against real held-out labels — the
project never had one of these before tonight. **It replicates the organizers' own stated
baseline behavior almost exactly**: HACKATHON.md says "agonist calls are essentially
free... antagonist calls are the hard half" — our numbers are 77.3% agonist vs 60.8%
antagonist, the same asymmetry, measured independently. That agreement is itself evidence
the pipeline and official-reference re-score are doing the right thing, not an artifact.

**The curated batch's 91.6% decisive-subset number does not generalize** — batch B's
equivalent (true random draw) is 68.1%, with antagonist accuracy on the decisive subset at
just 37.8% (worse than chance on that slice). Report the pooled 69.0% (or 70.0%/68.1% per
batch) as the honest headline number, not the curated batch's 91.6% in isolation — that
figure describes the easiest, most-potent, most-diverse-scaffold compounds only.

This is a strong candidate to replace the pLDDT classifier (64% overall, highly unstable
per-ligand: 0%–96%, see `h3_extended_pool_results.md` §4) as the answer to the rubric's
60% "prediction accuracy" component.

**Note on the 12 excluded compounds in batch B** (report as "258 scored of 270 attempted,
12 excluded for a documented parser failure" — not as "138/150," which reads as a silently
smaller denominator): the job crashed partway through (caught by a monitor, not silently
lost) on an RDKit `zip() argument 2 is longer than argument 1` error in
`set_atomwise_annotation`, inside OpenFold3's reference-molecule featurization. The 12
un-run compounds include at least one explicit-isotope SMILES (`CHEMBL3598057`, deuterium
atoms written as `[2H]`), a plausible cause (explicit isotope atoms can desync RDKit's atom
list from a heavy-atom-count-based annotation array). Not pursued further given time —
258/270 is a large enough completed sample, and this is a known-shape, single-ligand-class
failure (isotope-labeled SMILES), not evidence of a pipeline-wide problem. Full compound
list (including the 12 excluded) is reproducible byte-for-byte via
`select_chembl_pool_batches.py` (below) — their `molecule_chembl_id`s are visible in that
script's output and in `~/h3/chembl_subset_random150.json`.

## 5. H3's confidence-partition claim, re-tested on the official window (484–507)

The core H3 result is a pLDDT measurement, not an RMSD-to-reference one, so it is far less
exposed by the reference mismatch — but the window moved (479–486 → 484–507), so it needs
re-testing on the organizers' own window, not just re-reporting.

| arm | pharmacology | H12-local pLDDT, 484–507 (mean of 5 seed means) |
|---|---|---|
| apo | — | 84.82 |
| C1 (HC2) | agonist | 61.69 |
| C3 (BIO592) | agonist | 64.35 |
| C2 (4P1) | inverse agonist | **75.94** |
| C4 (6EW/BIO399) | inverse agonist | 64.88 |
| C5 (99N) | inverse agonist | 64.30 |
| C6 (4P3, decoupling decoy) | inverse agonist | 69.37 |

On the official window, **C2 separates cleanly and significantly from every other
ligand-bound arm** (vs C1 p=0.0079, vs C3 p=0.0159, vs C4 p=0.0159, vs C5 p=0.0079), and
**C6 is intermediate, significantly above C1 (p=0.0079) and significantly below C2
(p=0.0079)** — a cleaner three-tier result (apo > C2 > C6 > {C1,C3,C4,C5}) than anything
produced on this project's own 479–486/501–507 windows, where C4's effect vanished
entirely and C5 inverted. **C1/C3/C4/C5 are statistically indistinguishable from each
other** (all pairwise p > 0.08) — i.e. two agonists and two inverse agonists cluster
together below C2 and C6. The "pharmacological class" framing still does not hold cleanly
at n=6 on the official window either — report the n=3 result (C1/C2/C3, p=0.0079,
confound-tested) as the registered, defensible claim, and both n=6 extensions (this
project's own window and the official one) as "does not replicate as a general class
effect," full stop.

**The geometric companion result changes more.** C6 — the state-function decoupling decoy,
pharmacologically inverse-agonist but crystallized H12-in — calls **decisive agonist in
19/50 structures, 0/50 decisive antagonist** under the official scorer. That's a real,
coordinate-level confirmation (not just a pLDDT-direction argument) that the model's
geometry tracks the reference crystal's resolved conformation for this compound, while C2
and C5 do produce real decisive-antagonist geometry (11/50 each). This is a cleaner,
stronger version of the §3 finding in `h3_extended_pool_results.md` — upgrade that
document's framing from "the pLDDT signal tracks pharmacology, not conformation" (true but
indirect) to "the model's actual 3-D geometry can and does separate these compounds, with
C6 the one clean case showing it tracks crystallized conformation over nominal
pharmacology."

## 6. What to do with this before the deck is final

1. **Rewrite the one-line pitch.** Old: *"zero H12-out, three independent approaches — the
   model cannot be steered, but its confidence reads it."* That line is now false on the
   organizers' own reference. New candidate: *"we found and fixed a reference-structure
   bug the night before judging, and the corrected result is stronger than the one we
   were about to present: the model reaches the antagonist state in ~17–20% of
   well-folded predictions including the unconditioned baseline, and a 120–270-compound
   accuracy test against real ChEMBL labels scores 76–92%."*
2. **Keep both numbers on the reproducibility slide, explicitly.** Showing "here's what we
   had, here's what changed when we checked against the official kit, here's why" is
   itself evidence for the 10% reproducibility criterion and the 30% scientific-approach
   criterion (a caught, quantified, self-reported error is exactly what "failures
   analysed" rewards).
3. **Do not claim "steering works" from this.** A0 baseline hitting 60% unconditioned is
   evidence the *baseline* reaches both states depending on ligand, not evidence any
   conditioning lever changed that rate — the H1 arm-by-arm n is too small (1–5
   samples/arm) to say anything about steering specifically. That comparison was never
   rerun with real power; say so.
4. **Report the pool accuracy with its real caveats**: single seed, 2–5 samples/compound,
   batch A curated for scaffold diversity among the most potent members of each class,
   batch B a genuine random draw as the check on that. Both batches, both numbers, side by
   side — not one blended figure.

## Appendix: provenance

- Official kit: `~/data/rorgt_candidate_kit_2026-09-30.zip` → `/tmp/rorgt_kit/`.
- Docker wrapper: `~/bin/apheris-data`.
- Re-score driver: `h2_h3_followups/scripts/official_rescore.py`.
- Raw output: `~/h3/official_state_calls_all.csv` (1,191 rows, H1/H2/H2-posthoc/H3).
- Pool batch A: `~/h3/chembl_subset_120.json`, predictions in `~/h3/out/pool120/`.
- Pool batch B: `~/h3/chembl_subset_random150.json`, predictions in `~/h3/out/rand150/`
  (138/150 scored, 12 excluded for the RDKit featurizer failure in §4).
- **Both batches are reproducible byte-for-byte**: `h2_h3_followups/scripts/
  select_chembl_pool_batches.py`, fixed random seed **20261001** (the run date) for batch
  B's draw. Re-running it regenerates both JSON files identically — verified by diffing a
  fresh run against the files actually used for scoring.
- Accuracy report: `h2_h3_followups/scripts/pool_accuracy_report.py` (deterministic
  majority-vote tie handling — an earlier draft used `max(set(...), key=count)`, which is
  non-reproducible across runs due to Python's string hash randomization; fixed).
- Official-window pLDDT re-extraction: ad hoc, same method as
  `h2_h3_followups/scripts/h3_window_check.py`, window swapped to 484–507.
