# Multi-state steering kit — Apheris Foundry Hackathon, 1–2 Oct 2026

Evaluation and experiment-grid tooling for the challenge: *steering co-folding
models toward under-represented conformational states via templating, pocket
conditioning and fine-tuning, and finding out which works best.*

Everything here is tested and runs offline. `python test_kit.py` and
`python test_pipeline.py` both pass on a clean machine with numpy + gemmi.

---

## The 10-minute version

```bash
pip install gemmi numpy            # the only dependencies
python test_kit.py                 # 27 checks, all should pass
bash fetch_refs.sh                 # needs venue wifi
python score_run.py validate --config targets.json --target abl1_imatinib
python make_inputs.py --config targets.json --target abl1_imatinib --out runs/
# ...run the models on runs/abl1_imatinib/*.yaml ...
python score_run.py score --config targets.json --target abl1_imatinib \
    --runs outputs/ --out results.csv
```

The last command prints the table that goes straight on the slide.

---

## Why this framing can win

Three things are established in the literature as of now:

1. **Co-folding models mode-collapse onto the dominant state.** KinConfBench
   (npj Drug Discovery, 2026) benchmarked Boltz-2, Chai-1, Protenix and RFAA
   across 2,225 kinase chains: ~65–75% conformational classification accuracy,
   with severe mode collapse, and apo-state memorisation shared across all
   four models rather than unique to one.

2. **Pose accuracy and state recovery are different things.** The same paper
   found geometric success metrics — ligand RMSD, lDDT-PLI — do *not* correlate
   strongly with recovering the correct conformational state. A model can put
   the ligand roughly right while leaving the protein in the wrong state.

3. **MSA perturbation fixes alternative-state sampling — for apo prediction.**
   Subsampling, clustering and column masking all improve AF3's sampling of
   alternative states, and the AF2 lineage (AF-Cluster, SPEACH_AF, AFsample2,
   SF-Cluster) has been at this for years.

**The gap:** the diagnosis work (1, 2) tests only default inference. The fix
work (3) is about apo structures scored by TM-score, with no ligand in the
picture. **Nobody has put the co-folding-specific levers — templating, pocket
conditioning, MSA depth — head to head on ligand-induced state recovery, and
nobody has asked whether they compound or fight each other.** That is exactly
the question the hackathon asked, and it is genuinely open.

### The specific idea worth pitching: ligand-implied pocket conditioning

Boltz-2's `pocket` constraint declares that a ligand must bind within a zone
defined by a list of residues. For a type II inhibitor, part of the molecule
sits in the allosteric back pocket **that only exists in DFG-out**.

So: don't tell the model to be DFG-out. Tell it *this ligand has to reach
here* — and the only way to satisfy that is to flip the DFG motif.

This is worth pitching because:

- It uses information a medicinal chemist actually has. They know their
  compound is type II; they know which residues line the back pocket. They do
  not have the answer structure.
- It is cheap. No fine-tuning, no GPU-hours, works at inference time.
- It has a clean falsifiable test, with controls that are already built into
  the grid below.

`make_inputs.py` derives the back-pocket residue set from your two reference
structures rather than hardcoding it from a paper — contacts in the target
state minus contacts in the other state. Correct by construction, and easy to
defend when a judge asks where the numbers came from.

### The controls that make it credible

Two arms in the grid exist purely to stop you fooling yourself:

- **`A2_pocket_shared`** conditions on contacts present in *both* states. If
  this flips the state too, you are not steering — you are just perturbing.
- **`A4_template_wrong`** gives a template of the *wrong* state. If the model
  still lands on the target state, it is ignoring your input and you are
  measuring memorisation of the PDB.

If Apheris provide proprietary or post-cutoff structures, use them. Nothing
kills the memorisation objection as cleanly.

---

## What each file does

| File | Purpose |
|---|---|
| `structure_io.py` | PDB/mmCIF loading. gemmi if present, pure-Python PDB fallback if not. |
| `kinase_state.py` | DFG-in/out/inter via the Modi–Dunbrack KinCore D1/D2 criteria. Auto-detects the motifs and prints what it found. |
| `state_recovery.py` | Dual-reference state recovery, plus ligand centroid distance and pocket-contact recall. |
| `make_inputs.py` | Derives the state-specific pocket, then writes all 13 arms as Boltz-2 YAML / OpenFold3 JSON. |
| `score_run.py` | `validate` your references, then `score` a whole grid into a CSV and a summary table. |
| `targets.json` | ABL1 and p38α configs, with honest verified/unverified flags on every PDB ID. |
| `fetch_refs.sh` | Downloads the reference structures. |
| `test_kit.py` | 27 unit checks on synthetic structures with known geometry. |
| `test_pipeline.py` | End-to-end rehearsal: validate → generate → score. |

---

## The two metrics, and why they are separate

**State recovery** (the headline). Superpose the prediction on the *invariant
core*, then measure RMSD over the *state-defining region only*, against both
reference states. Whichever is closer is the verdict; the gap is the margin.

The superpose-on-core trick matters. Global superposition plus global RMSD lets
the shared fold — 95% of the atoms, identical between states — swamp the local
rearrangement you are trying to detect.

A margin under 1 Å means the two states were not really distinguished. The
scorer flags this rather than reporting a confident-looking verdict.

**The DFG order parameter** (the independent check). D1/D2 per KinCore, which
knows nothing about your reference structures. The summary table reports how
often the two metrics agree. If agreement is low, put that on the slide — it
means the state call is shaky, and saying so is a finding, not a weakness.

**Pose accuracy** (reported alongside, never instead). Ligand centroid distance
and pocket-contact recall. Given finding (2) above, an arm that improves pose
while leaving the state wrong is a *result*, not a success.

### Thresholds you should sanity-check

`kinase_state.DFG_THRESHOLDS` holds the D1/D2 cutoffs:

```
DFG-in     D1 < 11.0  and  D2 > 14.0
DFG-out    D1 > 11.0  and  D2 < 14.0
DFG-inter  D1 < 11.0  and  D2 < 11.0
```

Different papers quote these slightly differently and KinCore has revised its
criteria at least once. **Run `score_run.py validate` first.** If your known
DFG-out reference does not come back DFG-out, the thresholds or the motif
assignment are wrong for your target — not the model you are testing.

---

## The experiment grid

Each arm changes exactly one thing from baseline, except the last block.

| Arm | What it tests |
|---|---|
| `A0_baseline` | Default inference. Should fail visibly. |
| `A1_pocket_specific` | **The hypothesis.** Pocket conditioned on state-specific contacts. |
| `A2_pocket_shared` | Control: contacts shared by both states. Should *not* flip. |
| `A3_template_target` | Template of the target state. |
| `A4_template_wrong` | Control: wrong-state template. Should *not* land on target. |
| `A5_msa_free` | Removes the evolutionary prior toward the dominant state. |
| `A6_pocket_plus_template` | Do the levers compound? |
| `A7_pocket_msa_free` | Pocket conditioning without the MSA fighting it. |
| `A8_conflict` | Target pocket vs wrong template. Which lever wins? |
| `B0`–`B3` | The same story in OpenFold3, for the cross-model claim. |

Run the cheap levers first (A0–A5 are all inference-time). Only start
fine-tuning overnight if the inference-time arms are already scored.

---

## Day 1 order of work

| Time | Do this |
|---|---|
| 11:00–12:00 | `test_kit.py`, fetch refs, **validate**. Get the baseline running on one target. Do not touch the grid yet. |
| 12:00–13:00 | Baseline scored on ABL1. You should now be able to say "the model collapses to DFG-in N/M times". That alone is a result. |
| 13:00–15:30 | A1–A5. One target, all five levers. Score as they land. |
| 15:30–16:30 | Second target (p38α) on whichever levers won. This is what turns a quirk into a method. |
| 16:30–17:30 | Combination arms, and the conflict arm if time allows. Fine-tuning goes on overnight *only* if everything above is scored. |
| 17:30–18:00 | Build the figure. One panel: state recovery by arm, with the two controls visible. |

Day 2 is only 2.5 hours and demos start at 09:30, so the figure must exist
before you leave on Day 1.

---

## The pitch

> Co-folding models mode-collapse onto the state they saw most, and the field's
> standard metric — ligand RMSD — doesn't detect it. We measured state recovery
> directly, then asked which steering lever actually fixes it. Conditioning on
> the sub-pocket that only exists in the target state recovered it in N of M
> cases, where the baseline got 0. Conditioning on the shared pocket did not,
> so this is steering, not perturbation.

Then the honest line, which is worth more than an extra percentage point:

> We did not test X. The wrong-template control tells us how much of this is
> memorisation, and here is what it showed.

---

## Sources

- KinConfBench / cofolding kinase conformational states — https://www.nature.com/articles/s44386-026-00068-z
- MSA-perturbing methods enhance AF3 alternative-state sampling — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC13588791/
- Persistent bias in multi-state proteins — https://www.biorxiv.org/content/10.64898/2026.07.10.737860.full.pdf
- Kincore (D1/D2 criteria, state labels) — https://dunbrack.fccc.edu/kincore/
- Modi & Dunbrack, kinase conformation nomenclature — https://pubmed.ncbi.nlm.nih.gov/30867294/
- DFG-out conformational analysis, D1/D2 cutoffs — https://pubs.acs.org/doi/10.1021/jm501603h
- SF-Cluster, frustration-guided MSA subsampling — https://arxiv.org/html/2607.00180v1
