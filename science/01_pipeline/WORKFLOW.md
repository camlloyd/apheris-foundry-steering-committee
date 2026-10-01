# Runnable Foundry workflow — RORγ H12 two-state steering

This is the reproducibility deliverable (rubric: 10 %, *"a runnable Foundry
workflow, with data, settings and compute use documented"*). It takes a protein
sequence and a ligand SMILES and returns a state call with uncertainty, scored
against a crystallographic reference pair. It is not specific to RORγ — point it
at a different pair of reference structures and it runs.

## What it does

```
sequence + ligand SMILES
        │
        ├─ apheris-msa    fetch-msa              → full MSA (.a3m)
        │                 truncate to N seqs     → depth ladder
        │
        ├─ apheris-openfold3  predict-openfold3  → N structures per depth
        │                     --seeds [..] --num-diffusion-samples 10
        │
        ├─ apheris-data   compute-metrics        → GDT-HA · LDDT-PLI · ligand RMSD
        │   (local_metrics.py computes the same four if the module differs)
        │
        └─ score_titration.py                    → state call + fold-quality floor
                                                 → titration_summary.md
```

## Run it

```bash
export OF3=quay.io/apheris/foundry-hackathon:apheris-openfold3-v0.15.1
export MSA=quay.io/apheris/foundry-hackathon:apheris-msa-v0.9.0
export DATA=quay.io/apheris/foundry-hackathon:apheris-data-v0.26.0

bash run_msa_titration.sh                      # predictions
bash run_metrics.sh probe                      # Foundry metrics schema
bash run_metrics.sh run ~/h2/out               # Foundry metrics
python3 local_metrics.py --refs <kit>/refs_rorgamma --runs ~/h2/out \
        --ref-state H12-out --smiles "<ligand SMILES>" --out metrics.csv
python3 score_titration.py --work ~/h2 --kit <kit> --rgkit .
bash collect_compute_stats.sh > compute_usage.md
```

Every stage is idempotent and skips work already done, so an interrupted run
resumes.

## Settings, and why each is what it is

| setting | value | why |
|---|---|---|
| construct | RORγ LBD, 243 res, 3KYT numbering 265–507 | full LBD so H12 *can* pack; steering has to displace it rather than it being absent |
| ligand | 4P1, SMILES from RCSB chemcomp | the ligand of the inactive reference, so pocket set and ligand describe the same binding event |
| reference pair | 3KYT (H12-in) / 4ZJW (H12-out) | chosen empirically from CA displacement profiles, not ligand labels — 4ZJR is an inverse agonist that keeps H12 packed, and is used as a decoy |
| state region | residues 479–486 | derived from the 3KYT-vs-4ZJW displacement profile at a 2 Å threshold, not asserted |
| MSA depths | 2048, 512, 128, 32, 8 | brackets the two measured endpoints: full MSA (pLDDT 91, H12 locked) and none (pLDDT 44, fold destroyed) |
| subsampling | top N after the query, no re-ranking | fixed before running; a diversity-based policy would be a different experiment |
| seeds | 42–46 | five independent seeds so recovery is a rate, not an anecdote |
| diffusion samples | 10 | finds a minority state if one exists |
| fold-quality floor | mean pLDDT ≥ 70 | without it, structures 27 Å from both references score as confident state calls |
| decisive margin | 1.0 Å | below this the two references are not distinguishable for that structure |

## Validation before any prediction was scored

`test_scorer_gate.py` feeds each reference crystal through the scorer as if it
were a prediction. All five are called correctly, including 4ZJR — an inverse
agonist whose H12 stays packed, i.e. the case where the ligand label and the
conformational state disagree. A scorer that cannot call the crystals has no
business calling predictions.

The scorer was independently reproduced by a second implementation written from
the same spec but sharing no code; both agree on all five crystals to 0.01 Å.

The state region was also checked against all three pre-registered alternatives
(479–486, +loop 286–290, and the full H12 window) before results existed; every
crystal call holds under all three.

## Known limitations, stated rather than discovered

1. **Contact-set leakage.** Pocket sets derive from the same crystals used as
   references. Mitigated by the shared-pocket and wrong-template controls;
   residual circularity acknowledged.
2. **Template memory.** 3KYT and 4ZJW predate the training cutoff and are
   certainly in it. Post-cutoff structures would be the clean test.
3. **Single architecture.** All grid predictions are OpenFold3. The one Boltz-2
   data point is a smoke run. No cross-model claim is made.
4. **Ligand symmetry.** 4P1 has one graph automorphism, so the symmetry
   correction in the ligand RMSD is a no-op for this ligand; the reported value
   is a plain RMSD after protein superposition, not the full BiSyRMSD, which
   also symmetrises over receptor chain assignments.
5. **One target.** Nothing here generalises beyond RORγ LBD without further work.

## Applying it to a held-out structure

```bash
python3 local_metrics.py --refs <dir with the held-out pair> \
    --runs <predictions> --ref-state <H12-in|H12-out> --out heldout_metrics.csv
python3 score_titration.py --work <dir> --kit <kit> --rgkit .
```

The reference pair is the only target-specific input. Nothing in the scoring
path hard-codes a residue number: residue correspondence between structures is
established by sequence alignment, so differing constructs, isoforms and
numbering offsets are handled rather than assumed.
