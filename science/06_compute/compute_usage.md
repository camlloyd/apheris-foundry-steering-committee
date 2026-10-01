# Compute usage

Generated 2026-10-01T14:01:33Z on `hyd-production-vm-dav0tjf07u6f4nlqffv0`.

## Hardware

```
name, memory.total [MiB], driver_version
NVIDIA RTX PRO 6000 Blackwell Server Edition, 97887 MiB, 595.58.03
```

CPU: 16 cores · RAM: 78 GB

## Images (pinned by digest — this is what makes the run repeatable)

| image | tag | digest |
|---|---|---|
| `quay.io/apheris/foundry-hackathon` | `apheris-openfold3-v0.15.1` | `sha256:2a8b6249315e2ca80a025939111ef266488a6f37d7325888645678b934116ad4` |
| `quay.io/apheris/foundry-hackathon` | `apheris-msa-v0.9.0` | `sha256:f80851ec0263fe04bb80bb654f00bf53302e9db82fc7faa0738c4ebca5d130ae` |
| `quay.io/apheris/foundry-hackathon` | `apheris-data-v0.26.0` | `sha256:e93e8f8cc167b3cd0966298e78a78cd053ced6522488c0525b6f7ba2984e3ae1` |

## Jobs

| stage | wall time | structures |
|---|---|---|
| `depth_128` | see log | 50 |
| `depth_2048` | see log | 0 |
| `depth_32` | see log | 50 |
| `depth_512` | see log | 50 |
| `depth_8` | see log | 50 |
| `driver` | /home/lyceum/h2/out/depth_8 | 0 |
| `fetch_msa` | see log | 0 |
| **total** | | **200** |

## Reproducing this run

```bash
export OF3=quay.io/apheris/foundry-hackathon:apheris-openfold3-v0.15.1
export MSA=quay.io/apheris/foundry-hackathon:apheris-msa-v0.9.0
export DATA=quay.io/apheris/foundry-hackathon:apheris-data-v0.26.0
bash run_msa_titration.sh
python3 score_titration.py --work $HOME/h2 --kit $HOME/apheris_kit_rorgamma --rgkit .
python3 local_metrics.py --refs $HOME/apheris_kit_rorgamma/refs_rorgamma \
    --runs $HOME/h2/out --ref-state H12-out --smiles '<4P1 SMILES>' --out metrics.csv
```

Seeds, depths and diffusion-sample count are fixed in `run_msa_titration.sh`;
the success criterion is fixed in `H2_PREREGISTRATION.md` and was timestamped
before the first job ran.
