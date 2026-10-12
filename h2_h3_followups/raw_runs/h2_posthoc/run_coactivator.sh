#!/usr/bin/env bash
# POST-HOC — not pre-registered. Coactivator co-fold: full LBD + SRC2-2 LXXLL
# peptide (KHKILHRLLQDS, the exact sequence resolved in chain C of 3KYT) +
# inverse agonist 4P1. If the model docks the peptide into the AF-2 groove
# while the inverse agonist is bound, that structure is biologically
# impossible (H12 must be displaced for the groove to be open) -- a
# functional readout rather than a purely geometric one.
set -euo pipefail
WORK=$HOME/h2_posthoc
SEEDS_FLAGS=(--seeds 42 --seeds 43 --seeds 44 --seeds 45 --seeds 46)

o="$WORK/out/coactivator"
mkdir -p "$o"
echo "-- coactivator co-fold ($(date -u +%H:%M:%SZ))"
docker run --gpus all --rm \
  --user "$(id -u):$(id -g)" \
  -v /etc/passwd:/etc/passwd:ro -v /etc/group:/etc/group:ro \
  -e HOME=/tmp -e XDG_CACHE_HOME=/tmp/.cache -e TRITON_CACHE_DIR=/tmp/.triton \
  -v "$WORK/in_coact:/in" -v "$HOME/h2/msa_full:/msa" -v "$o:/out" \
  "$OF3" module-run predict-openfold3 \
    --input /in/query.json \
    --output /out \
    --weights /app/weights/of3-ob-2025-06-30-174k.pt \
    --msa-dir /msa \
    "${SEEDS_FLAGS[@]}" \
    --num-diffusion-samples 10 \
  2>&1 | tee "$WORK/logs/coactivator.log"
echo "   done $(date -u +%H:%M:%SZ) -> $o"
find "$o" -name "*.cif" | wc -l
