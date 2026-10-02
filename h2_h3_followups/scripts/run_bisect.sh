#!/usr/bin/env bash
# POST-HOC — not pre-registered. Depths chosen after seeing the H2 result, to
# check the fold-transition region (32-128) more closely. Same query, weights,
# seeds, diffusion-sample count as the pre-registered titration.
set -euo pipefail
WORK=$HOME/h2_posthoc
IN=$HOME/h2/in/query.json
DEPTHS="48 64 96"
SEEDS_FLAGS=(--seeds 42 --seeds 43 --seeds 44 --seeds 45 --seeds 46)

for N in $DEPTHS; do
  o="$WORK/out/depth_$N"
  if [ -d "$o" ] && [ -n "$(ls -A "$o" 2>/dev/null)" ]; then
    echo "-- depth $N already has output, skipping"; continue
  fi
  mkdir -p "$o"
  echo "-- depth $N ($(date -u +%H:%M:%SZ))"
  docker run --gpus all --rm \
    --user "$(id -u):$(id -g)" \
    -v /etc/passwd:/etc/passwd:ro -v /etc/group:/etc/group:ro \
    -e HOME=/tmp -e XDG_CACHE_HOME=/tmp/.cache -e TRITON_CACHE_DIR=/tmp/.triton \
    -v "$HOME/h2/in:/in" -v "$WORK/msa_$N:/msa" -v "$o:/out" \
    "$OF3" module-run predict-openfold3 \
      --input /in/query.json \
      --output /out \
      --weights /app/weights/of3-ob-2025-06-30-174k.pt \
      --msa-dir /msa \
      "${SEEDS_FLAGS[@]}" \
      --num-diffusion-samples 10 \
    2>&1 | tee "$WORK/logs/depth_$N.log"
  echo "   done $(date -u +%H:%M:%SZ) -> $o"
done
echo "bisect jobs complete."
find "$WORK/out" -name "*.cif" | wc -l
