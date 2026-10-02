#!/usr/bin/env bash
set -euo pipefail
WORK=$HOME/h3
SEEDS_FLAGS=(--seeds 42 --seeds 43 --seeds 44 --seeds 45 --seeds 46)
for arm in apo c1 c2 c3; do
  o="$WORK/out/$arm"
  if [ -d "$o" ] && [ -n "$(ls -A "$o" 2>/dev/null)" ]; then
    echo "-- $arm already has output, skipping"; continue
  fi
  mkdir -p "$o"
  echo "-- $arm ($(date -u +%H:%M:%SZ))"
  docker run --gpus all --rm \
    --user "$(id -u):$(id -g)" \
    -v /etc/passwd:/etc/passwd:ro -v /etc/group:/etc/group:ro \
    -e HOME=/tmp -e XDG_CACHE_HOME=/tmp/.cache -e TRITON_CACHE_DIR=/tmp/.triton \
    -v "$WORK/in_$arm:/in" -v "$HOME/h2/msa_full:/msa" -v "$o:/out" \
    "$OF3" module-run predict-openfold3 \
      --input /in/query.json \
      --output /out \
      --weights /app/weights/of3-ob-2025-06-30-174k.pt \
      --msa-dir /msa \
      "${SEEDS_FLAGS[@]}" \
      --num-diffusion-samples 10 \
    2>&1 | tee "$WORK/logs/$arm.log"
  echo "   done $(date -u +%H:%M:%SZ) -> $o"
done
echo "=== H3 ALL ARMS DONE ==="
find "$WORK/out" -name "*.cif" | wc -l
