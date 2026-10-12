#!/usr/bin/env bash
# POST-HOC — not pre-registered. Template x depth interaction: does a
# correct-state template bite at depth 128, where the MSA is thin enough that
# the model isn't fully committed? Two jobs: correct-state template (4ZJW,
# H12-out) and the wrong-state control (3KYT, H12-in), both at depth 128,
# same seeds/weights/diffusion-sample count as the pre-registered run.
set -euo pipefail
WORK=$HOME/h2_posthoc
REFS=$HOME/metrics/refs_patched   # assembly-table-patched CIFs; see notes
MSA128=$HOME/h2/msa_128
SEEDS_FLAGS=(--seeds 42 --seeds 43 --seeds 44 --seeds 45 --seeds 46)

declare -A TEMPLATES=( [4zjw]="4ZJW_clean.cif" [3kyt]="3KYT_clean.cif" )

for tag in 4zjw 3kyt; do
  o="$WORK/out/template_${tag}"
  if [ -d "$o" ] && [ -n "$(ls -A "$o" 2>/dev/null)" ]; then
    echo "-- template $tag already has output, skipping"; continue
  fi
  mkdir -p "$o"
  echo "-- template $tag ($(date -u +%H:%M:%SZ))"
  docker run --gpus all --rm \
    --user "$(id -u):$(id -g)" \
    -v /etc/passwd:/etc/passwd:ro -v /etc/group:/etc/group:ro \
    -e HOME=/tmp -e XDG_CACHE_HOME=/tmp/.cache -e TRITON_CACHE_DIR=/tmp/.triton \
    -v "$HOME/h2/in:/in" -v "$MSA128:/msa" -v "$REFS:/templates:ro" -v "$o:/out" \
    "$OF3" module-run predict-openfold3 \
      --input /in/query.json \
      --output /out \
      --weights /app/weights/of3-ob-2025-06-30-174k.pt \
      --msa-dir /msa \
      --template-paths "/templates/${TEMPLATES[$tag]}" \
      "${SEEDS_FLAGS[@]}" \
      --num-diffusion-samples 10 \
    2>&1 | tee "$WORK/logs/template_${tag}.log"
  echo "   done $(date -u +%H:%M:%SZ) -> $o"
done
echo "template jobs complete."
find "$WORK/out/template_4zjw" "$WORK/out/template_3kyt" -name "*.cif" 2>/dev/null | wc -l
