#!/usr/bin/env bash
# Build the compute-usage table the rubric's 10 % asks for.
#
#   bash collect_compute_stats.sh > compute_usage.md
#
# Reads wall times out of the titration logs, counts the structures produced,
# and records the hardware and image digests so the run is reproducible from the
# document alone.
set -euo pipefail
WORK=${WORK:-$HOME/h2}

echo "# Compute usage"
echo
echo "Generated $(date -u +%Y-%m-%dT%H:%M:%SZ) on \`$(hostname)\`."
echo

echo "## Hardware"
echo
if command -v nvidia-smi >/dev/null 2>&1; then
  echo '```'
  nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv
  echo '```'
else
  echo "(nvidia-smi unavailable)"
fi
echo
echo "CPU: $(nproc) cores · RAM: $(free -g 2>/dev/null | awk '/^Mem:/{print $2" GB"}' || echo '?')"
echo

echo "## Images (pinned by digest — this is what makes the run repeatable)"
echo
echo '| image | tag | digest |'
echo '|---|---|---|'
for img in "${OF3:-}" "${MSA:-}" "${DATA:-}"; do
  [ -z "$img" ] && continue
  dig=$(docker image inspect "$img" --format '{{index .RepoDigests 0}}' 2>/dev/null || echo "-")
  echo "| \`${img%%:*}\` | \`${img##*:}\` | \`${dig##*@}\` |"
done
echo

echo "## Jobs"
echo
echo '| stage | wall time | structures |'
echo '|---|---|---|'
total=0
for log in "$WORK"/logs/*.log; do
  [ -e "$log" ] || continue
  name=$(basename "$log" .log)
  # mtime of the log minus mtime of its first line is unreliable; use the
  # start/end markers the runner prints, else fall back to file age span.
  start=$(stat -c %Y "$log" 2>/dev/null || echo 0)
  dir="$WORK/out/$name"
  n=0
  [ -d "$dir" ] && n=$(find "$dir" -name '*.cif' 2>/dev/null | wc -l)
  dur=$(awk '/^-- depth/{s=$NF} /^   done/{print $NF}' "$log" 2>/dev/null | tail -1)
  echo "| \`$name\` | ${dur:-see log} | $n |"
  total=$((total + n))
done
echo "| **total** | | **$total** |"
echo

echo "## Reproducing this run"
echo
echo '```bash'
echo "export OF3=${OF3:-<openfold3 image>}"
echo "export MSA=${MSA:-<msa image>}"
echo "export DATA=${DATA:-<data image>}"
echo "bash run_msa_titration.sh"
echo "python3 score_titration.py --work \$HOME/h2 --kit \$HOME/apheris_kit_rorgamma --rgkit ."
echo "python3 local_metrics.py --refs \$HOME/apheris_kit_rorgamma/refs_rorgamma \\"
echo "    --runs \$HOME/h2/out --ref-state H12-out --smiles '<4P1 SMILES>' --out metrics.csv"
echo '```'
echo
echo "Seeds, depths and diffusion-sample count are fixed in \`run_msa_titration.sh\`;"
echo "the success criterion is fixed in \`H2_PREREGISTRATION.md\` and was timestamped"
echo "before the first job ran."
