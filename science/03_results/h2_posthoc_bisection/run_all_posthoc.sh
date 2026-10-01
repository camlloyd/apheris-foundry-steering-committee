#!/usr/bin/env bash
set -euo pipefail
cd ~/h2_posthoc
until ! pgrep -f run_bisect.sh >/dev/null; do sleep 15; done
echo "=== bisect done, starting template jobs ==="
bash run_template.sh
echo "=== template done, starting coactivator job ==="
bash run_coactivator.sh
echo "=== ALL POSTHOC JOBS DONE ==="
