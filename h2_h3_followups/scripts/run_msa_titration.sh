#!/usr/bin/env bash
# H2 MSA-depth titration — Apheris Foundry, team-3 VM
#
#   bash run_msa_titration.sh
#
# Five depths x five seeds x 10 diffusion samples = 5 GPU jobs, 250 structures.
# Criterion is fixed in H2_PREREGISTRATION.md — read it before running, and do
# not edit it afterwards.
#
# Prereqs on the VM:
#   export OF3=quay.io/apheris/foundry-hackathon:apheris-openfold3-v0.15.1
#   export MSA=quay.io/apheris/foundry-hackathon:apheris-msa-v0.9.0
#   export DATA=quay.io/apheris/foundry-hackathon:apheris-data-v0.26.0
set -euo pipefail

WORK=${WORK:-$HOME/h2}
DEPTHS=${DEPTHS:-"2048 512 128 32 8"}
SEEDS=${SEEDS:-"[42,43,44,45,46]"}
NSAMP=${NSAMP:-10}

: "${OF3:?set OF3 to the openfold3 image}"
: "${MSA:?set MSA to the msa image}"

mkdir -p "$WORK"/{in,msa_full,out,logs}

# ---------------------------------------------------------------- 1. the query
# Inverse-agonist direction: full LBD construct (243 res, 3KYT numbering
# 265-507) + 4P1. Identical to the H1 B0_of3_baseline input, so the titration is
# comparable to the H1 baseline by construction.
if [ ! -f "$WORK/in/query.json" ]; then
cat > "$WORK/in/query.json" <<'JSON'
{
  "queries": [
    {
      "query_id": "H2_msa_titration",
      "chains": [
        {
          "chain_id": "A",
          "polymer_type": "protein",
          "sequence": "ASLTEIEHLVQSVCKSYRETCQLRLEDLLRQRSNIFSREEVTGYQRKSMWEMWERCAHHLTEAIQYVVEFAKRLSGFMELCQNDQIVLLKAGAMEVVLVRMCRAYNADNRTVFFEGKYGGMELFRALGCSELISSIFDFSHSLSALHFSEDEIALYTALVLINAHRPGLQEKRKVEQLQYNLELAFHHHLCKTHRQSILAKLPPKGKLRSLCSQHVERLQIFQHLHPIVVQAAFPPLYKELFS"
        },
        {
          "chain_id": "B",
          "smiles": "CNC(=O)c1ccc(c(c1)c2ccc3c(c2)CCCN3C(=O)c4c(cccc4Cl)F)Cl"
        }
      ]
    }
  ]
}
JSON
echo "wrote $WORK/in/query.json"
fi

# ------------------------------------------------------- 2. fetch the full MSA
if [ -z "$(ls -A "$WORK/msa_full" 2>/dev/null)" ]; then
  echo "== fetching full MSA (CPU, a few minutes) =="
  docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v /etc/passwd:/etc/passwd:ro -v /etc/group:/etc/group:ro \
    -e HOME=/tmp -e XDG_CACHE_HOME=/tmp/.cache \
    -v "$WORK/in:/in" -v "$WORK/msa_full:/out" \
    "$MSA" module-run fetch-msa --input /in --output /out \
      --input-format foundry-v1 \
      --server.type colabfold --server.url https://api.colabfold.com \
    2>&1 | tee "$WORK/logs/fetch_msa.log"
else
  echo "== full MSA already present, skipping fetch =="
fi

A3M=$(find "$WORK/msa_full" -name "*.a3m" | head -1)
if [ -z "$A3M" ]; then
  echo "ERROR: no .a3m produced. Check $WORK/logs/fetch_msa.log — the rest of" >&2
  echo "the titration cannot run without it." >&2
  exit 1
fi
FULL=$(grep -c '^>' "$A3M")
echo "full MSA: $A3M  ($FULL sequences)"
echo "$FULL" > "$WORK/full_depth.txt"

# ------------------------------------------------- 3. build the depth ladder
# Depth = number of sequences retained AFTER the query, taken from the top of
# the file with no re-ranking. This policy is fixed in the pre-registration.
for N in $DEPTHS; do
  d="$WORK/msa_$N"
  mkdir -p "$d"
  base=$(basename "$A3M")
  awk -v keep=$((N+1)) '
    /^>/ { n++ }
    { if (n <= keep) print; else exit }
  ' "$A3M" > "$d/$base"
  got=$(grep -c '^>' "$d/$base")
  echo "  depth $N -> $d/$base ($got sequences incl. query)"
done

# --------------------------------------------------------- 4. launch the grid
# --seeds takes one integer per flag occurrence (confirmed against the deployed
# apheris-openfold3 v0.15.1 CLI contract; a JSON-array string is rejected).
SEEDS_FLAGS=()
for s in $(echo "$SEEDS" | tr -d '[]' | tr ',' ' '); do
  SEEDS_FLAGS+=(--seeds "$s")
done

echo
echo "== launching $(echo $DEPTHS | wc -w) jobs, seeds $SEEDS, $NSAMP diffusion samples each =="
for N in $DEPTHS; do
  if [ "$N" -ge "$FULL" ]; then
    echo "-- depth $N >= full MSA depth ($FULL sequences); skipping, would duplicate the full-MSA endpoint"
    continue
  fi
  o="$WORK/out/depth_$N"
  if [ -d "$o" ] && [ -n "$(ls -A "$o" 2>/dev/null)" ]; then
    echo "-- depth $N already has output, skipping"
    continue
  fi
  mkdir -p "$o"
  echo "-- depth $N  ($(date -u +%H:%M:%SZ))"
  docker run --gpus all --rm \
    --user "$(id -u):$(id -g)" \
    -v /etc/passwd:/etc/passwd:ro -v /etc/group:/etc/group:ro \
    -e HOME=/tmp -e XDG_CACHE_HOME=/tmp/.cache -e TRITON_CACHE_DIR=/tmp/.triton \
    -v "$WORK/in:/in" -v "$WORK/msa_$N:/msa" -v "$o:/out" \
    "$OF3" module-run predict-openfold3 \
      --input /in/query.json \
      --output /out \
      --weights /app/weights/of3-ob-2025-06-30-174k.pt \
      --msa-dir /msa \
      "${SEEDS_FLAGS[@]}" \
      --num-diffusion-samples "$NSAMP" \
    2>&1 | tee "$WORK/logs/depth_$N.log"
  echo "   done $(date -u +%H:%M:%SZ)  -> $o"
done

echo
echo "all depths complete. Structures:"
find "$WORK/out" -name "*.cif" | wc -l
echo
echo "next:  python3 score_titration.py --work $WORK --kit /path/to/apheris_kit_rorgamma"
