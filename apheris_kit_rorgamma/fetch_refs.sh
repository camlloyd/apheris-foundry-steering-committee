#!/usr/bin/env bash
# Fetch the reference structures. Run this on the venue wifi -- it needs
# ordinary internet access to RCSB.
#
#   bash fetch_refs.sh
#
# Then immediately run:
#   python score_run.py validate --config targets.json --target abl1_imatinib
#
# Do not skip the validate step. It is the only thing standing between you and
# spending the afternoon measuring the wrong residues.

set -euo pipefail
mkdir -p refs
cd refs

# ABL1 -- verified DFG-out holo references
IDS_VERIFIED="1IEP 2HYY 3OXZ 2HIW 6XRG 6XR6"
# ABL1 -- candidates for the DFG-in holo reference, state NOT confirmed
IDS_CHECK="2GQG 1FPU"
# p38 alpha
IDS_P38="1KV2 3HEC 3HEG 1A9U 2EWA"

fetch () {
  for id in $1; do
    if [ -f "${id}.cif" ]; then
      echo "  ${id}.cif already here"
      continue
    fi
    echo "  fetching ${id}"
    if ! curl -sSfL -o "${id}.cif" "https://files.rcsb.org/download/${id}.cif"; then
      echo "  ! ${id} failed -- fetch it by hand from rcsb.org/structure/${id}"
      rm -f "${id}.cif"
    fi
  done
}

echo "ABL1 (verified):"
fetch "$IDS_VERIFIED"
echo "ABL1 (UNVERIFIED -- validate these before trusting them):"
fetch "$IDS_CHECK"
echo "p38 alpha:"
fetch "$IDS_P38"

echo
echo "done. Now run:"
echo "  python score_run.py validate --config targets.json --target abl1_imatinib"
echo "  python score_run.py validate --config targets.json --target p38a_birb796"
