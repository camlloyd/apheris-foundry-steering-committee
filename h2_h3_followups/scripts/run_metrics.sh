#!/usr/bin/env bash
# Foundry's own accuracy metrics — the 60 % of the rubric.
#
#   bash run_metrics.sh probe              # capture the compute-metrics schema
#   bash run_metrics.sh run <preds_dir>    # run it over a directory of predictions
#
# The rubric scores accuracy as H12 RMSD plus the mean of GDT-HA, LDDT-PLI and
# ligand BiSyRMSD. Prefer Foundry's numbers if this works — using the platform's
# own metrics module is also reproducibility evidence for the separate 10 %.
# `local_metrics.py` computes the same four locally either way, so we are not
# blocked on this.
set -euo pipefail

: "${DATA:?set DATA=quay.io/apheris/foundry-hackathon:apheris-data-v0.26.0}"
WORK=${WORK:-$HOME/metrics}
KIT=${KIT:-$HOME/apheris_kit_rorgamma}
mkdir -p "$WORK"/{in,out,probe}

case "${1:-probe}" in

probe)
  echo "== apheris-data tasks =="
  docker run --rm "$DATA" module-run --help 2>&1 | tee "$WORK/probe/tasks.txt"
  echo
  for t in compute-metrics preprocess; do
    echo "== $t describe =="
    if docker run --rm "$DATA" module-run "$t" describe --json \
         > "$WORK/probe/$t.json" 2>"$WORK/probe/$t.err"; then
      echo "   -> $WORK/probe/$t.json  ($(wc -l < "$WORK/probe/$t.json") lines)"
    else
      echo "   (no such task, or describe failed — see $WORK/probe/$t.err)"
    fi
  done
  echo
  echo "Now read the flags it accepts:"
  echo "  python3 - <<'EOF'"
  echo "  import json; d=json.load(open('$WORK/probe/compute-metrics.json'))"
  echo "  print('FLAGS:', d['parameters'].get('x-foundry-flags'))"
  echo "  print(json.dumps(d['parameters'].get('properties',{}), indent=1)[:3000])"
  echo "  EOF"
  echo
  echo "Then edit the 'run' block below to match, and re-run:  bash $0 run <preds_dir>"
  ;;

run)
  PREDS=${2:?usage: bash run_metrics.sh run <predictions_dir>}
  # Real schema (from `probe`): --input takes prediction files directly (repeatable,
  # but multiple --input flags are scored as ONE ensemble) or a YAML manifest whose
  # `samples` list gives independent prediction/reference/id triples, scored and
  # reported separately in result.json keyed by `id`. --reference is a single file
  # path, not a directory of references, so we build one manifest entry per predicted
  # structure against the H12-out reference (4ZJW) this H2 run is testing for.
  echo "staging predictions + manifest into $WORK/in"
  rm -rf "$WORK/in"; mkdir -p "$WORK/in"
  PREDS_ABS=$(cd "$PREDS" && pwd)
  # The kit's refs_rorgamma/{3KYT,4ZJW}_clean.cif have a stray
  # _pdbx_struct_assembly_gen row referencing chains dropped during cleaning
  # (B/E not in _struct_asym); apheris-data's CIF reader rejects that outright.
  # Patched copies (assembly row intersected with the real chain set, no atoms/
  # residues touched) live in $WORK/refs_patched — see conversation notes.
  REF_CONTAINER="/kit/4ZJW_clean.cif"
  {
    echo "samples:"
    find "$PREDS_ABS" -name '*.cif' -not -path '*__MACOSX*' | sort | while read -r f; do
      rel=${f#"$PREDS_ABS"/}
      id=$(echo "$rel" | tr '/' '_' | sed 's/\.cif$//')
      echo "  - id: $id"
      echo "    prediction: /preds/$rel"
      echo "    reference: $REF_CONTAINER"
    done
  } > "$WORK/in/manifest.yaml"
  n=$(grep -c '^  - id:' "$WORK/in/manifest.yaml")
  echo "  $n predictions in manifest, reference = $REF_CONTAINER (H12-out, 4ZJW)"

  echo
  echo "== compute-metrics =="
  docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v /etc/passwd:/etc/passwd:ro -v /etc/group:/etc/group:ro \
    -e HOME=/tmp \
    -v "$WORK/in:/in" -v "$WORK/refs_patched:/kit:ro" -v "$PREDS_ABS:/preds:ro" -v "$WORK/out:/out" \
    "$DATA" module-run compute-metrics \
      --input /in/manifest.yaml \
      --output /out \
    2>&1 | tee "$WORK/compute_metrics.log"

  echo
  echo "outputs:"; ls -la "$WORK/out"
  ;;

*)
  sed -n '2,8p' "$0"; exit 1 ;;
esac
