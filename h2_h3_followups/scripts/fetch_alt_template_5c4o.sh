#!/usr/bin/env bash
# Fetches and patches PDB 5C4O (RORyt + allosteric inverse agonist MRL-871,
# Scheepstra et al. Nat Commun 2015) for use as a --template-paths input to
# predict-openfold3. Unlike 4ZJW (gapped, residues 487-500 unresolved),
# 5C4O resolves H12 fully (267-507, no gaps) in a genuinely alternative,
# reoriented pose -- verified sequence-identical to the project construct
# over 479-507 before use (see h2_h3_followups/results/H3_RESULTS_REPORT.md
# appendix and INTEGRITY_NOTES.md).
#
# 5C4O has a real _pdbx_struct_assembly_gen / _struct_asym mismatch bug (same
# class as the one found in the kit's own refs_rorgamma/{3KYT,4ZJW}_clean.cif
# -- stray chains in the assembly row that aren't in struct_asym), which
# makes the deployed apheris-openfold3 CIF reader reject it outright. This
# script downloads the raw file and rewrites only the assembly_gen row to
# intersect with the real chain set -- no atoms/residues/coordinates touched.
set -euo pipefail
OUT=${1:-$HOME/h2_posthoc/alt_templates}
mkdir -p "$OUT"

curl -sL -o "$OUT/5C4O.cif" https://files.rcsb.org/download/5C4O.cif

python3 - "$OUT/5C4O.cif" "$OUT/5C4O_patched.cif" <<'PYEOF'
import sys
src, dst = sys.argv[1], sys.argv[2]
lines = open(src).readlines()
valid = None
for i, l in enumerate(lines):
    if l.strip() == "_struct_asym.id":
        j = i + 1
        while j < len(lines) and lines[j].startswith('_struct_asym.'):
            j += 1
        asym_chains = []
        while j < len(lines) and lines[j].strip() and not lines[j].startswith('_') and not lines[j].startswith('loop_') and not lines[j].startswith('#'):
            asym_chains.append(lines[j].split()[0]); j += 1
        valid = set(asym_chains)
        break
print("struct_asym chains:", sorted(valid))
final = []
i = 0
while i < len(lines):
    l = lines[i]; final.append(l)
    if l.strip().startswith("_pdbx_struct_assembly_gen.asym_id_list"):
        parts = l.split()
        chains_in_row = parts[-1].split(',')
        kept = [c for c in chains_in_row if c in valid]
        if kept != chains_in_row:
            final[-1] = ' '.join(parts[:-1] + [','.join(kept)]) + '\n'
            print("rewrote assembly row:", chains_in_row, "->", kept)
    i += 1
open(dst, 'w').writelines(final)
PYEOF

echo "patched template ready at: $OUT/5C4O.cif (use the ORIGINAL file directly --"
echo "predict-openfold3's own CIF reader tolerated the mismatch in testing; the"
echo "patched copy is provided in case a stricter reader is used downstream)."
