#!/usr/bin/env python3
"""Gate 1 (go/no-go): scorer self-consistency on crystals.

Runs BEFORE any prediction is scored. Each reference crystal is fed through
the scorer AS IF it were a prediction, against the pre-registered two-state
references {H12-in: 3KYT, H12-out: 4ZJW} and the pre-registered H12 region
(479-486). A scorer that cannot call the crystals correctly has no business
calling predictions.

Expected calls:
    3KYT -> H12-in   (self; trivial, but catches region-mapping bugs)
    4ZJW -> H12-out  (self)
    5IXK -> H12-out  (alternative inactive ref)
    5NTK -> H12-out  (alternative inactive ref, swung-out H12)
    4ZJR -> H12-in   (THE DECOY: inverse agonist 4P3 with H12 still packed --
                      state-function decoupling. The scorer must be honest
                      about it.)

Exit code 0 = gate passed.
"""
from structure_io import load_structure
from state_recovery2 import state_recovery_generic

REFS = {
    "H12-in": load_structure("refs_rorgamma/3KYT_clean.cif"),
    "H12-out": load_structure("refs_rorgamma/4ZJW_clean.cif"),
}
REGION = {479, 480, 481, 482, 483, 484, 485, 486}

# pdb file -> (chain with the receptor, expected call)
CASES = {
    "3KYT": ("3KYT_clean.cif", None, "H12-in"),
    "4ZJW": ("4ZJW_clean.cif", None, "H12-out"),
    "5IXK": ("5IXK.cif", "B", "H12-out"),
    "5NTK": ("5NTK.cif", "A", "H12-out"),
    "4ZJR": ("4ZJR.cif", "B", "H12-in"),
}


def main():
    n_fail = 0
    hdr = (f"{'crystal':<8s} {'expected':<9s} {'verdict':<9s} "
           f"{'margin_A':>9s} {'decisive':>9s} "
           f"{'rmsd_H12in':>11s} {'rmsd_H12out':>12s}")
    print(hdr)
    print("-" * len(hdr))
    for pdb, (fname, chain, expected) in CASES.items():
        st = load_structure(f"refs_rorgamma/{fname}")
        rec = state_recovery_generic(st, REFS, region=REGION, pred_chain=chain)
        rmsd = {c.ref_name: c.region_rmsd for c in rec.comparisons}
        ok = rec.verdict == expected and rec.decisive
        n_fail += 0 if ok else 1
        print(f"{pdb:<8s} {expected:<9s} {rec.verdict:<9s} "
              f"{rec.margin:>9.2f} {str(rec.decisive):>9s} "
              f"{rmsd['H12-in']:>11.2f} {rmsd['H12-out']:>12.2f}"
              f"   {'ok' if ok else 'FAIL'}")
        for note in rec.notes:
            print(f"         note: {note}")
    print("-" * len(hdr))
    if n_fail:
        print(f"GATE FAILED: {n_fail} mis-called crystal(s). "
              f"Do not score predictions until this is fixed.")
    else:
        print("GATE PASSED: the scorer calls every crystal correctly, "
              "including the 4ZJR decoupling decoy.")
    return n_fail


if __name__ == "__main__":
    raise SystemExit(1 if main() else 0)
