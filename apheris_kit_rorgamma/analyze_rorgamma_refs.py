#!/usr/bin/env python3
"""Pick the RORgamma two-state reference pair empirically.

Loads the downloaded RORgamma LBD references, computes pairwise CA
displacement profiles against the active reference 3KYT, reports the
H12-window signal, and derives the state region for each candidate pair.
Outputs a JSON summary + a PNG profile plot for the 5-minute human review.

Run from the kit root:  python analyze_rorgamma_refs.py
"""
from __future__ import annotations

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from structure_io import load_structure
from state_recovery2 import derive_state_region, displacement_profile

ACTIVE = "3KYT"
CANDIDATES = ["5IXK", "4ZJW", "4ZJR", "5NTK"]
H12_WINDOW = (470, 510)   # RORgamma numbering; agonist lock H479-Y502-F506
OUT_JSON = "refs_rorgamma/ref_pair_analysis.json"
OUT_PNG = "refs_rorgamma/displacement_profiles.png"


def receptor_chain(st):
    """Longest polymer chain = receptor (peptide coactivators are short)."""
    return max(st.chain_ids(), key=lambda c: len(st.chain(c)))


def main():
    refs, meta = {}, {}
    for pdb in [ACTIVE] + CANDIDATES:
        st = load_structure(f"refs_rorgamma/{pdb}.cif")
        ch = receptor_chain(st)
        res = st.chain(ch)
        ligs = [(r.name, r.chain, r.seqid) for r in st.ligands()]
        refs[pdb] = (st, ch)
        meta[pdb] = {
            "chain": ch,
            "n_res": len(res),
            "seq_range": [min(r.seqid for r in res), max(r.seqid for r in res)],
            "ligands": ligs,
        }
        print(f"{pdb}: chain {ch}, {len(res)} res, "
              f"{meta[pdb]['seq_range'][0]}-{meta[pdb]['seq_range'][1]}, "
              f"ligands: {sorted(set(n for n, _, _ in ligs))}")

    summary = {"active": ACTIVE, "h12_window": H12_WINDOW, "candidates": {}}
    fig, ax = plt.subplots(figsize=(9, 4.5))
    # Okabe-Ito colorblind-safe palette + distinct linestyles
    style = {"5IXK": ("#0072B2", "-"), "4ZJW": ("#E69F00", "-"),
             "4ZJR": ("#009E73", "--"), "5NTK": ("#CC79A7", "-")}
    fig, ax = plt.subplots(figsize=(9, 4.5))

    for pdb in CANDIDATES:
        st_a, ch_a = refs[ACTIVE]
        st_b, ch_b = refs[pdb]
        profile = displacement_profile(st_a, st_b, ch_a, ch_b)
        if not profile:
            print(f"{pdb}: NO correspondence with {ACTIVE} -- skipped")
            continue
        der = derive_state_region(st_a, st_b, a_chain=ch_a, b_chain=ch_b)

        disps = np.array([d for _, d in profile])
        seqs = np.array([s for s, _ in profile])
        h12_mask = (seqs >= H12_WINDOW[0]) & (seqs <= H12_WINDOW[1])
        h12_disps = disps[h12_mask]
        region_in_h12 = sum(1 for s in der["region"] if H12_WINDOW[0] <= s <= H12_WINDOW[1])

        cand = {
            "n_aligned": len(profile),
            "max_disp": round(float(disps.max()), 2),
            "median_disp": round(float(np.median(disps)), 2),
            "h12_max_disp": round(float(h12_disps.max()), 2) if h12_disps.size else None,
            "h12_mean_disp": round(float(h12_disps.mean()), 2) if h12_disps.size else None,
            "region_size": len(der["region"]),
            "region_segments": [ [int(x) for x in sg] for sg in der["segments"] ],
            "region_residues_in_h12_window": region_in_h12,
            "threshold": round(der["threshold"], 2) if der["threshold"] else None,
            "warnings": der["warnings"],
        }
        summary["candidates"][pdb] = cand
        print(f"{pdb}: aligned={len(profile)} max={disps.max():.2f} "
              f"H12max={cand['h12_max_disp']} H12mean={cand['h12_mean_disp']} "
              f"region={len(der['region'])} res ({region_in_h12} in H12 window) "
              f"warn={len(der['warnings'])}")
        for w in der["warnings"]:
            print(f"    WARN: {w}")

        color, ls = style.get(pdb, ("0.4", "-"))
        ax.plot(seqs, disps, lw=1.4, label=pdb, alpha=0.9, color=color,
                linestyle=ls)

    ax.axvspan(H12_WINDOW[0], H12_WINDOW[1], color="0.85", zorder=0,
               label="H12 window")
    ax.set_xlabel("Residue number (3KYT numbering)")
    ax.set_ylabel(r"CA displacement vs 3KYT ($\AA$)")
    ax.set_title("RORgamma LBD reference displacement profiles")
    ax.legend(frameon=False, ncol=5)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=200)
    print(f"wrote {OUT_PNG}")

    with open(OUT_JSON, "w") as fh:
        json.dump({"meta": meta, "summary": summary}, fh, indent=2)
    print(f"wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
