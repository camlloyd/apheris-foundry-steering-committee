#!/usr/bin/env python3
"""Per-residue discriminating analysis: is a pLDDT effect H12-specific, or a
general ligand-contact-region effect?

    python3 per_residue_discriminator.py --work ~/h3 \
        --arms c1:agonist c2:inverse_agonist c3:agonist \
        --out discriminator_report.csv

For each named arm (a subdirectory of --work/out/<arm>/**/*.cif with a ligand
in chain B), computes the real per-structure ligand-protein contact footprint
(any protein residue with a heavy atom within --cutoff of a ligand heavy
atom), pools a stable >=80%-of-structures contact footprint per arm, unions
across arms, and compares seed-level mean pLDDT at:
  (a) the H12 state region (479-486, construct numbering)
  (b) the union footprint EXCLUDING H12 ("non-H12 contact residues")
between every pair of arms with different labels, via an exact permutation
test on 5-vs-5 seed-level means (assumes 5 seeds/arm, same as every H1-H3
run; override with --seeds if different).

Decision rule (pre-specified, not fit after seeing results):
  - if the non-H12 contact effect size (Cohen's d) is comparable to or larger
    than the H12 effect size, and same direction/significance -> CONTACT-
    PERTURBATION model: ligand engagement suppresses confidence broadly, not
    specifically at the state-defining residues.
  - if the H12 effect is present and the non-H12 contact effect is absent or
    much smaller -> H12-SPECIFIC, pharmacology-linked model.

This is exactly the test that arbitrated the H3 follow-up: on RORgamma LBD
with HC2/4P1/BIO592, the result was CONTACT-PERTURBATION (non-H12 effect
size equal to or larger than the H12 effect size) -- see H3_RESULTS_REPORT.md
section 4/5 (confidence-accuracy decoupling) for the write-up. Re-run this on
the deposited pool to check whether that finding generalizes beyond n=3
ligands.
"""

from __future__ import annotations

import argparse
import glob
import itertools
import json
import re
import sys
from collections import Counter, defaultdict

import gemmi
import numpy as np

OFFSET = 264  # 3KYT/construct numbering 265-507 -> model output numbering 1-243
H12 = set(range(479 - OFFSET, 487 - OFFSET))


def ligand_contacts(cif_path: str, cutoff: float) -> set[int]:
    st = gemmi.read_structure(cif_path)
    m = st[0]
    names = [c.name for c in m]
    if "B" not in names:
        return set()
    lig_atoms = np.array([[a.pos.x, a.pos.y, a.pos.z]
                           for res in m["B"] for a in res if not a.is_hydrogen()])
    if len(lig_atoms) == 0:
        return set()
    contacts = set()
    for res in m["A"]:
        prot_atoms = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in res if not a.is_hydrogen()])
        if len(prot_atoms) == 0:
            continue
        d = np.min(np.linalg.norm(prot_atoms[:, None, :] - lig_atoms[None, :, :], axis=-1))
        if d <= cutoff:
            contacts.add(res.seqid.num)
    return contacts


def per_residue_plddt(scores_path: str) -> dict[int, float]:
    d = json.load(open(scores_path))
    atoms = d["atom_predicted_local_distance_difference_test"]
    toks = d["atom_tokens"]
    sums = defaultdict(float)
    counts = defaultdict(int)
    for a, t in zip(atoms, toks):
        sums[t] += a
        counts[t] += 1
    return {t: sums[t] / counts[t] for t in sums}


def seed_of(path: str) -> int:
    return int(re.search(r"seed[_-](\d+)", path).group(1))


def permutation_test(a: np.ndarray, b: np.ndarray):
    pooled = np.concatenate([a, b])
    n_a = len(a)
    obs = np.mean(a) - np.mean(b)
    count = total = 0
    for idx in itertools.combinations(range(len(pooled)), n_a):
        idxs = set(idx)
        pa = pooled[[i in idxs for i in range(len(pooled))]]
        pb = pooled[[i not in idxs for i in range(len(pooled))]]
        if abs(np.mean(pa) - np.mean(pb)) >= abs(obs) - 1e-9:
            count += 1
        total += 1
    sd = np.sqrt((np.var(a, ddof=1) + np.var(b, ddof=1)) / 2)
    cohend = (np.mean(a) - np.mean(b)) / sd if sd > 0 else float("nan")
    gt = sum(1 for x in a for y in b if x > y)
    lt = sum(1 for x in a for y in b if x < y)
    cliffs = (gt - lt) / (len(a) * len(b))
    return obs, count / total, cohend, cliffs


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True)
    ap.add_argument("--arms", nargs="+", required=True,
                     help="arm:label pairs, e.g. c1:agonist c2:inverse_agonist c3:agonist. "
                          "Arms sharing a label are not compared to each other.")
    ap.add_argument("--cutoff", type=float, default=4.5, help="contact distance cutoff, Angstrom")
    ap.add_argument("--footprint-frac", type=float, default=0.8,
                     help="fraction of structures a residue must contact in to join the footprint")
    args = ap.parse_args()

    arm_labels = dict(a.split(":", 1) for a in args.arms)
    arms = list(arm_labels)

    footprint = {}
    plddt_data = {arm: [] for arm in arms}
    for arm in arms:
        cifs = sorted(glob.glob(f"{args.work}/out/{arm}/**/*.cif", recursive=True))
        if not cifs:
            print(f"WARNING: no cifs found for arm '{arm}' under {args.work}/out/{arm}/", file=sys.stderr)
        counter = Counter()
        for cif in cifs:
            contacts = ligand_contacts(cif, args.cutoff)
            counter.update(contacts)
            scores = cif.replace("_model.cif", "_scores.json")
            pr = per_residue_plddt(scores)
            plddt_data[arm].append({"seed": seed_of(cif), "pr": pr})
        n = len(cifs)
        footprint[arm] = {r for r, c in counter.items() if n and c >= args.footprint_frac * n}
        print(f"{arm} ({arm_labels[arm]}): n={n}  footprint (3KYT numbering) = "
              f"{sorted(r + OFFSET for r in footprint[arm])}")

    union_contacts = set().union(*footprint.values()) if footprint else set()
    non_h12 = sorted(union_contacts - H12)
    print(f"\nUnion non-H12 contact residues (3KYT numbering): {[r + OFFSET for r in non_h12]}")

    def seed_means(arm, residue_set):
        by_seed = defaultdict(list)
        for row in plddt_data[arm]:
            vals = [row["pr"][t] for t in residue_set if t in row["pr"]]
            if vals:
                by_seed[row["seed"]].append(np.mean(vals))
        seeds = sorted(by_seed)
        return np.array([np.mean(by_seed[s]) for s in seeds])

    print("\n=== H12 region (479-486) ===")
    h12_means = {arm: seed_means(arm, H12) for arm in arms}
    for arm in arms:
        print(f"  {arm} ({arm_labels[arm]}): {np.round(h12_means[arm], 2)}  mean={np.mean(h12_means[arm]):.2f}")

    print("\n=== Non-H12 ligand-contact residues ===")
    nonh12_means = {arm: seed_means(arm, non_h12) for arm in arms}
    for arm in arms:
        print(f"  {arm} ({arm_labels[arm]}): {np.round(nonh12_means[arm], 2)}  mean={np.mean(nonh12_means[arm]):.2f}")

    print("\n=== Pairwise comparisons across different labels (permutation test) ===")
    pairs = [(a, b) for a, b in itertools.combinations(arms, 2) if arm_labels[a] != arm_labels[b]]
    verdicts = []
    for a, b in pairs:
        obs_h, p_h, d_h, c_h = permutation_test(h12_means[a], h12_means[b])
        obs_n, p_n, d_n, c_n = permutation_test(nonh12_means[a], nonh12_means[b])
        print(f"{a} vs {b}:")
        print(f"  H12:        diff={obs_h:+.2f}  p={p_h:.4f}  Cohen_d={d_h:+.2f}  Cliffs={c_h:+.2f}")
        print(f"  non-H12:    diff={obs_n:+.2f}  p={p_n:.4f}  Cohen_d={d_n:+.2f}  Cliffs={c_n:+.2f}")
        if abs(d_n) >= abs(d_h) * 0.7 and p_n < 0.05:
            verdict = "CONTACT-PERTURBATION (non-H12 effect comparable or larger)"
        elif p_n >= 0.05 and p_h < 0.05:
            verdict = "H12-SPECIFIC (effect absent away from H12)"
        else:
            verdict = "AMBIGUOUS (recheck manually)"
        print(f"  -> {verdict}")
        verdicts.append(verdict)

    print("\n=== Overall verdict ===")
    if verdicts and all("CONTACT-PERTURBATION" in v for v in verdicts):
        print("CONTACT-PERTURBATION model supported across all comparisons.")
    elif verdicts and all("H12-SPECIFIC" in v for v in verdicts):
        print("H12-SPECIFIC model supported across all comparisons.")
    else:
        print("Mixed/ambiguous -- report per-pair verdicts above, do not collapse to one line.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
