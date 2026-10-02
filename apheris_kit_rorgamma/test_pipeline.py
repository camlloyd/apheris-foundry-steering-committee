"""
End-to-end rehearsal of the whole pipeline on synthetic data.

Builds fake reference structures and fake "model outputs", then runs the exact
commands you will run on the day. If this passes, the only unknowns left are
the real structures and the real model -- the plumbing is proven.

    python test_pipeline.py
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from structure_io import write_pdb
from test_kit import build_kinase

HERE = os.path.dirname(os.path.abspath(__file__))


def perturb(st, seed, noise=0.15):
    """A model output: same structure, random pose in space, a little noise."""
    rng = np.random.default_rng(seed)
    R, _ = np.linalg.qr(rng.normal(0, 1, (3, 3)))
    if np.linalg.det(R) < 0:
        R[:, 0] *= -1
    t = rng.normal(0, 20, 3)
    for res in st.residues + st.hetero:
        for a in list(res.atoms):
            res.atoms[a] = res.atoms[a] @ R + t + rng.normal(0, noise, 3)
    return st


def main():
    work = tempfile.mkdtemp(prefix="apheris_pipeline_")
    print(f"working in {work}\n")

    # ---------------------------------------------------- reference structures
    refs = os.path.join(work, "refs")
    os.makedirs(refs)

    # DFG-out reference, ligand reaching into the back pocket
    ref_out = build_kinase(d1=13.0, d2=10.0, seed=0, start_resid=242,
                           with_ligand=True, ligand_at=[6.0, 4.0, 0.0])
    # DFG-in reference, ligand confined to the front pocket
    ref_in = build_kinase(d1=8.0, d2=16.0, seed=0, start_resid=242,
                          with_ligand=True, ligand_at=[10.0, 1.0, 2.0])

    # give both references full side chains near the site so contacts are real
    for st, centre in ((ref_out, [6.0, 4.0, 0.0]), (ref_in, [10.0, 1.0, 2.0])):
        rng = np.random.default_rng(11)
        for res in st.residues:
            ca = res.atoms.get("CA")
            if ca is None:
                continue
            if np.linalg.norm(ca - np.array(centre)) < 14:
                res.atoms.setdefault("CB", ca + rng.normal(0, 1.0, 3))
                res.atoms.setdefault("CG", ca + rng.normal(0, 1.8, 3))

    p_out = os.path.join(refs, "ref_dfg_out.pdb")
    p_in = os.path.join(refs, "ref_dfg_in.pdb")
    write_pdb(ref_out, p_out)
    write_pdb(ref_in, p_in)
    print(f"built references:\n  {p_out}\n  {p_in}\n")

    # ------------------------------------------------------------- the config
    seq = "".join(r.one_letter for r in ref_out.residues)
    cfg = {
        "synthetic": {
            "name": "synthetic test kinase",
            "sequence": seq,
            "smiles": "CC(=O)Nc1ccc(O)cc1",
            "target_state": "DFG-out",
            "other_state": "DFG-in",
            "ref_target_state": p_out,
            "ref_other_state": p_in,
            "template_target_state": p_out,
            "template_other_state": p_in,
            "msa_path": "/path/to/msa.a3m",
            "state_region_seqids": [400, 401, 402, 403, 404, 315, 336],
        }
    }
    cfg_path = os.path.join(work, "targets.json")
    with open(cfg_path, "w") as fh:
        json.dump(cfg, fh, indent=2)

    # ------------------------------------------------------ 1. validate refs
    print("=" * 70)
    print("STEP 1  validate references")
    print("=" * 70)
    r = subprocess.run(
        [sys.executable, os.path.join(HERE, "score_run.py"), "validate",
         "--config", cfg_path, "--target", "synthetic"],
        capture_output=True, text=True, cwd=HERE,
    )
    print(r.stdout)
    if r.stderr:
        print("stderr:", r.stderr[-1500:])
    assert "ok" in r.stdout, "reference validation did not report ok"

    # ------------------------------------------------- 2. generate the inputs
    print("=" * 70)
    print("STEP 2  generate the experiment grid")
    print("=" * 70)
    runs_root = os.path.join(work, "runs")
    r = subprocess.run(
        [sys.executable, os.path.join(HERE, "make_inputs.py"),
         "--config", cfg_path, "--target", "synthetic", "--out", runs_root],
        capture_output=True, text=True, cwd=HERE,
    )
    print(r.stdout)
    if r.stderr:
        print("stderr:", r.stderr[-2000:])
    assert r.returncode == 0, "make_inputs failed"

    arm_dir = os.path.join(runs_root, "synthetic")
    yamls = [f for f in os.listdir(arm_dir) if f.endswith(".yaml")]
    jsons = [f for f in os.listdir(arm_dir) if f.endswith(".json")
             and f != "manifest.json"]
    print(f"generated {len(yamls)} Boltz YAMLs and {len(jsons)} OpenFold3 JSONs")
    assert len(yamls) >= 8 and len(jsons) >= 4

    print("\n--- example: A1_pocket_specific.yaml ---")
    with open(os.path.join(arm_dir, "A1_pocket_specific.yaml")) as fh:
        print(fh.read())

    # --------------------------------------- 3. fake model outputs, then score
    print("=" * 70)
    print("STEP 3  score a simulated set of model outputs")
    print("=" * 70)

    # baseline collapses to DFG-in; the pocket-conditioned arm mostly recovers
    # DFG-out; the shared-pocket control behaves like baseline.
    behaviour = {
        "A0_baseline":         ["in", "in", "in", "in"],
        "A1_pocket_specific":  ["out", "out", "out", "in"],
        "A2_pocket_shared":    ["in", "in", "in", "out"],
        "A3_template_target":  ["out", "out", "out", "out"],
        "A4_template_wrong":   ["in", "in", "in", "in"],
        "A5_msa_free":         ["out", "in", "out", "in"],
    }

    out_root = os.path.join(work, "outputs")
    for arm, states in behaviour.items():
        d = os.path.join(out_root, arm)
        os.makedirs(d, exist_ok=True)
        for i, s in enumerate(states):
            src = (build_kinase(d1=13.0, d2=10.0, seed=0, start_resid=1,
                                with_ligand=True, ligand_at=[6.0, 4.0, 0.0])
                   if s == "out" else
                   build_kinase(d1=8.0, d2=16.0, seed=0, start_resid=1,
                                with_ligand=True, ligand_at=[10.0, 1.0, 2.0]))
            write_pdb(perturb(src, seed=hash((arm, i)) % 10000),
                      os.path.join(d, f"sample_{i}.pdb"))
    shutil.copy(os.path.join(arm_dir, "manifest.json"),
                os.path.join(out_root, "manifest.json"))

    r = subprocess.run(
        [sys.executable, os.path.join(HERE, "score_run.py"), "score",
         "--config", cfg_path, "--target", "synthetic",
         "--runs", out_root, "--out", os.path.join(work, "results.csv")],
        capture_output=True, text=True, cwd=HERE,
    )
    print(r.stdout)
    if r.stderr:
        print("stderr:", r.stderr[-2000:])
    assert r.returncode == 0, "scoring failed"

    # The pipeline must reproduce the behaviour we planted.
    assert "A1_pocket_specific" in r.stdout
    lines = {l.split()[0]: l for l in r.stdout.splitlines()
             if l and l.split() and l.split()[0].startswith("A")}
    base = lines.get("A0_baseline", "")
    pocket = lines.get("A1_pocket_specific", "")
    print("\nplanted signal recovered?")
    print(f"  baseline : {base.strip()}")
    print(f"  pocket   : {pocket.strip()}")
    assert " 0%" in base, "baseline should show 0% target-state recovery"
    assert "75%" in pocket, "pocket arm should show 75% target-state recovery"

    print("\n" + "=" * 70)
    print("PIPELINE OK -- validate, generate, score all work end to end")
    print("=" * 70)
    shutil.rmtree(work)


if __name__ == "__main__":
    main()
