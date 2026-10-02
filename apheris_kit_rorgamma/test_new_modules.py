"""
Self-test for the hackathon extension modules (state_recovery2, triage,
foundry_emit, figure). Same philosophy as test_kit.py: synthetic structures
with KNOWN geometry, so a failure here means the code is wrong, not the data.

    python test_new_modules.py
"""
from __future__ import annotations

import copy
import json
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from test_kit import build_kinase
from structure_io import write_pdb
from state_recovery2 import (
    derive_state_region, displacement_profile, map_region, state_recovery_generic,
)
from triage import (
    cluster_conformers, emit_targets_json, group_by_sequence, rank_targets,
    triage_directory,
)
from foundry_emit import (
    arm_levers, arm_to_model_params, collect_results, emit_foundry_jobs,
    match_levers,
)
from figure import make_figure, summary_markdown

PASS, FAIL = "  ok  ", " FAIL "
_failures = []


def check(name, condition, detail=""):
    print(f"[{PASS if condition else FAIL}] {name}" + (f"  ({detail})" if detail else ""))
    if not condition:
        _failures.append(name)


def displace_segment(st, lo, hi, offset):
    """Copy st and rigidly translate residues lo..hi -- a planted state change."""
    out = copy.deepcopy(st)
    for r in out.residues:
        if lo <= r.seqid <= hi:
            for a in list(r.atoms):
                r.atoms[a] = r.atoms[a] + np.asarray(offset, dtype=float)
    return out


# ------------------------------------------------- state_recovery2

def test_region_derivation():
    print("\n--- auto state-region derivation ---")
    ref_a = build_kinase(d1=8.0, d2=16.0, seed=0)
    ref_b = displace_segment(ref_a, 100, 119, [6.0, 0.0, 0.0])

    der = derive_state_region(ref_a, ref_b)
    region = der["region"]
    print(f"    threshold={der['threshold']:.2f} A, region={len(region)} residues, "
          f"segments={[(s[0], s[-1]) for s in der['segments']]}")

    core_of_segment = set(range(103, 117))
    check("planted moving segment is inside the derived region",
          core_of_segment <= region,
          f"missing: {sorted(core_of_segment - region)}")
    check("static core is NOT in the region",
          30 not in region and 170 not in region)
    check("no warnings for a genuinely two-state pair",
          not der["warnings"], f"{der['warnings']}")

    profile = displacement_profile(ref_a, ref_b)
    check("displacement profile covers the chain",
          len(profile) == 200, f"{len(profile)} residues")

    mapped = map_region(region, ref_a, ref_b)
    check("region maps 1:1 onto identical numbering", mapped == region)

    same = derive_state_region(ref_a, copy.deepcopy(ref_a))
    check("identical references trigger the SAME-state warning",
          any("SAME state" in w for w in same["warnings"]),
          f"{same['warnings']}")


def test_generic_recovery():
    print("\n--- target-agnostic state recovery ---")
    ref_a = build_kinase(d1=8.0, d2=16.0, seed=0)
    ref_b = displace_segment(ref_a, 100, 119, [6.0, 0.0, 0.0])

    rng = np.random.default_rng(11)
    pred = displace_segment(ref_a, 100, 119, [6.0, 0.0, 0.0])
    R_true, _ = np.linalg.qr(rng.normal(0, 1, (3, 3)))
    if np.linalg.det(R_true) < 0:
        R_true[:, 0] *= -1
    for r in pred.residues:
        for a in list(r.atoms):
            r.atoms[a] = r.atoms[a] @ R_true + np.array([20.0, 5.0, -9.0]) \
                + rng.normal(0, 0.05, 3)

    res = state_recovery_generic(pred, {"stateA": ref_a, "stateB": ref_b})
    print("   ", res.summary().replace("\n", "\n    "))
    check("verdict is the planted state through a rigid transform",
          res.verdict == "stateB", f"got '{res.verdict}'")
    check("margin is decisive", res.decisive, f"margin={res.margin:.2f} A")
    check("notes record the auto-derived region",
          any("auto-derived" in n for n in res.notes))


# ------------------------------------------------------------- triage

def test_triage():
    print("\n--- pharma data triage ---")
    x_a = build_kinase(d1=8.0, d2=16.0, seed=0, with_ligand=True,
                       ligand_at=[4.0, 2.0, 0.0])
    x_b = displace_segment(
        build_kinase(d1=8.0, d2=16.0, seed=0, with_ligand=True,
                     ligand_at=[4.0, 2.0, 10.0]),
        90, 129, [8.0, 0.0, 0.0])
    rng = np.random.default_rng(3)
    x_b2 = copy.deepcopy(x_b)
    for r in x_b2.residues:
        for a in list(r.atoms):
            r.atoms[a] = r.atoms[a] + rng.normal(0, 0.1, 3)
    y = build_kinase(d1=8.0, d2=16.0, seed=1)
    # a genuinely different target: mutate ~1/3 of the sequence (motif
    # positions are never ALA, so they are untouched)
    for i, r in enumerate(y.residues):
        if r.name == "ALA" and i % 3 == 0:
            r.name = "GLY"

    with tempfile.TemporaryDirectory() as d:
        for name, st in (("x_stateA.pdb", x_a), ("x_stateB.pdb", x_b),
                         ("x_stateB_dup.pdb", x_b2), ("y_other.pdb", y)):
            write_pdb(st, os.path.join(d, name))

        rows = triage_directory(d)
        check("all four structures ingested", len(rows) == 4, f"{len(rows)}")
        check("ligands detected on the X structures",
              all(r.get("max_ligand_atoms", 0) >= 6 for r in rows
                  if "x_" in r["id"]))

        groups = group_by_sequence(rows)
        check("two target groups found", len(groups) == 2, f"{len(groups)}")

        clusters = cluster_conformers(groups, rmsd_cut=2.5)
        x_gid = next(g for g, m in groups.items() if len(m) == 3)
        y_gid = next(g for g, m in groups.items() if len(m) == 1)
        check("target X splits into two conformational states",
              len(clusters[x_gid]) == 2, f"{len(clusters[x_gid])}")
        check("target Y stays one cluster", len(clusters[y_gid]) == 1)

        ranked = rank_targets(groups, clusters)
        check("two-state target ranks above single-state target",
              ranked[0]["target_group"] == x_gid,
              f"top score {ranked[0]['score']} ({ranked[0]['reasons']})")
        check("ranking reasons are printed for the judge",
              "states present" in ranked[0]["reasons"])

        entry = emit_targets_json(groups, clusters, x_gid)
        key = next(iter(entry))
        cfg = entry[key]
        for field in ("target_state", "other_state", "ref_target_state",
                      "ref_other_state", "sequence", "smiles"):
            check(f"targets.json entry has '{field}'", field in cfg)
        check("the two reference states differ",
              cfg["ref_target_state"] != cfg["ref_other_state"])


# -------------------------------------------------------- foundry_emit

def test_foundry_emit():
    print("\n--- Foundry emitter ---")
    mock_schema = {"models": [{"id": "openfold3",
                               "params": {"templates": {}, "seeds": {}}}]}
    keys = match_levers(mock_schema)
    check("schema probe matches the template lever", "template" in keys,
          f"{keys}")
    check("schema probe matches the seeds lever", "seeds" in keys)
    check("pocket lever correctly NOT matched", "pocket" not in keys)

    arm_tpl = {"id": "A3_template_target", "engine": "boltz", "desc": "t",
               "kwargs": {"msa_path": None, "template_path": "refs/tpl.cif"}}
    arm_pkt = {"id": "A1_pocket_specific", "engine": "boltz", "desc": "p",
               "kwargs": {"msa_path": None, "pocket_residues": [10, 11, 12]}}
    arm_msa = {"id": "A5_msa_free", "engine": "boltz", "desc": "m",
               "kwargs": {"msa_path": "empty"}}

    params, unsup = arm_to_model_params(arm_tpl, keys)
    check("template arm maps onto the schema key",
          any("template" in k for k in params), f"{params}")
    check("template arm has nothing unsupported", not unsup)

    params, unsup = arm_to_model_params(arm_pkt, keys)
    check("pocket arm is flagged unsupported when schema lacks it",
          unsup == ["pocket"], f"unsupported={unsup}")

    levers = arm_levers(arm_msa)
    check("MSA-free arm normalises to msa='none'", levers["msa"] == "none")

    with tempfile.TemporaryDirectory() as d:
        runs = os.path.join(d, "runs")
        os.makedirs(runs)
        manifest = {
            "target": "synthetic",
            "config": {"sequence": "ACDEFGHIK", "smiles": "CCO"},
            "arms": [dict(a, input="x.yaml") for a in (arm_tpl, arm_pkt, arm_msa)],
        }
        with open(os.path.join(runs, "manifest.json"), "w") as fh:
            json.dump(manifest, fh)

        out = os.path.join(d, "foundry_jobs")
        written = emit_foundry_jobs(os.path.join(runs, "manifest.json"), out,
                                    schema=mock_schema)
        check("one job per arm written", len(written) == 3)

        with open(os.path.join(out, "A3_template_target.request.json")) as fh:
            req = json.load(fh)
        check("request JSON carries the sequence and SMILES",
              req["sequences"][0]["protein"]["sequence"] == "ACDEFGHIK"
              and req["sequences"][1]["ligand"]["smiles"] == "CCO")

        with open(os.path.join(out, "A5_msa_free.request.json")) as fh:
            req = json.load(fh)
        check("MSA-free request carries msa mode none",
              req["sequences"][0]["protein"].get("msa", {}).get("mode") == "none")

        with open(os.path.join(out, "submit_all.sh")) as fh:
            sub = fh.read()
        check("submit script batches all arms",
              sub.count("workflows run") == 3)
        check("unsupported levers are annotated in the submit script",
              "NOT in the probed schema" in sub)

        jobs = os.path.join(d, "jobs")
        os.makedirs(os.path.join(jobs, "A3_template_target"))
        st = build_kinase(d1=13.0, d2=10.0)
        write_pdb(st, os.path.join(jobs, "A3_template_target", "model_0.pdb"))
        n = collect_results(jobs, runs)
        check("collect_results gathers structure artifacts",
              n == 1 and os.path.exists(
                  os.path.join(runs, "A3_template_target", "model_0.pdb")))


# ------------------------------------------------------------- figure

def test_figure():
    print("\n--- Day-1 figure ---")
    with tempfile.TemporaryDirectory() as d:
        csv_path = os.path.join(d, "results.csv")
        with open(csv_path, "w") as fh:
            fh.write("arm,sample,target_state_recovered,dfg_agrees_with_target,"
                     "decisive,ligand_centroid_dist_A\n")
            cases = [("A0_baseline", 0, 5), ("A1_pocket_specific", 4, 5),
                     ("A2_pocket_shared", 1, 5), ("A3_template_target", 5, 5)]
            for arm, hits, n in cases:
                for i in range(n):
                    fh.write(f"{arm},s{i},{int(i < hits)},{int(i < hits)},1,"
                             f"{1.5 if i < hits else 8.0}\n")

        out_prefix = os.path.join(d, "fig")
        png = make_figure(csv_path, out_prefix, target_state="DFG-out")
        check("figure PNG is written and non-trivial",
              os.path.exists(png) and os.path.getsize(png) > 10000,
              f"{os.path.getsize(png)} bytes")
        check("figure SVG is written",
              os.path.exists(out_prefix + ".svg"))

        md = summary_markdown(csv_path, "DFG-out")
        check("markdown table names the arms",
              "A1_pocket_specific" in md and "A0_baseline" in md)
        check("markdown table has the right recovery rate",
              "80%" in md, "A1 should read 4/5 = 80%")


def test_motif_robustness():
    print("\n--- motif autodetect: real-world failure modes ---")
    # Failure mode 1 (1IEP): DFG-out / alphaC-out breaks the true K-E salt
    # bridge, and a decoy hydrophobic Lys with an in-window Glu sits downstream.
    st = build_kinase(d1=13.0, d2=10.0, seed=0)
    # break the true salt bridge: push the aC-Glu away from the b3-Lys
    st.residues[90].atoms["OE1"] = st.residues[90].atoms["OE1"] + np.array([0, 14.0, 0])
    # decoy: Lys at index 120 preceded by three hydrophobics, Glu at +17
    st.residues[117].name = "VAL"
    st.residues[118].name = "ALA"
    st.residues[119].name = "VAL"
    st.residues[120].name = "LYS"
    st.residues[120].atoms["NZ"] = st.residues[120].atoms["CA"] + np.array([1.5, 2.0, 0.0])
    st.residues[137].name = "GLU"
    st.residues[137].atoms["OE1"] = st.residues[120].atoms["NZ"] + np.array([3.0, 11.0, 0.0])

    from kinase_state import find_kinase_motifs, classify_dfg
    m = find_kinase_motifs(st)
    check("broken salt bridge: true b3-Lys beats the downstream decoy",
          m.b3_lys == 74, f"got {m.b3_lys}, want 74")
    check("broken salt bridge: aC-Glu still the true partner",
          m.ac_glu == 91, f"got {m.ac_glu}, want 91")
    res = classify_dfg(st, m)
    check("broken salt bridge: DFG-out classification survives",
          res.label == "DFG-out", f"got '{res.label}'")

    # Failure mode 2 (1KV2): a motif residue (the DFG Gly) is unmodeled, so
    # the parsed sequence has no contiguous DFG triplet.
    st2 = build_kinase(d1=13.0, d2=10.0, seed=0)
    del st2.residues[162]  # the DFG Gly, seqid 163
    m2 = find_kinase_motifs(st2)
    check("unmodeled DFG Gly: Phe still found via the Asp-Phe anchor",
          m2.dfg_phe == 162, f"got {m2.dfg_phe}, want 162")
    res2 = classify_dfg(st2, m2)
    check("unmodeled DFG Gly: classification still works",
          res2.label == "DFG-out", f"got '{res2.label}'")

    # Overrides: explicit motifs must bypass autodetect entirely.
    from kinase_state import motifs_from_overrides
    st3 = build_kinase(d1=8.0, d2=16.0, seed=0)  # DFG-in geometry
    ov = motifs_from_overrides(st3, {"chain": "A", "b3_lys": 74,
                                     "ac_glu": 91, "ac_glu_p4": 95,
                                     "dfg_asp": 161, "dfg_phe": 162})
    res3 = classify_dfg(st3, ov)
    check("motif_overrides bypass autodetect and classify correctly",
          res3.label == "DFG-in", f"got '{res3.label}'")


if __name__ == "__main__":
    print("=" * 66)
    print("extension modules -- self test")
    print("=" * 66)

    test_region_derivation()
    test_generic_recovery()
    test_triage()
    test_foundry_emit()
    test_figure()
    test_motif_robustness()

    print("\n" + "=" * 66)
    if _failures:
        print(f"{len(_failures)} FAILED: " + ", ".join(_failures))
        sys.exit(1)
    print("all checks passed -- the extension modules are working")
    print("=" * 66)
