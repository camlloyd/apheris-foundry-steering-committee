import json
from pathlib import Path

import pytest

import build_grid

GRID_DIR = Path(build_grid.__file__).resolve().parent


# ---------------------------------------------------------------------------
# ligand_chain
# ---------------------------------------------------------------------------

def test_ligand_chain_none():
    assert build_grid.ligand_chain(None, sequences={}) is None
    assert build_grid.ligand_chain("none", sequences={}) is None


def test_ligand_chain_component_id():
    sequences = {"ligands": {"25-HC": {"component_id": "HC2"}}}
    assert build_grid.ligand_chain("25-HC", sequences) == {"component_id": "HC2"}


def test_ligand_chain_smiles():
    sequences = {"ligands": {"inverse_agonist": {"smiles": "CCO"}}}
    assert build_grid.ligand_chain("inverse_agonist", sequences) == {"smiles": "CCO"}


def test_ligand_chain_prefers_component_id_over_smiles():
    sequences = {"ligands": {"x": {"component_id": "ABC", "smiles": "CCO"}}}
    assert build_grid.ligand_chain("x", sequences) == {"component_id": "ABC"}


# ---------------------------------------------------------------------------
# pocket_constraint
# ---------------------------------------------------------------------------

POCKET_SETS = {
    "binder_chain_index": 1,
    "max_distance": 6.0,
    "active_pocket": {"binding_residues": [{"chain_index": 0, "residue_index": 10}]},
    "inactive_pocket": {"binding_residues": [{"chain_index": 0, "residue_index": 20}]},
    "shared_pocket": {"binding_residues": [{"chain_index": 0, "residue_index": 30}]},
}


def test_pocket_constraint_none_for_unconditioned():
    assert build_grid.pocket_constraint("none", POCKET_SETS) is None
    assert build_grid.pocket_constraint("wrong_template", POCKET_SETS) is None


@pytest.mark.parametrize(
    "conditioning,expected_key",
    [
        ("active_pocket", "active_pocket"),
        ("inactive_pocket", "inactive_pocket"),
        ("shared_pocket", "shared_pocket"),
    ],
)
def test_pocket_constraint_shapes(conditioning, expected_key):
    constraint = build_grid.pocket_constraint(conditioning, POCKET_SETS)
    assert constraint["binder_chain_index"] == 1
    assert constraint["binding_residues"] == POCKET_SETS[expected_key]["binding_residues"]
    assert constraint["max_distance"] == 6.0


def test_pocket_constraint_omits_max_distance_when_absent():
    pocket_sets = {
        "binder_chain_index": 1,
        "max_distance": None,
        "active_pocket": {"binding_residues": [{"chain_index": 0, "residue_index": 1}]},
    }
    constraint = build_grid.pocket_constraint("active_pocket", pocket_sets)
    assert "max_distance" not in constraint


# ---------------------------------------------------------------------------
# build_query
# ---------------------------------------------------------------------------

SEQUENCES = {
    "protein_sequence": "MQIFVK",
    "ligands": {
        "25-HC": {"component_id": "HC2"},
        "inverse_agonist": {"smiles": "CCO"},
    },
    "coactivator_peptide": {"sequence": "LHRLL"},
}


def test_build_query_apo_arm_has_only_protein_chain():
    arm = {"query_id_prefix": "a7_apo", "ligand": "none", "conditioning": "none"}
    query = build_grid.build_query(arm, SEQUENCES, POCKET_SETS)
    assert query["query_id"] == "a7_apo"
    assert len(query["chains"]) == 1
    assert query["chains"][0]["chain_id"] == "A"
    assert query["chains"][0]["sequence"] == "MQIFVK"


def test_build_query_ligand_arm_with_no_conditioning_has_no_pocket_constraint():
    arm = {"query_id_prefix": "a0_baseline_agonist", "ligand": "25-HC", "conditioning": "none"}
    query = build_grid.build_query(arm, SEQUENCES, POCKET_SETS)
    assert len(query["chains"]) == 2
    lig_chain = query["chains"][1]
    assert lig_chain["chain_id"] == "B"
    assert lig_chain["component_id"] == "HC2"
    assert "pocket_constraint" not in lig_chain


def test_build_query_conditioned_arm_carries_pocket_constraint_on_ligand_chain():
    arm = {"query_id_prefix": "a1_active_steer", "ligand": "25-HC", "conditioning": "active_pocket"}
    query = build_grid.build_query(arm, SEQUENCES, POCKET_SETS)
    lig_chain = query["chains"][1]
    assert lig_chain["pocket_constraint"]["binding_residues"] == POCKET_SETS["active_pocket"]["binding_residues"]


def test_build_query_peptide_cofold_adds_third_chain():
    arm = {
        "query_id_prefix": "b3_peptide_cofold",
        "ligand": "25-HC",
        "conditioning": "active_pocket",
        "extra_chain": "coactivator_peptide",
    }
    query = build_grid.build_query(arm, SEQUENCES, POCKET_SETS)
    assert len(query["chains"]) == 3
    assert query["chains"][2]["chain_id"] == "C"
    assert query["chains"][2]["sequence"] == "LHRLL"


# ---------------------------------------------------------------------------
# group_key
# ---------------------------------------------------------------------------

def test_group_key_no_template():
    assert build_grid.group_key({"seeds": 2, "template_paths": None}) == (2, None)


def test_group_key_with_template_is_hashable_tuple():
    key = build_grid.group_key({"seeds": 2, "template_paths": ["ref.cif"]})
    assert key == (2, ("ref.cif",))
    hash(key)  # must be hashable to use as a dict key


def test_group_key_distinguishes_different_seed_counts():
    a = build_grid.group_key({"seeds": 2, "template_paths": None})
    b = build_grid.group_key({"seeds": 5, "template_paths": None})
    assert a != b


# ---------------------------------------------------------------------------
# main() - integration test against the real committed placeholder data
# ---------------------------------------------------------------------------

def test_main_writes_three_groups_from_committed_grid_spec(tmp_path, monkeypatch):
    out_dir = tmp_path / "runs_input"
    monkeypatch.setattr(
        "sys.argv", ["build_grid.py", "--out", str(out_dir)]
    )

    build_grid.main()

    group_dirs = sorted(p for p in out_dir.iterdir() if p.is_dir())
    assert len(group_dirs) == 3

    manifest_path = out_dir / "arm_manifest.csv"
    assert manifest_path.exists()
    rows = manifest_path.read_text().strip().splitlines()
    assert len(rows) == 1 + 10  # header + 10 arms (A0-A8, B3)

    commands_text = (out_dir / "commands.txt").read_text()
    assert commands_text.count("nextflow run") == 3


def test_main_template_group_copies_wrong_template_cif(tmp_path, monkeypatch):
    out_dir = tmp_path / "runs_input"
    monkeypatch.setattr("sys.argv", ["build_grid.py", "--out", str(out_dir)])

    build_grid.main()

    template_group = next(
        p for p in out_dir.iterdir() if p.is_dir() and "template" in p.name
    )
    assert (template_group / "wrong_template_ref.cif").exists()

    request = json.loads((template_group / "request.json").read_text())
    assert len(request["queries"]) == 1
    assert request["queries"][0]["query_id"] == "a6_wrong_template"


def test_main_seeds5_group_has_the_three_inverse_agonist_arms(tmp_path, monkeypatch):
    out_dir = tmp_path / "runs_input"
    monkeypatch.setattr("sys.argv", ["build_grid.py", "--out", str(out_dir)])

    build_grid.main()

    seeds5_group = next(p for p in out_dir.iterdir() if p.name.endswith("seeds5"))
    request = json.loads((seeds5_group / "request.json").read_text())
    query_ids = {q["query_id"] for q in request["queries"]}
    assert query_ids == {"a3_baseline_inverse", "a4_money_arm", "a5_cross_control"}


def test_main_is_safe_to_rerun_and_overwrites_stale_output(tmp_path, monkeypatch):
    out_dir = tmp_path / "runs_input"
    monkeypatch.setattr("sys.argv", ["build_grid.py", "--out", str(out_dir)])

    build_grid.main()
    stale_marker = out_dir / "stale.txt"
    stale_marker.write_text("leftover from a previous run")

    build_grid.main()

    assert not stale_marker.exists()
    assert (out_dir / "arm_manifest.csv").exists()
