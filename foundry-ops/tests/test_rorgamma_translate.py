import re

import rorgamma_translate as rt


# ---------------------------------------------------------------------------
# translate_arm
# ---------------------------------------------------------------------------

def _request(name="A0_baseline", protein_id="A", protein_seq="ASLTEI", ligand_id="B", smiles="CCO"):
    return {
        "name": name,
        "sequences": [
            {"protein": {"id": protein_id, "sequence": protein_seq}},
            {"ligand": {"id": ligand_id, "smiles": smiles}},
        ],
    }


def test_translate_arm_builds_chains_from_protein_and_ligand():
    query, template_paths, msa_free = rt.translate_arm(_request(), {})
    assert query["query_id"] == "A0_baseline"
    assert query["chains"] == [
        {"chain_id": "A", "polymer_type": "protein", "sequence": "ASLTEI"},
        {"chain_id": "B", "smiles": "CCO"},
    ]
    assert template_paths == ()
    assert msa_free is False


def test_translate_arm_contact_offset_is_one_based_to_zero_based():
    params = {"pocket": {"binder": "B", "contacts": [["A", 98], ["A", 216]]}}
    query, _, _ = rt.translate_arm(_request(), params)
    residues = query["pocket_constraint"]["binding_residues"]
    assert residues == [
        {"chain_index": 0, "residue_index": 97},
        {"chain_index": 0, "residue_index": 215},
    ]


def test_translate_arm_pocket_constraint_is_at_query_level_not_chain_level():
    params = {"pocket": {"binder": "B", "contacts": [["A", 1]]}}
    query, _, _ = rt.translate_arm(_request(), params)
    assert "pocket_constraint" in query
    for chain in query["chains"]:
        assert "pocket_constraint" not in chain


def test_translate_arm_no_pocket_means_no_constraint():
    query, _, _ = rt.translate_arm(_request(), {})
    assert "pocket_constraint" not in query


def test_translate_arm_binder_chain_index_resolved_by_letter():
    # protein is "A" (index 0), ligand is "B" (index 1); binder references "B"
    params = {"pocket": {"binder": "B", "contacts": [["A", 1]]}}
    query, _, _ = rt.translate_arm(_request(), params)
    assert query["pocket_constraint"]["binder_chain_index"] == 1


def test_translate_arm_contacts_chain_index_resolved_by_letter():
    params = {"pocket": {"binder": "B", "contacts": [["A", 1]]}}
    query, _, _ = rt.translate_arm(_request(), params)
    assert query["pocket_constraint"]["binding_residues"][0]["chain_index"] == 0


def test_translate_arm_template_cif_becomes_template_paths():
    params = {"template_cif": "refs_rorgamma/4ZJW_clean.cif"}
    _, template_paths, _ = rt.translate_arm(_request(), params)
    assert template_paths == ("refs_rorgamma/4ZJW_clean.cif",)


def test_translate_arm_msa_mode_none_sets_msa_free():
    _, _, msa_free = rt.translate_arm(_request(), {"msa_mode": "none"})
    assert msa_free is True


def test_translate_arm_default_is_not_msa_free():
    _, _, msa_free = rt.translate_arm(_request(), {})
    assert msa_free is False


def test_translate_arm_pocket_plus_template_both_present():
    params = {
        "template_cif": "refs_rorgamma/3KYT_clean.cif",
        "pocket": {"binder": "B", "contacts": [["A", 5]]},
    }
    query, template_paths, msa_free = rt.translate_arm(_request(), params)
    assert template_paths == ("refs_rorgamma/3KYT_clean.cif",)
    assert "pocket_constraint" in query
    assert msa_free is False


# ---------------------------------------------------------------------------
# patch_stale_assembly_gen
# ---------------------------------------------------------------------------

def test_patch_3kyt_filters_mixed_row_instead_of_dropping():
    text = (
        "_pdbx_struct_assembly_gen.assembly_id\n"
        "_pdbx_struct_assembly_gen.oper_expression\n"
        "_pdbx_struct_assembly_gen.asym_id_list\n"
        "1 1 A,B,C,D,E\n"
        "\n"
        "loop_\n"
    )
    patched = rt.patch_stale_assembly_gen(text, valid_chains={"A", "C", "D"})
    assert "1 1 A,C,D" in patched.split("\n")
    assert "1 1 A,B,C,D,E" not in patched


def test_patch_4zjw_drops_entirely_stale_row_keeps_valid_row():
    text = (
        "_pdbx_struct_assembly_gen.assembly_id\n"
        "_pdbx_struct_assembly_gen.oper_expression\n"
        "_pdbx_struct_assembly_gen.asym_id_list\n"
        "1 1 A,C,D\n"
        "2 1 B,E\n"
        "\n"
        "loop_\n"
    )
    patched = rt.patch_stale_assembly_gen(text, valid_chains={"A", "C", "D"})
    lines = patched.split("\n")
    assert "1 1 A,C,D" in lines
    assert not any(line.startswith("2 1") for line in lines)
    assert not any("B" in line or "E" in line for line in lines if re_data_row(line))


def re_data_row(line):
    return bool(re.match(r"^\d+ \d+ [A-Za-z,]+$", line))


def test_patch_row_with_all_invalid_chains_is_dropped():
    text = "1 1 X,Y,Z\n"
    patched = rt.patch_stale_assembly_gen(text, valid_chains={"A", "C", "D"})
    assert patched == ""


def test_patch_row_with_all_valid_chains_passes_through_unchanged():
    text = "1 1 A,C,D\n"
    patched = rt.patch_stale_assembly_gen(text, valid_chains={"A", "C", "D"})
    assert patched == "1 1 A,C,D\n"


def test_patch_non_matching_lines_are_untouched():
    text = "loop_\n_some_other_tag.value\nnot a data row at all\n"
    patched = rt.patch_stale_assembly_gen(text, valid_chains={"A"})
    assert patched == text


# ---------------------------------------------------------------------------
# group_arms
# ---------------------------------------------------------------------------

def test_group_arms_same_signature_lands_in_one_group():
    arms = [
        ("A0", (), False),
        ("A1", (), False),
        ("A2", (), False),
    ]
    groups = rt.group_arms(arms)
    assert len(groups) == 1
    assert groups[((), False)] == ["A0", "A1", "A2"]


def test_group_arms_differing_template_paths_split():
    arms = [
        ("A3", ("3KYT_clean.cif",), False),
        ("A4", ("4ZJW_clean.cif",), False),
    ]
    groups = rt.group_arms(arms)
    assert len(groups) == 2
    assert groups[(("3KYT_clean.cif",), False)] == ["A3"]
    assert groups[(("4ZJW_clean.cif",), False)] == ["A4"]


def test_group_arms_differing_msa_free_split_even_with_same_template():
    arms = [
        ("A3", ("3KYT_clean.cif",), False),
        ("B3", ("3KYT_clean.cif",), True),
    ]
    groups = rt.group_arms(arms)
    assert len(groups) == 2
    assert groups[(("3KYT_clean.cif",), False)] == ["A3"]
    assert groups[(("3KYT_clean.cif",), True)] == ["B3"]


def test_group_arms_no_arm_dropped_or_duplicated():
    arms = [
        ("A0", (), False),
        ("A3", ("3KYT_clean.cif",), False),
        ("A5", (), True),
        ("B1", ("4ZJW_clean.cif",), False),
    ]
    groups = rt.group_arms(arms)
    all_grouped = [arm_id for ids in groups.values() for arm_id in ids]
    assert sorted(all_grouped) == sorted(a for a, _, _ in arms)
