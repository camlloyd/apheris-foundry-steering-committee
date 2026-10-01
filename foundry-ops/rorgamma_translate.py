#!/usr/bin/env python3
"""Translate Lily's foundry_jobs_{invago,25hc} request/params (her best-guess
Hub schema) into the real Nextflow predict-workflow query schema, and patch
stale `_pdbx_struct_assembly_gen` rows in the RCSB-derived "_clean" template
CIFs so the module's CIF parser accepts them.

Pocket `contacts` in her params.json are already 1-based chain-relative
positions (verified against targets_rorgamma.json: author_seqid - 264 equals
the contact value), so residue_index (0-based) = contact_value - 1.
pocket_constraint belongs at the QUERY level (sibling of "chains"), not on
the binder chain -- see build_grid.py's history for why.
"""
import re


def translate_arm(request, params):
    """Convert one arm's (request.json, params.json) pair into
    (query, template_paths, msa_free)."""
    chains = []
    id_to_index = {}
    for s in request["sequences"]:
        if "protein" in s:
            p = s["protein"]
            id_to_index[p["id"]] = len(chains)
            chains.append({"chain_id": p["id"], "polymer_type": "protein", "sequence": p["sequence"]})
        elif "ligand" in s:
            l = s["ligand"]
            id_to_index[l["id"]] = len(chains)
            chains.append({"chain_id": l["id"], "smiles": l["smiles"]})

    query = {"query_id": request["name"], "chains": chains}

    pocket = params.get("pocket")
    if pocket:
        binding_residues = [
            {"chain_index": id_to_index[chain_letter], "residue_index": resnum - 1}
            for chain_letter, resnum in pocket["contacts"]
        ]
        query["pocket_constraint"] = {
            "binder_chain_index": id_to_index[pocket["binder"]],
            "binding_residues": binding_residues,
        }

    template_paths = tuple([params["template_cif"]] if "template_cif" in params else [])
    msa_free = params.get("msa_mode") == "none"
    return query, template_paths, msa_free


def patch_stale_assembly_gen(text, valid_chains):
    """Fix `_pdbx_struct_assembly_gen` data rows that reference chains the
    "_clean" step dropped. A row is filtered to only its valid chains if any
    remain, or dropped entirely if none do (an empty-but-present row would
    leave an invalid, data-less CIF loop_)."""
    valid_chains = set(valid_chains)
    lines = text.split("\n")
    kept = []
    for line in lines:
        m = re.match(r"^(\d+) (\d+) ([A-Za-z,]+)$", line)
        if m:
            present = [c for c in m.group(3).split(",") if c in valid_chains]
            if not present:
                continue  # row is entirely dropped chains
            kept.append(f"{m.group(1)} {m.group(2)} {','.join(present)}")
        else:
            kept.append(line)
    return "\n".join(kept)


def group_arms(arms):
    """Group (arm_id, template_paths, msa_free) tuples by their
    (template_paths, msa_free) signature -- these can't be mixed in one
    Nextflow input dir: template_paths is a batch-wide model param, and
    msa_free/full-MSA arms can't share an msa_dir (the module auto-matches
    an MSA by sequence content, with no per-chain opt-out, and rejects
    multiple a3m files with identical sequence content in one msa_dir).

    Returns a dict keyed by (template_paths, msa_free) -> list of arm_ids,
    in first-seen order, with every input arm_id appearing in exactly one
    group."""
    groups = {}
    for arm_id, template_paths, msa_free in arms:
        key = (template_paths, msa_free)
        groups.setdefault(key, []).append(arm_id)
    return groups
