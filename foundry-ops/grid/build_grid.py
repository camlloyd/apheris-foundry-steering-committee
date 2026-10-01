#!/usr/bin/env python3
"""Build Nextflow-ready request.json batches for the RORgamma steering grid.

Groups arms by (seeds, template_paths) -- these are batch-wide model params
for foundry-predict-workflow, not per-query -- and writes one input directory
per group, each holding a request.json with one query per (arm, seed-index)
slot, plus a top-level arm_manifest.csv for traceability.
"""
import argparse
import csv
import json
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORKFLOW_DIR = HERE / "../../../foundry-prototypes-and-demos/foundry-predict-workflow"


def load(path):
    with open(path) as f:
        return json.load(f)


def ligand_chain(ligand_name, sequences):
    if ligand_name in (None, "none"):
        return None
    lig = sequences["ligands"][ligand_name]
    if "component_id" in lig:
        return {"component_id": lig["component_id"]}
    return {"smiles": lig["smiles"]}


def pocket_constraint(conditioning, pocket_sets):
    key = {
        "active_pocket": "active_pocket",
        "inactive_pocket": "inactive_pocket",
        "shared_pocket": "shared_pocket",
    }.get(conditioning)
    if key is None:
        return None
    pocket = pocket_sets[key]
    constraint = {
        "binder_chain_index": pocket_sets["binder_chain_index"],
        "binding_residues": pocket["binding_residues"],
    }
    if pocket_sets.get("max_distance") is not None:
        constraint["max_distance"] = pocket_sets["max_distance"]
    return constraint


def build_query(arm, sequences, pocket_sets):
    chains = [
        {
            "chain_id": "A",
            "polymer_type": "protein",
            "sequence": sequences["protein_sequence"],
        }
    ]

    lig = ligand_chain(arm["ligand"], sequences)
    if lig is not None:
        lig_chain = dict(lig)
        lig_chain["chain_id"] = "B"
        constraint = pocket_constraint(arm["conditioning"], pocket_sets)
        if constraint is not None:
            lig_chain["pocket_constraint"] = constraint
        chains.append(lig_chain)

    if arm.get("extra_chain") == "coactivator_peptide":
        chains.append(
            {
                "chain_id": "C",
                "polymer_type": "protein",
                "sequence": sequences["coactivator_peptide"]["sequence"],
            }
        )

    return {"query_id": arm["query_id_prefix"], "chains": chains}


def group_key(arm):
    template_paths = arm["template_paths"]
    return (arm["seeds"], tuple(template_paths) if template_paths else None)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(HERE / "runs_input"))
    args = parser.parse_args()

    grid_spec = load(HERE / "grid_spec.json")
    pocket_sets = load(HERE / "pocket_sets.json")
    sequences = load(HERE / "sequences.json")

    out_root = Path(args.out)
    if out_root.exists():
        shutil.rmtree(out_root)
    out_root.mkdir(parents=True)

    groups = {}
    for arm in grid_spec["arms"]:
        groups.setdefault(group_key(arm), []).append(arm)

    manifest_rows = []
    commands = []

    for i, (key, arms) in enumerate(sorted(groups.items(), key=lambda kv: str(kv[0]))):
        seeds, template_paths = key
        group_name = f"group{i}_seeds{seeds}" + ("_template" if template_paths else "")
        group_dir = out_root / group_name
        group_dir.mkdir(parents=True)

        queries = []
        for arm in arms:
            # seeds is a batch-wide --seeds list applied to every query in
            # the run, so one query per arm here, run with a --seeds list
            # of length `seeds`, gives exactly `seeds` predictions for that
            # arm (arms in a group share the same seed count by
            # construction of group_key).
            query = build_query(arm, sequences=sequences, pocket_sets=pocket_sets)
            queries.append(query)
            manifest_rows.append(
                {
                    "query_id": query["query_id"],
                    "arm_id": arm["arm_id"],
                    "seed_count": seeds,
                    "conditioning": arm["conditioning"],
                    "ligand": arm["ligand"],
                    "group_dir": str(group_dir.relative_to(out_root)),
                }
            )

        request = {"queries": queries}
        (group_dir / "request.json").write_text(json.dumps(request, indent=2) + "\n")

        model_params = {"num_diffusion_samples": 1, "seeds": list(range(42, 42 + seeds))}
        if template_paths:
            model_params["template_paths"] = list(template_paths)
            src = Path(sequences["wrong_template_cif"]["path"])
            if not src.is_absolute():
                src = HERE / src
            shutil.copy(src, group_dir / template_paths[0])

        abs_input = group_dir.resolve()
        abs_output = (group_dir / "output").resolve()
        cmd = (
            f"cd {WORKFLOW_DIR.resolve()} && mkdir -p {abs_output} && "
            f"nextflow run workflow/predict/main.nf -c ../local.config "
            f"--input {abs_input} --output {abs_output} "
            f"--modelParams '{json.dumps(model_params)}'"
        )
        commands.append(cmd)

    with open(out_root / "arm_manifest.csv", "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["query_id", "arm_id", "seed_count", "conditioning", "ligand", "group_dir"]
        )
        writer.writeheader()
        writer.writerows(manifest_rows)

    (out_root / "commands.txt").write_text("\n\n".join(commands) + "\n")

    print(f"Wrote {len(groups)} group(s) under {out_root}/")
    print(f"arm_manifest.csv: {len(manifest_rows)} queries")
    print()
    print("Submission commands (also in runs_input/commands.txt):")
    print()
    for cmd in commands:
        print(cmd)
        print()


if __name__ == "__main__":
    main()
