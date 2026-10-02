"""
Generate Boltz-2 YAML and OpenFold3 JSON inputs for a steering experiment grid.

The point of this file is that you define the experiment ONCE and it emits
every arm, so on the day you are running and scoring rather than hand-editing
config files.

The important function here is derive_pocket_residues(). Rather than hardcoding
a list of "back pocket" residues from a paper (which would be wrong for half
your targets and impossible to defend to a judge), it works the pocket out from
your own reference structures:

    residues the ligand touches in the TARGET state
        minus
    residues the ligand touches in the OTHER state
        =
    the sub-pocket that only exists in the target state

Conditioning on that set is the interesting experiment: you are not telling the
model "be DFG-out", you are telling it "this ligand has to reach here", and the
only way to satisfy that is to adopt the target state. That is a constraint a
medicinal chemist genuinely has (they know their compound is type II), which is
what makes it a usable method rather than a benchmark trick.

Usage:
    python make_inputs.py --target abl1 --out runs/
    python make_inputs.py --target abl1 --show-pocket      # just print the pocket
"""
from __future__ import annotations

import argparse
import json
import os

from structure_io import load_structure
from state_recovery import pocket_contacts, residue_correspondence


# --------------------------------------------------------------- the pocket

def derive_pocket_residues(
    target_state_ref: str,
    other_state_ref: str,
    cutoff: float = 4.5,
    min_residues: int = 3,
) -> dict:
    """
    Residues that the ligand contacts in the target state but not in the other.

    Returns both the state-specific set (the interesting conditioning signal)
    and the shared set (the boring front-pocket contacts, which make a good
    negative control: conditioning on those should NOT flip the state).

    Numbering is that of target_state_ref.
    """
    tgt = load_structure(target_state_ref)
    oth = load_structure(other_state_ref)

    tgt_contacts = pocket_contacts(tgt, cutoff)
    oth_contacts_local = pocket_contacts(oth, cutoff)

    # map the other structure's contacts into target numbering
    pairs = residue_correspondence(tgt, oth)
    oth_to_tgt = {o.seqid: t.seqid for t, o in pairs}
    oth_contacts = {oth_to_tgt[s] for s in oth_contacts_local if s in oth_to_tgt}

    specific = sorted(tgt_contacts - oth_contacts)
    shared = sorted(tgt_contacts & oth_contacts)

    warnings = []
    if len(specific) < min_residues:
        warnings.append(
            f"only {len(specific)} state-specific contact residues found. "
            f"Either the two states' pockets overlap almost completely (so "
            f"pocket conditioning cannot discriminate them), or one reference "
            f"has no ligand. Check before running the grid."
        )
    if not tgt_contacts:
        warnings.append(f"no ligand contacts at all in {target_state_ref}")
    if not oth_contacts:
        warnings.append(f"no ligand contacts at all in {other_state_ref}")

    return {
        "state_specific": specific,
        "shared": shared,
        "all_target": sorted(tgt_contacts),
        "warnings": warnings,
    }


def map_residues_to_construct(ref_path: str, construct_seq: str,
                              ref_seqids: list) -> list:
    """
    Translate reference-numbered residues into 1-based construct positions,
    which is what Boltz-2 and OpenFold3 expect in their constraint blocks.

    Getting this wrong is the single easiest way to waste an afternoon: your
    pocket constraint silently points at the wrong residues and the arm looks
    like it "didn't work".
    """
    from structure_io import Residue, Structure

    ref = load_structure(ref_path)
    fake = Structure(path="<construct>")
    for i, aa in enumerate(construct_seq, start=1):
        three = next((k for k, v in _ONE_TO_THREE.items() if v == aa), "ALA")
        fake.residues.append(
            Residue(seqid=i, name=three, chain="A",
                    atoms={"CA": __import__("numpy").zeros(3)})
        )

    pairs = residue_correspondence(ref, fake)
    ref_to_construct = {r.seqid: f.seqid for r, f in pairs}

    out, missing = [], []
    for s in ref_seqids:
        if s in ref_to_construct:
            out.append(ref_to_construct[s])
        else:
            missing.append(s)
    if missing:
        print(f"  ! {len(missing)} pocket residues did not map onto the "
              f"construct sequence: {missing[:10]}")
    return out


_ONE_TO_THREE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


# ------------------------------------------------------------ Boltz-2 input

def boltz_yaml(
    sequence: str,
    smiles: str,
    msa_path: str | None = None,
    pocket_residues: list | None = None,
    template_path: str | None = None,
    protein_id: str = "A",
    ligand_id: str = "B",
    max_distance: float | None = None,
) -> str:
    """
    Boltz-2 input YAML. Written out by hand rather than via pyyaml so the file
    is readable and diffable between arms -- when an arm misbehaves you want to
    see at a glance what was different about it.
    """
    lines = ["sequences:"]
    lines += [
        "  - protein:",
        f"      id: {protein_id}",
        f"      sequence: {sequence}",
    ]
    if msa_path == "empty":
        lines.append("      msa: empty        # MSA-free inference")
    elif msa_path:
        lines.append(f"      msa: {msa_path}")

    lines += [
        "  - ligand:",
        f"      id: {ligand_id}",
        f"      smiles: '{smiles}'",
    ]

    if template_path:
        lines += [
            "",
            "templates:",
            f"  - cif: {template_path}",
        ]

    if pocket_residues:
        lines += [
            "",
            "constraints:",
            "  - pocket:",
            f"      binder: {ligand_id}",
            "      contacts:",
        ]
        for r in pocket_residues:
            lines.append(f"        - [{protein_id}, {r}]")
        if max_distance is not None:
            lines.append(f"      max_distance: {max_distance}")

    return "\n".join(lines) + "\n"


# --------------------------------------------------------- OpenFold3 input

def openfold3_query(
    name: str,
    sequence: str,
    smiles: str,
    template_cifs: list | None = None,
    msa_free: bool = False,
    max_rows: int | None = None,
) -> dict:
    """
    OpenFold3 query JSON. Schema follows the CLI's query-json format:
    entities under 'sequences', templates given as CIFs in CIF-direct mode.

    NOTE: OpenFold3's exact key names have moved between releases. Check these
    against the version Apheris gives you on the day before running 40 jobs --
    'run_openfold predict --help' and the repo's example query will settle it
    in two minutes. The structure of the experiment does not change either way.
    """
    protein: dict = {"id": "A", "sequence": sequence}

    if msa_free:
        protein["msa"] = {"mode": "none"}
    elif max_rows is not None:
        protein["msa"] = {"max_rows": max_rows,
                          "max_rows_paired": max(1, max_rows // 2)}

    if template_cifs:
        protein["templates"] = [{"cif": c} for c in template_cifs]

    return {
        "name": name,
        "sequences": [
            {"protein": protein},
            {"ligand": {"id": "B", "smiles": smiles}},
        ],
    }


# ------------------------------------------------------------ the arm grid

def build_arms(cfg: dict, pocket: dict, construct_pocket: dict) -> list:
    """
    The experiment grid. Each arm changes exactly ONE thing from the baseline
    (except the combination arms, which are the point of the last block).

    Negative controls are not optional here. Without the decoy-pocket and
    wrong-template arms you cannot tell steering apart from the model simply
    retrieving a structure it memorised, and that is the first thing a good
    judge will ask.
    """
    seq = cfg["sequence"]
    smiles = cfg["smiles"]
    msa = cfg.get("msa_path")
    tgt_tpl = cfg.get("template_target_state")
    wrong_tpl = cfg.get("template_other_state")

    specific = construct_pocket["state_specific"]
    shared = construct_pocket["shared"]

    arms = [
        # ---- baseline
        dict(id="A0_baseline", engine="boltz",
             desc="default inference, full MSA, no steering",
             kwargs=dict(msa_path=msa)),

        # ---- lever 1: pocket conditioning
        dict(id="A1_pocket_specific", engine="boltz",
             desc="pocket conditioned on state-SPECIFIC contacts (the hypothesis)",
             kwargs=dict(msa_path=msa, pocket_residues=specific)),
        dict(id="A2_pocket_shared", engine="boltz",
             desc="NEGATIVE CONTROL: pocket conditioned on contacts shared by "
                  "both states -- should not flip the state",
             kwargs=dict(msa_path=msa, pocket_residues=shared)),

        # ---- lever 2: templating
        dict(id="A3_template_target", engine="boltz",
             desc="template of the target state",
             kwargs=dict(msa_path=msa, template_path=tgt_tpl)),
        dict(id="A4_template_wrong", engine="boltz",
             desc="NEGATIVE CONTROL: template of the WRONG state -- if this "
                  "also lands on the target state, the model is ignoring the "
                  "template and you are measuring memorisation",
             kwargs=dict(msa_path=msa, template_path=wrong_tpl)),

        # ---- lever 3: MSA depth
        dict(id="A5_msa_free", engine="boltz",
             desc="MSA-free: removes the evolutionary prior that pulls the "
                  "model to the dominant state",
             kwargs=dict(msa_path="empty")),

        # ---- combinations: do the levers add up or fight each other?
        dict(id="A6_pocket_plus_template", engine="boltz",
             desc="state-specific pocket + target template",
             kwargs=dict(msa_path=msa, pocket_residues=specific,
                         template_path=tgt_tpl)),
        dict(id="A7_pocket_msa_free", engine="boltz",
             desc="state-specific pocket, MSA-free",
             kwargs=dict(msa_path="empty", pocket_residues=specific)),
        dict(id="A8_conflict", engine="boltz",
             desc="DIAGNOSTIC: target-state pocket vs WRONG-state template. "
                  "Which lever wins tells you their relative strength",
             kwargs=dict(msa_path=msa, pocket_residues=specific,
                         template_path=wrong_tpl)),

        # ---- OpenFold3 mirror arms, for the cross-model claim
        dict(id="B0_of3_baseline", engine="openfold3",
             desc="OpenFold3 default",
             kwargs=dict()),
        dict(id="B1_of3_template", engine="openfold3",
             desc="OpenFold3 with target-state template (CIF direct mode)",
             kwargs=dict(template_cifs=[tgt_tpl] if tgt_tpl else None)),
        dict(id="B2_of3_msa_free", engine="openfold3",
             desc="OpenFold3 MSA-free",
             kwargs=dict(msa_free=True)),
        dict(id="B3_of3_template_msa_free", engine="openfold3",
             desc="OpenFold3 target template, MSA-free",
             kwargs=dict(template_cifs=[tgt_tpl] if tgt_tpl else None,
                         msa_free=True)),
    ]

    for a in arms:
        if a["engine"] == "boltz":
            a["content"] = boltz_yaml(seq, smiles, **a["kwargs"])
            a["filename"] = f"{a['id']}.yaml"
        else:
            a["content"] = json.dumps(
                openfold3_query(a["id"], seq, smiles, **a["kwargs"]), indent=2
            ) + "\n"
            a["filename"] = f"{a['id']}.json"
    return arms


def resolve_sequence(cfg: dict) -> str:
    """
    The construct sequence to fold.

    Set "sequence": "AUTO" in the config and it is read straight off the
    target-state reference structure, which is what you usually want: you are
    folding the same construct that was crystallised, so the numbering lines
    up and the pocket mapping is exact. Paste a sequence explicitly only if
    you deliberately want a different construct.
    """
    seq = cfg.get("sequence", "AUTO")
    if seq and seq != "AUTO":
        return seq

    ref = load_structure(cfg["ref_target_state"])
    seq = ref.sequence()
    if not seq:
        raise SystemExit(f"could not read a sequence from {cfg['ref_target_state']}")
    print(f"  sequence: AUTO -> {len(seq)} residues read from "
          f"{os.path.basename(cfg['ref_target_state'])}")
    if "X" in seq:
        print(f"  ! sequence contains {seq.count('X')} unknown residues (X). "
              f"Modified or non-standard residues in the crystal structure; "
              f"replace them by hand or fold the UniProt sequence instead.")
    return seq


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True,
                    help="JSON config for one target (see targets.json)")
    ap.add_argument("--target", required=True, help="target key inside the config")
    ap.add_argument("--out", default="runs", help="output directory")
    ap.add_argument("--show-pocket", action="store_true",
                    help="print the derived pocket and stop")
    args = ap.parse_args()

    with open(args.config) as fh:
        cfg = json.load(fh)[args.target]

    print(f"target: {args.target}")
    print(f"  target state : {cfg['target_state']}  ref={cfg['ref_target_state']}")
    print(f"  other state  : {cfg['other_state']}  ref={cfg['ref_other_state']}")

    cfg["sequence"] = resolve_sequence(cfg)

    pocket = derive_pocket_residues(cfg["ref_target_state"], cfg["ref_other_state"])
    for w in pocket["warnings"]:
        print(f"  ! {w}")

    print(f"\n  state-specific contacts ({len(pocket['state_specific'])}): "
          f"{pocket['state_specific']}")
    print(f"  shared contacts        ({len(pocket['shared'])}): "
          f"{pocket['shared']}")

    construct_pocket = {
        "state_specific": map_residues_to_construct(
            cfg["ref_target_state"], cfg["sequence"], pocket["state_specific"]),
        "shared": map_residues_to_construct(
            cfg["ref_target_state"], cfg["sequence"], pocket["shared"]),
    }
    print(f"\n  in construct numbering:")
    print(f"    state-specific: {construct_pocket['state_specific']}")
    print(f"    shared:         {construct_pocket['shared']}")

    if args.show_pocket:
        return

    out_dir = os.path.join(args.out, args.target)
    os.makedirs(out_dir, exist_ok=True)

    arms = build_arms(cfg, pocket, construct_pocket)
    manifest = []
    for a in arms:
        path = os.path.join(out_dir, a["filename"])
        with open(path, "w") as fh:
            fh.write(a["content"])
        manifest.append({"id": a["id"], "engine": a["engine"],
                         "desc": a["desc"], "input": path,
                         "kwargs": a["kwargs"]})
        print(f"  wrote {path:<48s} {a['desc']}")

    with open(os.path.join(out_dir, "manifest.json"), "w") as fh:
        json.dump({"target": args.target, "config": cfg, "arms": manifest},
                  fh, indent=2)
    print(f"\n{len(arms)} arms written to {out_dir}/")


if __name__ == "__main__":
    main()
