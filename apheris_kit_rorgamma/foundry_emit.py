"""
Foundry adapter: run the SAME experiment grid through the Apheris Foundry Hub
instead of (or alongside) raw Boltz-2 / OpenFold3.

Why this file exists. Only model inference runs inside Foundry; everything
else (this code, the scoring, the figure) runs on your machine against
downloaded artifacts. Iteration through a gateway is slow -- submit, queue,
governance, run, download -- so the grid must go up as ONE batch and be
scored as results land. Hand-editing JSON between arms is how afternoons die.

The catch: we don't know the Hub's exact parameter schema until the opening
talks. So this module PROBES it first:

    python foundry_emit.py probe                          # prints what the Hub supports
    python foundry_emit.py emit --config targets.json --target abl1_imatinib \
        --runs runs/abl1_imatinib --out foundry_jobs/     # one request+params per arm
    python foundry_emit.py collect --jobs foundry_jobs/ --runs runs/abl1_imatinib

`emit` reads the arm definitions straight out of the kit's manifest.json, so
the Boltz-2 YAMLs, OpenFold3 JSONs and Foundry requests can never drift apart.

If the schema doesn't expose a lever (e.g. no pocket constraints), that arm is
still written, with the lever listed under "unsupported" in the manifest --
that list is itself a slide: it is concrete product feedback for Apheris,
which is half of what a hosted hackathon is for.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import subprocess

# Candidate schema keys per abstract lever, most likely first. Matched
# case-insensitively, substring-against-flattened-schema-keys. Extend this
# table at 10:35 when the real schema is in front of you -- that is the
# entire adaptation cost, by design.
LEVER_CANDIDATES = {
    "template": ["template_cif", "templates", "template", "template_path"],
    "pocket":   ["pocket_contacts", "pocket_constraint", "pocket", "constraints"],
    "msa":      ["msa_mode", "msa", "max_msa_rows", "msa_rows", "msa_path"],
    "seeds":    ["seeds", "seed", "num_seeds", "n_seeds"],
    "samples":  ["diffusion_samples", "num_samples", "samples", "n_samples"],
}


# ------------------------------------------------------------------ probe

def probe_hub_schema(workflow: str = "predict",
                     cli: str = "apheris-foundry", timeout: int = 60) -> dict:
    """
    Ask the Hub what the workflow's models accept. Returns the parsed JSON,
    or {"_error": ...} on any failure -- probing must never crash the run.
    """
    try:
        proc = subprocess.run(
            [cli, "workflows", "models", "--workflow", workflow, "--json"],
            capture_output=True, text=True, timeout=timeout,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired) as exc:
        return {"_error": f"probe failed: {exc}"}
    if proc.returncode != 0:
        return {"_error": f"probe exited {proc.returncode}: {proc.stderr.strip()}"}
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError:
        return {"_error": f"probe returned non-JSON: {proc.stdout[:200]}"}


def _flatten_keys(obj, prefix="") -> list:
    """All key paths in a nested dict/list schema, lowercase."""
    keys = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            p = f"{prefix}.{k}" if prefix else str(k)
            keys.append(p.lower())
            keys.extend(_flatten_keys(v, p))
    elif isinstance(obj, list):
        for v in obj:
            keys.extend(_flatten_keys(v, prefix))
    return keys


def match_levers(schema: dict) -> dict:
    """Which abstract levers the schema supports, and under which key."""
    flat = _flatten_keys(schema)
    out = {}
    for lever, candidates in LEVER_CANDIDATES.items():
        hit = next((c for c in candidates
                    if any(c in key for key in flat)), None)
        if hit:
            # the actual schema key that matched (first match wins)
            out[lever] = next(k for k in flat if hit in k)
    return out


# ------------------------------------------------------- arm -> parameters

def arm_levers(arm: dict) -> dict:
    """
    Normalise a kit arm (from manifest.json) into abstract levers, whatever
    engine it was written for. msa is 'default' unless the arm removed or
    capped it.
    """
    kw = arm.get("kwargs", {})
    levers = {"template": None, "pocket": None, "msa": "default"}

    if kw.get("template_path"):
        levers["template"] = kw["template_path"]
    elif kw.get("template_cifs"):
        levers["template"] = kw["template_cifs"][0]

    if kw.get("pocket_residues"):
        levers["pocket"] = kw["pocket_residues"]

    if kw.get("msa_path") == "empty" or kw.get("msa_free"):
        levers["msa"] = "none"
    elif kw.get("max_rows") is not None:
        levers["msa"] = {"max_rows": kw["max_rows"]}

    return levers


def arm_to_model_params(arm: dict, schema_keys: dict) -> tuple:
    """
    Map one arm's levers onto the probed schema. Returns (params, unsupported).
    With an empty schema_keys (probe failed), everything is passed through
    under our best-guess key names -- the Hub's schema validation will tell
    you what's wrong, which is what it's for.
    """
    levers = arm_levers(arm)
    params, unsupported = {}, []

    def put(lever, value, guess_key):
        if lever in schema_keys:
            params[schema_keys[lever].split(".")[-1]] = value
        elif not schema_keys:
            params[guess_key] = value
        else:
            unsupported.append(lever)

    if levers["template"]:
        put("template", levers["template"], "template_cif")
    if levers["pocket"]:
        put("pocket", {"binder": "B",
                       "contacts": [["A", r] for r in levers["pocket"]]},
            "pocket")
    if levers["msa"] == "none":
        put("msa", "none", "msa_mode")
    elif isinstance(levers["msa"], dict):
        put("msa", levers["msa"]["max_rows"], "max_msa_rows")

    return params, unsupported


# ------------------------------------------------------------------- emit

def emit_foundry_jobs(manifest_path: str, out_dir: str,
                      schema: dict | None = None,
                      workflow: str = "predict") -> list:
    """
    Write one request.json + params.json per arm, plus submit_all.sh.

    request.json follows the OpenFold3-style query shape (entities under
    'sequences'), the most likely Hub schema -- confirm against
    `apheris-foundry workflows get --workflow predict` on the day and adjust
    here if needed. The experiment structure does not change either way.
    """
    with open(manifest_path) as fh:
        manifest = json.load(fh)
    cfg = manifest["config"]
    arms = manifest["arms"]

    schema_keys = match_levers(schema) if schema else {}
    os.makedirs(out_dir, exist_ok=True)

    written = []
    for arm in arms:
        params, unsupported = arm_to_model_params(arm, schema_keys)

        request = {
            "name": arm["id"],
            "sequences": [
                {"protein": {"id": "A", "sequence": cfg["sequence"]}},
                {"ligand": {"id": "B", "smiles": cfg["smiles"]}},
            ],
        }
        if arm_levers(arm)["msa"] == "none":
            request["sequences"][0]["protein"]["msa"] = {"mode": "none"}

        req_path = os.path.join(out_dir, f"{arm['id']}.request.json")
        par_path = os.path.join(out_dir, f"{arm['id']}.params.json")
        with open(req_path, "w") as fh:
            json.dump(request, fh, indent=2)
        with open(par_path, "w") as fh:
            json.dump(params, fh, indent=2)

        written.append({"id": arm["id"], "desc": arm["desc"],
                        "request": req_path, "params": par_path,
                        "unsupported": unsupported})

    submit = os.path.join(out_dir, "submit_all.sh")
    with open(submit, "w") as fh:
        fh.write("#!/bin/bash\n# Submit the whole grid as ONE batch. Do not "
                 "trickle jobs through a gateway.\nset -e\n")
        for w in written:
            fh.write(f"\n# {w['id']}: {w['desc']}\n")
            if w["unsupported"]:
                fh.write(f"#   ! levers NOT in the probed schema: "
                         f"{', '.join(w['unsupported'])}\n")
            fh.write(f"apheris-foundry workflows run --workflow {workflow} "
                     f"--input {w['request']} --model-params @{w['params']}\n")
    os.chmod(submit, 0o755)

    with open(os.path.join(out_dir, "foundry_manifest.json"), "w") as fh:
        json.dump({"source_manifest": manifest_path, "workflow": workflow,
                   "schema_keys_matched": schema_keys, "arms": written},
                  fh, indent=2)

    n_unsup = sum(1 for w in written if w["unsupported"])
    print(f"wrote {len(written)} Foundry jobs to {out_dir}/"
          + (f" ({n_unsup} arms have levers the schema doesn't advertise)"
             if n_unsup else ""))
    return written


# ---------------------------------------------------------------- collect

def collect_results(jobs_dir: str, runs_dir: str) -> int:
    """
    Copy downloaded job artifacts into the runs/<arm_id>/ layout that
    score_run.py expects. Looks for structure files under each arm's job
    directory (after `apheris-foundry jobs download`), any nesting depth.
    """
    n = 0
    for arm_dir in sorted(glob.glob(os.path.join(jobs_dir, "*"))):
        if not os.path.isdir(arm_dir):
            continue
        arm_id = os.path.basename(arm_dir)
        dest = os.path.join(runs_dir, arm_id)
        os.makedirs(dest, exist_ok=True)
        for f in glob.glob(os.path.join(arm_dir, "**", "*"), recursive=True):
            if f.lower().endswith((".cif", ".mmcif", ".pdb")):
                shutil.copy(f, os.path.join(dest, os.path.basename(f)))
                n += 1
    print(f"collected {n} structure files into {runs_dir}/")
    return n


# ------------------------------------------------------------------- CLI

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("probe", help="print what the Hub workflow supports")
    p.add_argument("--workflow", default="predict")
    p.add_argument("--cli", default="apheris-foundry")

    e = sub.add_parser("emit", help="write Foundry jobs from a kit manifest")
    e.add_argument("--runs", required=True,
                   help="kit output dir containing manifest.json "
                        "(e.g. runs/abl1_imatinib)")
    e.add_argument("--out", default="foundry_jobs")
    e.add_argument("--workflow", default="predict")
    e.add_argument("--schema", default=None,
                   help="probed schema JSON file; if omitted, probe live")

    c = sub.add_parser("collect", help="gather downloaded artifacts for scoring")
    c.add_argument("--jobs", required=True)
    c.add_argument("--runs", required=True)

    args = ap.parse_args()

    if args.cmd == "probe":
        schema = probe_hub_schema(args.workflow, cli=args.cli)
        print(json.dumps(schema, indent=2)[:4000])
        if "_error" not in schema:
            print("\nlever support:")
            for lever in LEVER_CANDIDATES:
                hit = match_levers(schema).get(lever)
                print(f"  {lever:<10s} {'-> ' + hit if hit else 'NOT exposed'}")

    elif args.cmd == "emit":
        manifest = os.path.join(args.runs, "manifest.json")
        if args.schema:
            with open(args.schema) as fh:
                schema = json.load(fh)
        else:
            schema = probe_hub_schema(args.workflow)
            if "_error" in schema:
                print(f"! probe failed ({schema['_error']}); writing jobs with "
                      f"best-guess parameter names -- the Hub's schema "
                      f"validation will reject what's wrong")
                schema = None
        emit_foundry_jobs(manifest, args.out, schema=schema,
                          workflow=args.workflow)
        print(f"submit with: bash {os.path.join(args.out, 'submit_all.sh')}")

    elif args.cmd == "collect":
        collect_results(args.jobs, args.runs)
        print(f"now score with: python score_run.py score --runs {args.runs} ...")


if __name__ == "__main__":
    main()
