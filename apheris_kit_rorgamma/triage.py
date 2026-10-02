"""
Triage the pharma data drop: what did we get, and which targets are worth
the afternoon?

The single biggest scientific risk in this challenge is picking a target
where the baseline already succeeds (famous public complexes like
ABL1-imatinib are almost certainly memorised). The whole point of the
pharma data is that the models have never seen it. This tool finds the
targets where that advantage is real:

    1. ingest whatever landed (PDB/mmCIF structures, a metadata CSV,
       FASTA+SMILES if that's all there is)
    2. group structures by target (sequence similarity)
    3. within each target, cluster by conformation (core-superposed CA RMSD)
       to find two-state pairs
    4. classify kinase DFG state where a DFG motif exists (bonus, not required)
    5. rank targets by "baseline likely to fail here" heuristics, with the
       reasons printed so you can defend the choice to a judge

Usage:
    python triage.py DATA_DIR --out triage_report.csv
    python triage.py DATA_DIR --metadata meta.csv --write-targets new_targets.json --top 3

Nothing here is ground truth -- it is a fast, transparent triage. Confirm the
top candidates by eye (open the structures) before committing the grid.
"""
from __future__ import annotations

import argparse
import csv
import difflib
import glob
import json
import os
import re
from collections import defaultdict

import numpy as np

from structure_io import load_structure
from kinase_state import classify_dfg, find_kinase_motifs
from state_recovery import apply_transform, kabsch, residue_correspondence

STRUCTURE_EXTS = (".pdb", ".cif", ".mmcif")


# ------------------------------------------------------------- ingestion

def parse_deposition_date(path: str) -> str | None:
    """
    Best-effort deposition/release date. PDB HEADER records carry a date in
    cols 51-59 (DD-MMM-YY); mmCIF carries _pdbx_database_status fields.
    Returns ISO-ish string or None. Not all internal pharma files will have
    one -- that is fine, None just means the memorisation check is skipped.
    """
    try:
        with open(path, "r", errors="ignore") as fh:
            head = fh.read(20000)
    except OSError:
        return None

    m = re.search(r"_pdbx_database_status\.recvd_initial_deposition_date\s+"
                  r"['\"]?([0-9]{4}-[0-9]{2}-[0-9]{2})", head)
    if m:
        return m.group(1)

    for line in head.splitlines():
        if line.startswith("HEADER"):
            raw = line[50:59].strip()
            m = re.match(r"(\d{2})-([A-Z]{3})-(\d{2,4})", raw)
            if m:
                months = {"JAN": 1, "FEB": 2, "MAR": 3, "APR": 4, "MAY": 5,
                          "JUN": 6, "JUL": 7, "AUG": 8, "SEP": 9, "OCT": 10,
                          "NOV": 11, "DEC": 12}
                day, mon, yr = m.groups()
                yr = int(yr)
                yr = yr + 2000 if yr < 100 else yr
                return f"{yr:04d}-{months[mon]:02d}-{int(day):02d}"
    return None


def triage_directory(path: str, metadata_csv: str | None = None) -> list:
    """
    Walk a directory and characterise every structure in it. Returns one row
    per file. Sequence-only entries (FASTA) are recorded with seq_only=True
    so you know the pharma side gave you queries, not references.
    """
    meta = {}
    if metadata_csv and os.path.exists(metadata_csv):
        with open(metadata_csv) as fh:
            for row in csv.DictReader(fh):
                key = row.get("id") or row.get("file") or row.get("pdb") or ""
                if key:
                    meta[os.path.splitext(os.path.basename(key))[0]] = row

    rows = []
    files = sorted(
        f for f in glob.glob(os.path.join(path, "**", "*"), recursive=True)
        if f.lower().endswith(STRUCTURE_EXTS + (".fa", ".fasta"))
    )
    for f in files:
        row = {"file": f, "id": os.path.splitext(os.path.basename(f))[0],
               "seq_only": False, "warnings": ""}
        try:
            if f.lower().endswith((".fa", ".fasta")):
                with open(f) as fh:
                    row["sequence"] = "".join(
                        l.strip() for l in fh if not l.startswith(">"))
                row["seq_only"] = True
                rows.append(row)
                continue

            st = load_structure(f)
            row["n_chains"] = len(st.chain_ids())
            row["n_residues"] = len(st.residues)
            row["sequence"] = st.sequence()
            ligs = st.ligands()
            row["ligands"] = ";".join(l.name for l in ligs) or ""
            row["max_ligand_atoms"] = max((len(l.heavy_atoms()) for l in ligs),
                                          default=0)
            row["deposition_date"] = parse_deposition_date(f)

            motifs = find_kinase_motifs(st)
            row["has_dfg"] = motifs.dfg_phe is not None
            if motifs.complete():
                res = classify_dfg(st, motifs)
                row["dfg_label"] = res.label
                row["D1"] = round(res.d1, 2) if res.d1 is not None else ""
                row["D2"] = round(res.d2, 2) if res.d2 is not None else ""
            else:
                row["dfg_label"] = ""
                row["D1"] = row["D2"] = ""
        except Exception as exc:  # a corrupt file must not kill the triage
            row["warnings"] = f"unreadable: {exc}"
        rows.append(row)

    # merge metadata columns we don't already have (deposition date etc.)
    for row in rows:
        m = meta.get(row["id"])
        if m:
            if not row.get("deposition_date"):
                row["deposition_date"] = (m.get("deposition_date")
                                          or m.get("date") or "")
            for k, v in m.items():
                row.setdefault(f"meta_{k}", v)
    return rows


# --------------------------------------------------------- target grouping

def group_by_sequence(rows: list, threshold: float = 0.85) -> dict:
    """
    Assign each row to a target group by sequence identity. Exact-match fast
    path, then difflib ratio for near-identical constructs (different
    boundaries, tags, point mutants). O(n^2) -- fine for a hackathon-sized
    drop, and worth printing a warning if n gets large.
    """
    structs = [r for r in rows if r.get("sequence") and not r.get("seq_only")]
    if len(structs) > 500:
        print(f"  ! {len(structs)} structures: pairwise grouping may be slow")

    groups: dict[int, list] = {}
    rep_seq: dict[int, str] = {}
    for r in structs:
        placed = False
        for gid, rep in rep_seq.items():
            ratio = (1.0 if r["sequence"] == rep else
                     difflib.SequenceMatcher(None, r["sequence"], rep,
                                             autojunk=False).quick_ratio())
            if ratio >= threshold:
                groups[gid].append(r)
                placed = True
                break
        if not placed:
            gid = len(groups)
            groups[gid] = [r]
            rep_seq[gid] = r["sequence"]
    for gid, members in groups.items():
        for r in members:
            r["target_group"] = gid
    return groups


# --------------------------------------------------- conformer clustering

def _ca_rmsd(path_a: str, path_b: str) -> float | None:
    """Core-superposed CA RMSD over aligned residues. None if unalignable."""
    try:
        a, b = load_structure(path_a), load_structure(path_b)
    except Exception:
        return None
    pairs = residue_correspondence(a, b)
    if len(pairs) < 20:
        return None
    pa = np.array([ra.atom("CA") for ra, _ in pairs])
    pb = np.array([rb.atom("CA") for _, rb in pairs])
    R, t = kabsch(pb, pa)
    fit = apply_transform(pb, R, t)
    return float(np.sqrt(((fit - pa) ** 2).sum(axis=1).mean()))


def cluster_conformers(groups: dict, rmsd_cut: float = 2.5) -> dict:
    """
    Within each target group, single-linkage cluster the structures by
    pairwise CA RMSD. Two clusters = two conformational states = a target
    worth steering toward. Returns {gid: {cluster_id: [rows]}}.
    """
    out = {}
    for gid, members in groups.items():
        if len(members) < 2:
            out[gid] = {0: members}
            continue
        n = len(members)
        parent = list(range(n))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for i in range(n):
            for j in range(i + 1, n):
                d = _ca_rmsd(members[i]["file"], members[j]["file"])
                if d is not None and d < rmsd_cut:
                    pi, pj = find(i), find(j)
                    if pi != pj:
                        parent[pi] = pj

        clusters: dict[int, list] = defaultdict(list)
        for i in range(n):
            clusters[find(i)].append(members[i])
        out[gid] = {cid: mem for cid, mem in enumerate(clusters.values())}
    return clusters and out


# --------------------------------------------------------------- ranking

def rank_targets(groups: dict, clusters: dict,
                 cutoff_date: str = "2023-01-01") -> list:
    """
    Score each target group by how likely the baseline is to FAIL on it --
    which is what you want, because a steering effect needs room to exist.

    Heuristics, each printed as a reason so the choice is defensible:
      +3  two or more conformational clusters (there IS a second state)
      +2  every state cluster has a liganded member (pocket derivable)
      +2  any member deposited after cutoff_date (memorisation-free)
      +1  a kinase with a DFG-out / DFG-inter member (rare state)
      +1  largest ligand >= 25 heavy atoms (type-II-like size)

    Set cutoff_date to the actual training cutoff of the model you run on
    the day -- ask the Apheris/OpenFold crew, it is a fair question and the
    answer matters more than any other number here.
    """
    ranked = []
    for gid, members in groups.items():
        cl = clusters.get(gid, {})
        score, reasons = 0, []

        if len(cl) >= 2:
            score += 3
            reasons.append(f"{len(cl)} conformational states present")

        liganded_clusters = sum(
            1 for mem in cl.values()
            if any(m.get("max_ligand_atoms", 0) >= 6 for m in mem))
        if len(cl) >= 2 and liganded_clusters == len(cl):
            score += 2
            reasons.append("every state has a liganded reference")

        dates = [m.get("deposition_date") for m in members
                 if m.get("deposition_date")]
        if any(d >= cutoff_date for d in dates):
            score += 2
            reasons.append(f"post-cutoff structure (>= {cutoff_date})")

        labels = {m.get("dfg_label") for m in members if m.get("dfg_label")}
        if labels & {"DFG-out", "DFG-inter"}:
            score += 1
            reasons.append(f"rare kinase state present: {sorted(labels)}")

        max_atoms = max((m.get("max_ligand_atoms", 0) for m in members),
                        default=0)
        if max_atoms >= 25:
            score += 1
            reasons.append(f"large ligand ({max_atoms} heavy atoms)")

        rep = members[0]
        ranked.append({
            "target_group": gid,
            "n_structures": len(members),
            "n_states": len(cl),
            "representative": os.path.basename(rep["file"]),
            "is_kinase": any(m.get("has_dfg") for m in members),
            "dfg_states": ";".join(sorted(labels)) if labels else "",
            "score": score,
            "reasons": "; ".join(reasons) or "no distinguishing features",
        })
    return sorted(ranked, key=lambda r: -r["score"])


# ------------------------------------------------------- targets.json out

def emit_targets_json(groups: dict, clusters: dict, target_group: int,
                      out_path: str | None = None) -> dict:
    """
    Write a kit-compatible targets.json entry for one target group.

    The two largest conformational clusters become the two states. For
    kinases the rare DFG state is made the target state; otherwise the
    SMALLER cluster is made the target (the under-represented state is the
    one the challenge asks about). CONFIRM this choice by eye before running
    the grid -- the tool cannot know which state your ligand wants.
    """
    cl = clusters[target_group]
    members = groups[target_group]

    ordered = sorted(cl.values(), key=len, reverse=True)
    labels = {}
    for mem in cl.values():
        labs = {m.get("dfg_label") for m in mem if m.get("dfg_label")}
        labels[id(mem)] = sorted(labs)[0] if labs else None

    def pick_ref(mem):
        lig = [m for m in mem if m.get("max_ligand_atoms", 0) >= 6]
        return (lig or mem)[0]

    rare = {"DFG-out", "DFG-inter"}
    if any(labels[id(m)] in rare for m in ordered):
        tgt_mem = next(m for m in ordered if labels[id(m)] in rare)
    else:
        tgt_mem = ordered[-1]  # smallest cluster = under-represented state
    oth_mem = next(m for m in ordered if m is not tgt_mem)

    tgt_ref, oth_ref = pick_ref(tgt_mem), pick_ref(oth_mem)
    tgt_state = labels[id(tgt_mem)] or "state_target"
    oth_state = labels[id(oth_mem)] or "state_other"

    key = f"pharma_g{target_group}"
    entry = {
        key: {
            "name": f"pharma target group {target_group} "
                    f"({len(members)} structures)",
            "uniprot": None,
            "target_state": tgt_state,
            "other_state": oth_state,
            "ref_target_state": tgt_ref["file"],
            "ref_other_state": oth_ref["file"],
            "template_target_state": tgt_ref["file"],
            "template_other_state": oth_ref["file"],
            "sequence": "AUTO",
            "smiles": None,
            "_smiles_note": "FILL THIS IN from the pharma metadata or the "
                            "ligand code's RCSB page before running the grid.",
            "msa_path": None,
            "state_region_seqids": None,
            "_state_region_note": "null = auto-derive from the two references "
                                  "via state_recovery2 (works for non-kinases)",
            "back_pocket_seqids": None,
        }
    }
    if out_path:
        with open(out_path, "w") as fh:
            json.dump(entry, fh, indent=2)
    return entry


# ------------------------------------------------------------------- CLI

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir")
    ap.add_argument("--metadata", default=None, help="optional metadata CSV")
    ap.add_argument("--out", default="triage_report.csv")
    ap.add_argument("--cutoff-date", default="2023-01-01",
                    help="model training cutoff, ISO date (ask on the day)")
    ap.add_argument("--rmsd-cut", type=float, default=2.5,
                    help="CA RMSD cutoff for same-state clustering")
    ap.add_argument("--write-targets", default=None,
                    help="write kit-compatible targets.json for the top groups")
    ap.add_argument("--top", type=int, default=3)
    args = ap.parse_args()

    print(f"triaging {args.data_dir} ...")
    rows = triage_directory(args.data_dir, args.metadata)
    n_seq = sum(1 for r in rows if r.get("seq_only"))
    print(f"  {len(rows)} files ({n_seq} sequence-only)")

    groups = group_by_sequence(rows)
    print(f"  {len(groups)} target groups by sequence")

    clusters = cluster_conformers(groups, rmsd_cut=args.rmsd_cut)
    ranked = rank_targets(groups, clusters, cutoff_date=args.cutoff_date)

    fieldnames = sorted({k for r in rows for k in r})
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote per-structure report to {args.out}\n")

    print("=" * 88)
    print("TARGET RANKING  (higher = baseline more likely to fail = better demo)")
    print("=" * 88)
    for r in ranked:
        print(f"  score {r['score']:>2d}  group {r['target_group']:<3d} "
              f"({r['n_structures']} structures, {r['n_states']} states, "
              f"rep {r['representative']})")
        print(f"         {r['reasons']}")
    print("-" * 88)
    print("Confirm the top candidates BY EYE before committing the grid.")

    if args.write_targets and ranked:
        entries = {}
        for r in ranked[:args.top]:
            if r["n_states"] < 2:
                print(f"  ! group {r['target_group']} has one state only, "
                      f"skipping for targets.json")
                continue
            entries.update(emit_targets_json(groups, clusters,
                                             r["target_group"]))
        with open(args.write_targets, "w") as fh:
            json.dump(entries, fh, indent=2)
        print(f"\nwrote {len(entries)} candidate target configs to "
              f"{args.write_targets}")
        print("NEXT: fill in each 'smiles' field, then run "
              "score_run.py validate on the references.")


if __name__ == "__main__":
    main()
