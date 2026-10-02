#!/usr/bin/env python3
"""Re-score existing H1/H2/H3 predictions against the ORGANIZERS' OFFICIAL
references (5VB7 agonist, 6T4I antagonist) and H12 window (484-507), using
the organizers' own lib/h12_state.py + apheris-data CLI -- not this
project's state_recovery2/3KYT/4ZJW/479-486 choice.

Requires the rorgt_candidate_kit on disk (unzipped from
~/data/rorgt_candidate_kit_2026-09-30.zip) and an `apheris-data` binary on
PATH (a thin docker wrapper works -- see ~/bin/apheris-data).

    PYTHONPATH=/tmp/rorgt_kit/rorgt_candidate_kit PATH=$HOME/bin:$PATH \
        python3 official_rescore.py --kit-root /tmp/rorgt_kit/rorgt_candidate_kit \
        --out official_state_calls.csv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def collect_structures(roots: dict[str, str]) -> dict[str, Path]:
    structures: dict[str, Path] = {}
    for label, pattern in roots.items():
        base = Path(pattern)
        for cif in sorted(base.rglob("*_model.cif")):
            rel = cif.relative_to(base)
            sid = f"{label}__{str(rel).replace('/', '_').replace('.cif', '')}"
            structures[sid] = cif
    return structures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kit-root", required=True, help="path to unzipped rorgt_candidate_kit")
    ap.add_argument("--out", default="official_state_calls.csv")
    ap.add_argument("--work", default=str(Path.home() / "h3" / "official_rescore_work"))
    args = ap.parse_args()

    sys.path.insert(0, args.kit_root)
    from lib import h12_state  # noqa: E402

    home = Path.home()
    roots = {
        "H1_25hc": str(home / "data/apheris-foundry-steering-committee/apheris_kit_rorgamma/runs/rorgamma_25hc"),
        "H1_invago": str(home / "data/apheris-foundry-steering-committee/apheris_kit_rorgamma/runs/rorgamma_invago"),
        "H2": str(home / "h2/out"),
        "H2_posthoc": str(home / "h2_posthoc/out"),
        "H3": str(home / "h3/out"),
    }
    # H3/out also contains pool120 (the ChEMBL subset) -- score it too, separately labeled.
    structures = collect_structures(roots)
    print(f"collected {len(structures)} structures to score", file=sys.stderr)

    work = Path(args.work)
    work.mkdir(parents=True, exist_ok=True)
    calls = h12_state.call_states(
        structures,
        work,
        default_offset=h12_state.OFFSET_TO_UNIPROT,
    )

    out_csv = Path(args.out)
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "file", "rmsd_to_agonist", "rmsd_to_antagonist", "margin", "state", "n_atoms"])
        for sid, path in structures.items():
            c = calls.get(sid)
            if c is None:
                w.writerow([sid, str(path), "", "", "", "missing", ""])
                continue
            w.writerow([sid, str(path), c.rmsd_to_agonist, c.rmsd_to_antagonist, c.margin, c.state, c.n_atoms])
    print(f"wrote {out_csv}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
