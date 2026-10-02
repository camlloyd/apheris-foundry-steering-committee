"""
Kinase conformational state classification (DFG-in / DFG-out / DFG-inter).

Implements the Modi-Dunbrack "KinCore" spatial criteria, which place the
DFG-Phe side chain using two distances:

    D1 = dist( alphaC-Glu(+4) CA , DFG-Phe CZ )
    D2 = dist( beta3-Lys CA     , DFG-Phe CZ )

    DFG-in     : D1 <  11.0 AND D2 >  14.0
    DFG-out    : D1 >  11.0 AND D2 <  14.0
    DFG-inter  : D1 <  11.0 AND D2 <  11.0
    otherwise  : unassigned

  Reference: Modi & Dunbrack, PNAS 2019; Kincore, NAR 2022 (D50:D654).
  https://dunbrack.fccc.edu/kincore/

IMPORTANT, read before you trust a number:
  Different papers quote these thresholds slightly differently, and KinCore
  has revised its criteria at least once. The thresholds live in DFG_THRESHOLDS
  below so you can change them in one place. On Day 1, the FIRST thing to do is
  run this on your reference crystal structures and check that the known
  DFG-out one comes back DFG-out. If it doesn't, the thresholds or the motif
  detection are wrong for your target, not the model you're testing.

Motif detection is automatic but always printed, so you can eyeball it. Any
motif can be overridden by hand in targets.yaml -- do that for anything
important rather than trusting the autodetect.
"""
from __future__ import annotations

import difflib
from dataclasses import dataclass, field

import numpy as np

from structure_io import Structure

# {label: (d1_min, d1_max, d2_min, d2_max)}, None means unbounded.
DFG_THRESHOLDS = {
    "DFG-in":    (None, 11.0, 14.0, None),
    "DFG-out":   (11.0, None, None, 14.0),
    "DFG-inter": (None, 11.0, None, 11.0),
}

# HRD-like motifs: the catalytic loop. Used to disambiguate which DFG is the
# real one when the sequence contains more than one.
HRD_PATTERNS = ["HRD", "HRE", "HGD", "YRD", "HCD", "HKD"]

# DFG and its natural variants. Only DFG itself gives a Phe CZ to measure.
DFG_PATTERNS = ["DFG", "DLG", "DWG", "DYG", "DMG", "SFG"]

# b3-strand anchor: the conserved Lys sits at the end of a hydrophobic
# triplet (VAVK in ABL1/p38a, VAIK in some other kinases). Sequence-anchored,
# so it survives the broken K-E salt bridge of DFG-out / alphaC-out
# structures, where a distance-only search latches onto decoy Lys/Glu pairs.
B3_HYDROPHOBIC = {"VAL", "ILE", "LEU", "ALA", "MET"}


@dataclass
class KinaseMotifs:
    """Residue numbers (PDB/author numbering) of the motifs we measure from."""
    dfg_asp: int | None = None
    dfg_phe: int | None = None
    dfg_gly: int | None = None
    b3_lys: int | None = None
    ac_glu: int | None = None
    ac_glu_p4: int | None = None      # alphaC-Glu + 4, the D1 anchor
    hrd_his: int | None = None
    chain: str = "A"
    warnings: list = field(default_factory=list)

    def complete(self) -> bool:
        return None not in (self.dfg_phe, self.b3_lys, self.ac_glu_p4)

    def describe(self, structure: Structure | None = None) -> str:
        def fmt(label, seqid):
            if seqid is None:
                return f"{label}=?"
            nm = ""
            if structure is not None:
                r = structure.residue_by_seqid(seqid, self.chain)
                if r is not None:
                    nm = r.name
            return f"{label}={nm}{seqid}"

        parts = [
            fmt("HRD-His", self.hrd_his),
            fmt("b3-Lys", self.b3_lys),
            fmt("aC-Glu", self.ac_glu),
            fmt("aC-Glu+4", self.ac_glu_p4),
            fmt("DFG-Asp", self.dfg_asp),
            fmt("DFG-Phe", self.dfg_phe),
        ]
        line = "  ".join(parts)
        if self.warnings:
            line += "\n    ! " + "\n    ! ".join(self.warnings)
        return line


@dataclass
class DFGResult:
    label: str
    d1: float | None
    d2: float | None
    motifs: KinaseMotifs
    source: str = ""

    def __str__(self) -> str:
        d1 = f"{self.d1:.2f}" if self.d1 is not None else "  n/a"
        d2 = f"{self.d2:.2f}" if self.d2 is not None else "  n/a"
        return f"{self.label:<12s} D1={d1} A  D2={d2} A"


def _find_motif_one_gapped(seq: str, pattern: str, max_gap: int = 2) -> list:
    """
    Occurrences of `pattern` allowing up to `max_gap` unknown (X) residues
    inside it. Handles structures where a motif residue is unmodeled -- e.g.
    1KV2 (p38a + BIRB-796) has Gly170 of DFG missing, so the parsed sequence
    reads ...VDIW... and a contiguous search for 'DFG' finds nothing.

    Returns [(start_index, matched_string, n_gaps), ...] in sequence order.
    """
    hits = []
    start = 0
    while True:
        i = seq.find(pattern[0], start)
        if i < 0 or i + len(pattern) > len(seq):
            break
        window = seq[i:i + len(pattern)]
        if window == pattern:
            hits.append((i, pattern, 0))
        elif (window[0] == pattern[0] and window[-1] == pattern[-1]
              and window[1:-1].count("X") <= max_gap
              and all(a == b or b == "X" for a, b in zip(pattern[1:-1],
                                                         window[1:-1]))):
            hits.append((i, window, window.count("X")))
        start = i + 1
    return hits


def _find_motif(seq: str, patterns: list, max_gap: int = 2) -> list:
    """All occurrences of any pattern, gap-tolerant, in sequence order."""
    hits = []
    for pat in patterns:
        hits.extend(_find_motif_one_gapped(seq, pat, max_gap))
    return sorted(hits)


def find_kinase_motifs(structure: Structure, chain: str | None = None) -> KinaseMotifs:
    """
    Locate the motifs needed for DFG classification.

    Strategy, in order of reliability:
      1. DFG   - sequence motif, disambiguated by requiring it to sit 15-50
                 residues after an HRD-like catalytic-loop motif.
      2. b3-Lys- the VAIK lysine, found as the Lys whose NZ is closest to a
                 Glu carboxylate (the conserved K-E salt bridge), searching
                 only N-terminal to the DFG.
      3. aC-Glu- the Glu partner of that salt bridge; aC-Glu+4 is 4 residues on.

    Step 2/3 are structural rather than sequence-based because the alphaC
    helix has no reliable sequence signature. They are done on whatever
    structure you pass -- so run them on a good CRYSTAL structure and transfer
    the result with transfer_motifs(), rather than autodetecting on a
    prediction whose alphaC may be displaced.
    """
    if chain is None:
        chain = structure.chain_ids()[0] if structure.chain_ids() else "A"

    residues = structure.chain(chain)
    if not residues:
        return KinaseMotifs(chain=chain, warnings=[f"chain {chain} has no residues"])

    seq = "".join(r.one_letter for r in residues)
    motifs = KinaseMotifs(chain=chain)

    # --- HRD (catalytic loop), used only to disambiguate DFG ---
    hrd_hits = _find_motif(seq, HRD_PATTERNS)
    hrd_idx = hrd_hits[0][0] if hrd_hits else None
    if hrd_idx is not None:
        motifs.hrd_his = residues[hrd_idx].seqid

    # --- DFG ---
    dfg_hits = _find_motif(seq, DFG_PATTERNS)
    if not dfg_hits:
        # Fallback: anchor on the two-residue Asp-Phe pair. Handles a motif
        # residue that is unmodeled in the deposition -- 1KV2 (p38a + BIRB-796)
        # has Gly170 missing, so the triplet 'DFG' never appears in the parsed
        # sequence. Classification only needs the Phe CZ, so this is safe.
        dfg_hits = [(i, pat[:2], 1)
                    for pat in ("DF", "DL", "DW", "DY", "DM", "SF")
                    for (i, _m, _g) in _find_motif(seq, [pat])]
        if dfg_hits:
            motifs.warnings.append(
                "no contiguous DFG-like triplet in the parsed sequence; "
                "anchored on the Asp-Phe pair (a motif residue, likely the "
                "DFG Gly, is unmodeled). Fine for D1/D2 -- they need the Phe."
            )
    if not dfg_hits:
        motifs.warnings.append("no DFG-like motif found in sequence")
        return motifs

    chosen = None
    if hrd_idx is not None:
        plausible = [h for h in dfg_hits if 15 <= (h[0] - hrd_idx) <= 50]
        if plausible:
            chosen = plausible[0]
    if chosen is None:
        chosen = dfg_hits[0]
        if len(dfg_hits) > 1:
            motifs.warnings.append(
                f"{len(dfg_hits)} DFG-like motifs found "
                f"({', '.join(f'{p}@{residues[i].seqid}' for i, p, _g in dfg_hits)}); "
                f"using the first. Override dfg_phe in targets.yaml if wrong."
            )

    dfg_i, dfg_pat, dfg_gaps = chosen
    motifs.dfg_asp = residues[dfg_i].seqid
    motifs.dfg_phe = residues[dfg_i + 1].seqid if dfg_i + 1 < len(residues) else None
    if dfg_gaps == 0 and dfg_i + 2 < len(residues):
        motifs.dfg_gly = residues[dfg_i + 2].seqid
    if residues[dfg_i + 1].name.upper() != "PHE":
        motifs.warnings.append(
            f"motif middle residue is {residues[dfg_i + 1].name}, not PHE; "
            f"D1/D2 need a Phe CZ atom -- this target needs a different "
            f"order parameter or a hand-set dfg_phe."
        )

    # --- b3-Lys: beta-3 strand anchor. A Lys preceded by three hydrophobic
    # residues (VAVK/VAIK-type). Rank candidates: intact salt bridge first
    # (Glu partner within 8 A), then any Glu partner in the +10..+28 window,
    # then earliest in sequence -- the conserved VAVK is the first of them.
    # Distance-only choice fails on DFG-out references (bridge broken);
    # last-hit choice fails on decoy hydrophobic Lys downstream (K291 in ABL1).
    lys_res = None
    b3_hits = [i for i in range(3, dfg_i)
               if residues[i].name.upper() == "LYS"
               and residues[i].atom("NZ") is not None
               and all(residues[i - k].name.upper() in B3_HYDROPHOBIC
                       for k in (1, 2, 3))]

    def _glu_dist(i):
        nz = residues[i].atom("NZ")
        return min((float(np.linalg.norm(nz - g.atom("OE1")))
                    for g in residues[:dfg_i]
                    if g.name.upper() == "GLU" and g.atom("OE1") is not None
                    and residues[i].seqid + 10 <= g.seqid <= residues[i].seqid + 28),
                   default=None)

    def _rank(i):
        d = _glu_dist(i)
        return (0 if d is not None and d <= 8.0 else
                1 if d is not None else 2, i)

    if b3_hits:
        best_i = min(b3_hits, key=_rank)
        lys_res = residues[best_i]
        motifs.b3_lys = lys_res.seqid
        d = _glu_dist(best_i)
        if len(b3_hits) > 1:
            motifs.warnings.append(
                f"{len(b3_hits)} beta-3 Lys candidates "
                f"({', '.join(str(residues[i].seqid) for i in b3_hits)}); chose "
                f"{lys_res.seqid}. Verify, or set motif_overrides."
            )
        if d is None:
            motifs.warnings.append(
                "b3-Lys has no Glu partner in the +10..+28 window -- verify, "
                "or set motif_overrides"
            )
    else:
        motifs.warnings.append(
            "no beta-3 Lys found before the DFG -- set b3_lys via "
            "motif_overrides in the target config"
        )
        return motifs

    # --- aC-Glu: the Glu of the conserved salt bridge, restricted to the
    # alphaC region C-terminal of the VAIK Lys (+10..+28). The distance test
    # alone fails exactly on DFG-out references, where the bridge is broken.
    glu_candidates = [
        r for r in residues[:dfg_i]
        if r.name.upper() == "GLU" and r.atom("OE1") is not None
        and lys_res.seqid + 10 <= r.seqid <= lys_res.seqid + 28
    ]
    best = None  # (distance, glu)
    for glu in glu_candidates:
        d = float(np.linalg.norm(lys_res.atom("NZ") - glu.atom("OE1")))
        if best is None or d < best[0]:
            best = (d, glu)

    if best is None:
        motifs.warnings.append(
            "no aC-Glu found in the VAIK+10..+28 window -- set ac_glu and "
            "ac_glu_p4 by hand via motif_overrides in the target config"
        )
        return motifs

    dist, glu = best
    motifs.ac_glu = glu.seqid
    if dist > 8.0:
        motifs.warnings.append(
            f"closest in-window Lys-Glu pair is {dist:.1f} A apart -- the salt "
            f"bridge is broken (alphaC-out; expected for DFG-out type II "
            f"complexes). If D1/D2 look wrong, set motif_overrides."
        )

    p4 = structure.residue_by_seqid(glu.seqid + 4, chain)
    if p4 is None:
        idx = next((i for i, r in enumerate(residues) if r.seqid == glu.seqid), None)
        if idx is not None and idx + 4 < len(residues):
            p4 = residues[idx + 4]
            motifs.warnings.append("aC-Glu+4 taken by index (numbering has gaps)")
    motifs.ac_glu_p4 = p4.seqid if p4 is not None else None

    return motifs


def motifs_from_overrides(structure: Structure, overrides: dict) -> KinaseMotifs:
    """
    Build motifs from explicit config values -- targets.json 'motif_overrides'.

    Autodetection fails in two known ways: the salt-bridge distance search
    picks decoy anchors when alphaC is out (1IEP), and unmodeled motif
    residues break the sequence search (1KV2's missing DFG Gly). For anything
    important, set the motifs by hand here; seqids are in the reference
    structure's own numbering. Example:

        "motif_overrides": {"chain": "A", "b3_lys": 271, "ac_glu": 286,
                            "ac_glu_p4": 290, "dfg_asp": 381, "dfg_phe": 382}
    """
    m = KinaseMotifs(chain=overrides.get("chain", "A"))
    for field in ("dfg_asp", "dfg_phe", "dfg_gly", "b3_lys",
                  "ac_glu", "ac_glu_p4", "hrd_his"):
        v = overrides.get(field)
        if v is not None:
            setattr(m, field, int(v))
    return m


def transfer_motifs(
    ref: Structure,
    ref_motifs: KinaseMotifs,
    target: Structure,
    ref_chain: str | None = None,
    target_chain: str | None = None,
) -> KinaseMotifs:
    """
    Move motif assignments from a reference structure onto a prediction.

    Predictions are numbered 1..N over the input construct; crystal structures
    use author numbering with gaps. This aligns the two sequences and remaps
    every motif residue, so you only ever autodetect once, on the reference.
    """
    ref_res = ref.chain(ref_chain)
    tgt_res = target.chain(target_chain)
    ref_seq = "".join(r.one_letter for r in ref_res)
    tgt_seq = "".join(r.one_letter for r in tgt_res)

    matcher = difflib.SequenceMatcher(None, ref_seq, tgt_seq, autojunk=False)
    index_map = {}
    for block in matcher.get_matching_blocks():
        for k in range(block.size):
            index_map[block.a + k] = block.b + k

    ref_pos = {r.seqid: i for i, r in enumerate(ref_res)}
    out = KinaseMotifs(
        chain=tgt_res[0].chain if tgt_res else (target_chain or "A"),
        warnings=list(ref_motifs.warnings),
    )

    coverage = len(index_map) / max(len(ref_seq), 1)
    if coverage < 0.5:
        out.warnings.append(
            f"only {coverage:.0%} of the reference sequence aligned to the "
            f"prediction -- check you are comparing the same construct"
        )

    for field_name in ("dfg_asp", "dfg_phe", "dfg_gly", "b3_lys",
                       "ac_glu", "ac_glu_p4", "hrd_his"):
        seqid = getattr(ref_motifs, field_name)
        if seqid is None:
            continue
        i = ref_pos.get(seqid)
        if i is None or i not in index_map:
            out.warnings.append(f"{field_name} ({seqid}) did not map onto the prediction")
            continue
        setattr(out, field_name, tgt_res[index_map[i]].seqid)

    return out


def classify_dfg(
    structure: Structure,
    motifs: KinaseMotifs | None = None,
    chain: str | None = None,
) -> DFGResult:
    """Compute D1/D2 and return the KinCore spatial label."""
    if motifs is None:
        motifs = find_kinase_motifs(structure, chain)

    if not motifs.complete():
        return DFGResult("unassigned", None, None, motifs, source=structure.path)

    ch = motifs.chain
    phe = structure.residue_by_seqid(motifs.dfg_phe, ch)
    lys = structure.residue_by_seqid(motifs.b3_lys, ch)
    glu4 = structure.residue_by_seqid(motifs.ac_glu_p4, ch)

    missing = [n for n, r in (("DFG-Phe", phe), ("b3-Lys", lys), ("aC-Glu+4", glu4))
               if r is None]
    if missing:
        motifs.warnings.append(f"missing residues in this file: {', '.join(missing)}")
        return DFGResult("unassigned", None, None, motifs, source=structure.path)

    cz = phe.atom("CZ")
    lys_ca = lys.atom("CA")
    glu4_ca = glu4.atom("CA")
    if cz is None:
        motifs.warnings.append(
            "DFG-Phe has no CZ atom (truncated side chain) -- cannot classify"
        )
        return DFGResult("unassigned", None, None, motifs, source=structure.path)
    if lys_ca is None or glu4_ca is None:
        motifs.warnings.append("b3-Lys or aC-Glu+4 has no CA atom")
        return DFGResult("unassigned", None, None, motifs, source=structure.path)

    d1 = float(np.linalg.norm(glu4_ca - cz))
    d2 = float(np.linalg.norm(lys_ca - cz))

    label = "unassigned"
    for name, (d1_lo, d1_hi, d2_lo, d2_hi) in DFG_THRESHOLDS.items():
        if d1_lo is not None and d1 <= d1_lo:
            continue
        if d1_hi is not None and d1 >= d1_hi:
            continue
        if d2_lo is not None and d2 <= d2_lo:
            continue
        if d2_hi is not None and d2 >= d2_hi:
            continue
        label = name
        break

    return DFGResult(label, d1, d2, motifs, source=structure.path)
