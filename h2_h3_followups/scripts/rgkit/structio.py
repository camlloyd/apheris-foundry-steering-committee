"""Structure I/O, residue correspondence and superposition helpers.

Design notes
------------
* Residue correspondence between two RORgamma LBD structures is established by
  SEQUENCE ALIGNMENT, never by assuming shared author numbering.  RORgamma and
  RORgamma-t are different isoforms whose LBD residues differ by a constant
  offset, and crystal constructs differ further.  Anything that hard-codes a
  residue number (e.g. "H479") is therefore treated as a label to be *verified*
  against the structure, not as an index to index with.
* Superposition is done on a user-supplied subset of positions (normally the
  complement of the state region) so that a hinge motion in one helix is not
  smeared across the whole domain by a global fit.
"""

from __future__ import annotations

import dataclasses
import gzip
import pathlib
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import gemmi
import numpy as np

# Residues we never treat as part of the protein chain.
_NON_POLYMER_SKIP = {"HOH", "DOD", "WAT"}


@dataclasses.dataclass
class ChainModel:
    """CA trace of one polymer chain, plus the het groups in the same model."""

    source: str                      # file path or PDB id
    chain_id: str
    seq: str                         # one-letter sequence, CA-complete residues only
    numbers: List[int]               # author seq ids, parallel to seq
    ca: np.ndarray                   # (N, 3) CA coordinates, parallel to seq
    residues: List[gemmi.Residue]    # parallel to seq, for all-atom work
    structure: gemmi.Structure       # the parent structure (model 0 retained)

    def __len__(self) -> int:
        return len(self.seq)

    def index_of_number(self, num: int) -> Optional[int]:
        try:
            return self.numbers.index(num)
        except ValueError:
            return None


def read_structure(path: str | pathlib.Path) -> gemmi.Structure:
    """Read a .cif/.pdb file, optionally gzipped, and keep only the first model."""
    p = pathlib.Path(path)
    if p.suffix == ".gz":
        text = gzip.decompress(p.read_bytes()).decode("utf-8", "replace")
        if ".cif" in p.name:
            doc = gemmi.cif.read_string(text)
            st = gemmi.make_structure_from_block(doc.sole_block())
        else:
            st = gemmi.read_pdb_string(text)
    else:
        st = gemmi.read_structure(str(p))
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    while len(st) > 1:
        del st[1]
    return st


def _one_letter(res: gemmi.Residue) -> str:
    info = gemmi.find_tabulated_residue(res.name)
    if info is None or not info.is_amino_acid():
        return ""
    code = gemmi.find_tabulated_residue(res.name).one_letter_code
    return code.upper() if code else "X"


def polymer_chains(st: gemmi.Structure, min_len: int = 20) -> List[str]:
    """Chain ids with at least `min_len` CA-bearing amino acids."""
    out = []
    for ch in st[0]:
        n = sum(1 for r in ch if _one_letter(r) and r.find_atom("CA", "*"))
        if n >= min_len:
            out.append(ch.name)
    return out


def load_chain(
    path: str | pathlib.Path,
    chain_id: Optional[str] = None,
    min_len: int = 20,
    source: Optional[str] = None,
) -> ChainModel:
    """Load the CA trace of one chain.

    If `chain_id` is None the longest polymer chain is used, which for an LBD
    entry is the receptor rather than a bound coactivator peptide.
    """
    st = read_structure(path)
    candidates = polymer_chains(st, min_len=min_len)
    if not candidates:
        raise ValueError(f"{path}: no polymer chain with >= {min_len} residues")

    if chain_id is None:
        best, best_n = None, -1
        for cid in candidates:
            ch = st[0][cid]
            n = sum(1 for r in ch if _one_letter(r) and r.find_atom("CA", "*"))
            if n > best_n:
                best, best_n = cid, n
        chain_id = best
    if chain_id not in [c.name for c in st[0]]:
        raise ValueError(f"{path}: chain {chain_id} not present")

    seq, numbers, coords, residues = [], [], [], []
    for res in st[0][chain_id]:
        if res.name in _NON_POLYMER_SKIP:
            continue
        letter = _one_letter(res)
        if not letter:
            continue
        ca = res.find_atom("CA", "*")
        if ca is None:
            continue
        seq.append(letter)
        numbers.append(res.seqid.num)
        coords.append([ca.pos.x, ca.pos.y, ca.pos.z])
        residues.append(res)

    return ChainModel(
        source=str(source or path),
        chain_id=chain_id,
        seq="".join(seq),
        numbers=numbers,
        ca=np.asarray(coords, dtype=float),
        residues=residues,
        structure=st,
    )


def het_residues(
    st: gemmi.Structure,
    exclude: Iterable[str] = (),
    min_atoms: int = 6,
) -> List[Tuple[str, gemmi.Residue]]:
    """Non-polymer residues plausibly representing a bound ligand.

    Returns (chain_id, residue) pairs.  Waters, ions and cryoprotectants that
    fall below `min_atoms` are dropped; the caller still has to decide which of
    the survivors is the ligand of interest.
    """
    skip = set(_NON_POLYMER_SKIP) | {s.upper() for s in exclude}
    out = []
    for ch in st[0]:
        for res in ch:
            if res.name in skip:
                continue
            if _one_letter(res):            # amino acid -> part of polymer
                continue
            info = gemmi.find_tabulated_residue(res.name)
            if info is not None and info.is_nucleic_acid():
                continue
            if len(res) < min_atoms:
                continue
            out.append((ch.name, res))
    return out


def _parse_cigar(cigar: str) -> List[Tuple[int, str]]:
    """Parse a gemmi cigar string such as '5M3I5M1D' into (length, op) pairs."""
    out, num = [], ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            out.append((int(num) if num else 1, ch))
            num = ""
    return out


def align_pair(a: ChainModel, b: ChainModel) -> List[Tuple[int, int]]:
    """Index pairs (i_in_a, j_in_b) for aligned, non-gap positions.

    Uses gemmi's global sequence alignment so that differing construct
    boundaries, isoform offsets and disordered loops are handled explicitly.
    """
    result = gemmi.align_string_sequences(list(a.seq), list(b.seq), [])
    pairs: List[Tuple[int, int]] = []
    i = j = 0
    for length, op in _parse_cigar(result.cigar_str()):
        if op == "M":
            for k in range(length):
                if a.seq[i + k] == b.seq[j + k]:
                    pairs.append((i + k, j + k))
            i += length
            j += length
        elif op in ("I", "S"):
            i += length
        elif op in ("D", "N"):
            j += length
        else:  # pragma: no cover - defensive
            raise ValueError(f"unexpected cigar op {op!r}")
    return pairs


def kabsch(mob: np.ndarray, ref: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Rotation R and translation t minimising |R*mob + t - ref|."""
    mc, rc = mob.mean(axis=0), ref.mean(axis=0)
    m, r = mob - mc, ref - rc
    u, _, vt = np.linalg.svd(m.T @ r)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rot, rc - rot @ mc


def superpose_on(
    mob: np.ndarray,
    ref: np.ndarray,
    fit_mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Apply the transform fitted on `fit_mask` rows to all rows of `mob`."""
    if fit_mask is None:
        fit_mask = np.ones(len(mob), dtype=bool)
    if fit_mask.sum() < 3:
        raise ValueError("need at least 3 positions to superpose")
    rot, t = kabsch(mob[fit_mask], ref[fit_mask])
    return (rot @ mob.T).T + t


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))
