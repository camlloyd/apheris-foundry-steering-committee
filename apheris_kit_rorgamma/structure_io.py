"""
Minimal structure I/O for the Apheris multi-state hackathon kit.

Design goal: run anywhere. Uses gemmi when available (handles mmCIF properly,
which is what Boltz-2 / OpenFold3 emit), and falls back to a pure-Python PDB
parser when it isn't installed. No other hard dependencies beyond numpy.

Everything downstream works on a simple Structure object so the metric code
never has to care which parser produced it.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field

import numpy as np

try:
    import gemmi

    HAVE_GEMMI = True
except ImportError:  # pragma: no cover - environment dependent
    HAVE_GEMMI = False


THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # common modified residues seen in kinase structures
    "MSE": "M", "SEP": "S", "TPO": "T", "PTR": "Y", "CSO": "C",
    "HIC": "H", "KCX": "K", "LLP": "K", "MLY": "K",
}

# Things that are in the file but are not the protein and not the ligand.
SOLVENT_AND_IONS = {
    "HOH", "DOD", "WAT", "SO4", "PO4", "GOL", "EDO", "PEG", "PGE", "MPD",
    "DMS", "ACT", "FMT", "CL", "NA", "K", "MG", "CA", "ZN", "MN", "NI",
    "CD", "CO", "FE", "CU", "BR", "IOD", "TRS", "EPE", "IMD", "NO3",
}


@dataclass
class Residue:
    seqid: int
    name: str
    chain: str
    atoms: dict = field(default_factory=dict)  # atom_name -> np.array([x, y, z])
    is_polymer: bool = True

    @property
    def one_letter(self) -> str:
        return THREE_TO_ONE.get(self.name.upper(), "X")

    def atom(self, name: str):
        """Coordinates of one atom, or None if absent (disordered side chains happen)."""
        return self.atoms.get(name)

    def heavy_atoms(self) -> np.ndarray:
        """All non-hydrogen atom coordinates, (n, 3)."""
        coords = [
            xyz for nm, xyz in self.atoms.items()
            if not nm.startswith("H") and not re.match(r"^\d+H", nm)
        ]
        return np.array(coords) if coords else np.zeros((0, 3))


@dataclass
class Structure:
    path: str
    residues: list = field(default_factory=list)   # polymer residues, in order
    hetero: list = field(default_factory=list)     # non-polymer residues (ligands etc.)

    def chain(self, chain_id: str | None = None) -> list:
        """Polymer residues for one chain. Defaults to the first chain present."""
        if chain_id is None:
            if not self.residues:
                return []
            chain_id = self.residues[0].chain
        return [r for r in self.residues if r.chain == chain_id]

    def chain_ids(self) -> list:
        seen = []
        for r in self.residues:
            if r.chain not in seen:
                seen.append(r.chain)
        return seen

    def sequence(self, chain_id: str | None = None) -> str:
        return "".join(r.one_letter for r in self.chain(chain_id))

    def residue_by_seqid(self, seqid: int, chain_id: str | None = None):
        for r in self.chain(chain_id):
            if r.seqid == seqid:
                return r
        return None

    def ca_coords(self, chain_id: str | None = None):
        """Returns (coords, seqids) for residues that actually have a CA."""
        coords, seqids = [], []
        for r in self.chain(chain_id):
            ca = r.atom("CA")
            if ca is not None:
                coords.append(ca)
                seqids.append(r.seqid)
        return np.array(coords) if coords else np.zeros((0, 3)), seqids

    def ligands(self, min_atoms: int = 6) -> list:
        """
        Non-polymer residues big enough to be a real ligand. Filters out water,
        buffer components and ions, which otherwise dominate the hetero list.
        """
        out = []
        for r in self.hetero:
            if r.name.upper() in SOLVENT_AND_IONS:
                continue
            if len(r.heavy_atoms()) < min_atoms:
                continue
            out.append(r)
        return out


def _load_gemmi(path: str) -> Structure:
    st = gemmi.read_structure(path)
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()

    out = Structure(path=path)
    model = st[0]
    for ch in model:
        for res in ch:
            atoms = {a.name: np.array([a.pos.x, a.pos.y, a.pos.z]) for a in res}
            is_poly = res.name.upper() in THREE_TO_ONE
            r = Residue(
                seqid=res.seqid.num,
                name=res.name,
                chain=ch.name,
                atoms=atoms,
                is_polymer=is_poly,
            )
            (out.residues if is_poly else out.hetero).append(r)
    return out


def _load_pdb_fallback(path: str) -> Structure:
    """Pure-Python PDB reader. Only used when gemmi isn't installed."""
    out = Structure(path=path)
    current: dict = {}

    with open(path) as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break  # first model only
            if not line.startswith(("ATOM  ", "HETATM")):
                continue

            altloc = line[16]
            if altloc not in (" ", "A"):
                continue

            name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21].strip() or "A"
            try:
                seqid = int(line[22:26])
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue

            if name.startswith("H") or re.match(r"^\d+H", name):
                continue

            key = (chain, seqid, resname)
            if key not in current:
                is_poly = resname.upper() in THREE_TO_ONE
                r = Residue(seqid=seqid, name=resname, chain=chain,
                            atoms={}, is_polymer=is_poly)
                current[key] = r
                (out.residues if is_poly else out.hetero).append(r)
            current[key].atoms[name] = np.array([x, y, z])

    return out


def load_structure(path: str) -> Structure:
    """Load a .pdb / .cif / .mmcif file into a Structure."""
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    if HAVE_GEMMI:
        return _load_gemmi(path)
    if path.lower().endswith((".cif", ".mmcif")):
        raise RuntimeError(
            f"{path} is mmCIF and gemmi isn't installed. "
            "Run: pip install gemmi   (or convert the file to PDB first)"
        )
    return _load_pdb_fallback(path)


def _pdb_atom_name(name: str) -> str:
    """
    PDB atom names occupy columns 13-16, and single-letter elements are
    indented by one so the element symbol sits in 13-14. Getting this wrong
    makes parsers read 'CA' as calcium instead of an alpha carbon.
    """
    if len(name) >= 4:
        return name[:4]
    if name[0].isdigit():
        return f"{name:<4s}"
    return f" {name:<3s}"


def write_pdb(structure: Structure, path: str) -> None:
    """Minimal PDB writer, used by the self-test to build synthetic cases."""
    serial = 1
    with open(path, "w") as fh:
        for res in structure.residues + structure.hetero:
            record = "ATOM  " if res.is_polymer else "HETATM"
            for aname, xyz in res.atoms.items():
                element = "".join(c for c in aname if c.isalpha())[:1]
                fh.write(
                    # cols: 1-6 record, 7-11 serial, 13-16 name, 17 altLoc,
                    #       18-20 resName, 22 chain, 23-26 resSeq, 31-54 xyz
                    f"{record}{serial:5d} {_pdb_atom_name(aname)} "
                    f"{res.name:>3s} {res.chain}{res.seqid:4d}    "
                    f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
                    f"{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n"
                )
                serial += 1
        fh.write("END\n")
