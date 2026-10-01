"""Pocket-set construction for state-specific conditioning.

Four sets are produced from the active/inactive reference pair:

* ``active``   - residues within `cutoff` of the ligand in the active reference
* ``inactive`` - the same for the inactive reference
* ``shared``   - the intersection (the A8 control: pocket constraints that carry
                 no state information)
* ``active_specific`` / ``inactive_specific`` - the set differences, i.e. the
  residues that actually distinguish the two states

All sets are returned in ACTIVE-reference numbering as well as in each
structure's own numbering, because the Foundry request has to name residues in
the numbering of the construct being folded, not of the crystal entry.
"""

from __future__ import annotations

import dataclasses
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import gemmi
import numpy as np

from .structio import ChainModel, align_pair, het_residues

DEFAULT_CUTOFF = 4.0

# Het codes that are almost always crystallisation additives rather than the
# ligand of interest.  Extend at the command line if an entry needs it.
COMMON_ADDITIVES = {
    "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "2PE", "P6G", "DMS", "MES",
    "TRS", "EPE", "ACT", "ACY", "FMT", "SO4", "PO4", "CL", "NA", "K", "MG",
    "CA", "ZN", "NI", "MN", "IOD", "BR", "CIT", "MPD", "BME", "IMD", "NH4",
}


@dataclasses.dataclass
class LigandPick:
    chain_id: str
    name: str
    n_atoms: int
    alternatives: List[Tuple[str, int]]

    def as_dict(self) -> dict:
        return {
            "chain_id": self.chain_id,
            "ccd_code": self.name,
            "n_atoms": self.n_atoms,
            "alternatives": [{"ccd_code": n, "n_atoms": k} for n, k in self.alternatives],
        }


def pick_ligand(
    chain: ChainModel,
    prefer: Optional[str] = None,
    extra_exclude: Iterable[str] = (),
) -> LigandPick:
    """Choose the bound ligand in a structure.

    `prefer` short-circuits the heuristic with an explicit CCD code, which is
    what you want once the code has been read off the entry.  Otherwise the
    largest non-additive het group wins, and every runner-up is reported so the
    choice can be audited rather than trusted.
    """
    cands = het_residues(chain.structure, exclude=set(COMMON_ADDITIVES) | set(extra_exclude))
    if not cands:
        raise ValueError(f"{chain.source}: no candidate ligand het group found")

    scored = sorted(cands, key=lambda cr: len(cr[1]), reverse=True)
    if prefer:
        for cid, res in scored:
            if res.name.upper() == prefer.upper():
                others = [(r.name, len(r)) for c, r in scored if r is not res]
                return LigandPick(cid, res.name, len(res), others)
        raise ValueError(
            f"{chain.source}: requested ligand {prefer} not found; "
            f"available: {[r.name for _, r in scored]}"
        )

    cid, res = scored[0]
    others = [(r.name, len(r)) for c, r in scored[1:]]
    return LigandPick(cid, res.name, len(res), others)


def contact_residues(
    chain: ChainModel,
    ligand: LigandPick,
    cutoff: float = DEFAULT_CUTOFF,
) -> List[int]:
    """Author seq ids of chain residues with any heavy atom within `cutoff`."""
    st = chain.structure
    lig_res = None
    for ch in st[0]:
        if ch.name != ligand.chain_id:
            continue
        for res in ch:
            if res.name == ligand.name and len(res) == ligand.n_atoms:
                lig_res = res
                break
    if lig_res is None:
        raise ValueError(f"{chain.source}: ligand {ligand.name} not re-found")

    lig_xyz = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in lig_res if a.element.name != "H"])
    hits: List[int] = []
    for res in chain.residues:
        xyz = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in res if a.element.name != "H"])
        if len(xyz) == 0:
            continue
        d = np.linalg.norm(xyz[:, None, :] - lig_xyz[None, :, :], axis=-1)
        if d.min() <= cutoff:
            hits.append(res.seqid.num)
    return sorted(hits)


def to_canonical(
    canonical: ChainModel,
    other: ChainModel,
    numbers_in_other: Sequence[int],
) -> List[int]:
    """Translate residue numbers from `other`'s numbering into `canonical`'s."""
    pairs = align_pair(canonical, other)
    rev = {other.numbers[j]: canonical.numbers[i] for i, j in pairs}
    return sorted({rev[n] for n in numbers_in_other if n in rev})


def from_canonical(
    canonical: ChainModel,
    other: ChainModel,
    numbers_in_canonical: Sequence[int],
) -> List[int]:
    pairs = align_pair(canonical, other)
    fwd = {canonical.numbers[i]: other.numbers[j] for i, j in pairs}
    return sorted({fwd[n] for n in numbers_in_canonical if n in fwd})


def build_pocket_sets(
    active: ChainModel,
    inactive: ChainModel,
    active_ligand: LigandPick,
    inactive_ligand: LigandPick,
    cutoff: float = DEFAULT_CUTOFF,
) -> Dict[str, object]:
    """Compute all five pocket sets plus provenance, in canonical numbering."""
    act_own = contact_residues(active, active_ligand, cutoff)
    ina_own = contact_residues(inactive, inactive_ligand, cutoff)

    act = set(act_own)                                   # already canonical
    ina = set(to_canonical(active, inactive, ina_own))

    shared = sorted(act & ina)
    act_spec = sorted(act - ina)
    ina_spec = sorted(ina - act)

    warnings: List[str] = []
    for label, s in (("active", act), ("inactive", ina)):
        if len(s) < 5:
            warnings.append(
                f"{label} pocket set has only {len(s)} residues at {cutoff} A - "
                f"check the ligand pick for that reference before conditioning on it"
            )
    if not act_spec and not ina_spec:
        warnings.append(
            "active and inactive pocket sets are identical: there is no "
            "state-specific conditioning to test, only the A8 shared-pocket control"
        )
    if len(shared) == 0:
        warnings.append(
            "active and inactive pocket sets are disjoint - suspicious for the "
            "same binding site; verify both ligand picks"
        )

    return {
        "cutoff_A": cutoff,
        "canonical_numbering": f"{active.source}:{active.chain_id}",
        "active": sorted(act),
        "inactive": sorted(ina),
        "shared": shared,
        "active_specific": act_spec,
        "inactive_specific": ina_spec,
        "provenance": {
            "active": {
                "source": active.source,
                "chain": active.chain_id,
                "ligand": active_ligand.as_dict(),
                "numbers_in_own_frame": act_own,
            },
            "inactive": {
                "source": inactive.source,
                "chain": inactive.chain_id,
                "ligand": inactive_ligand.as_dict(),
                "numbers_in_own_frame": ina_own,
            },
        },
        "sizes": {
            "active": len(act),
            "inactive": len(ina),
            "shared": len(shared),
            "active_specific": len(act_spec),
            "inactive_specific": len(ina_spec),
        },
        "warnings": warnings,
    }


def scrambled_set(
    chain: ChainModel,
    size: int,
    exclude: Iterable[int] = (),
    seed: int = 0,
) -> List[int]:
    """A size-matched decoy pocket set drawn from non-pocket residues.

    This is the honest negative control for pocket conditioning: same number of
    constrained residues, no pocket meaning.  It replaces a cross-family
    (kinase) reference, whose numbering cannot be mapped onto an LBD at all.
    """
    rng = np.random.default_rng(seed)
    pool = [n for n in chain.numbers if n not in set(exclude)]
    if len(pool) < size:
        raise ValueError("not enough residues outside the excluded set")
    return sorted(int(x) for x in rng.choice(pool, size=size, replace=False))
