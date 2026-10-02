"""
Dual-reference state recovery: which conformational basin did the model land in?

This is the metric the hackathon challenge actually asks about, and it is NOT
the same thing as ligand pose accuracy. KinConfBench (npj Drug Discovery, 2026)
found that geometric pose success -- ligand RMSD < 2 A, high lDDT-PLI -- does
not correlate strongly with recovering the correct protein conformational
state. A model can put the ligand in roughly the right place while leaving the
protein in the wrong state entirely.

So everything here reports state recovery and pose accuracy SEPARATELY, and
the headline number is state recovery.

The one trick that matters: superpose on the invariant core, then measure RMSD
on the state-defining region only. Superposing and measuring globally lets the
shared fold (which is ~95% of the atoms and identical between states) swamp the
local rearrangement you are trying to detect.
"""
from __future__ import annotations

import difflib
from dataclasses import dataclass, field

import numpy as np

from structure_io import Structure


# ---------------------------------------------------------------- alignment

def residue_correspondence(
    a: Structure,
    b: Structure,
    a_chain: str | None = None,
    b_chain: str | None = None,
) -> list:
    """
    Pair up residues between two structures of the same protein by sequence.

    Returns [(residue_in_a, residue_in_b), ...] for aligned positions where
    both sides have a CA atom. Handles the common case where one structure is
    a crystal structure with author numbering and gaps, and the other is a
    prediction numbered 1..N.
    """
    ra, rb = a.chain(a_chain), b.chain(b_chain)
    sa = "".join(r.one_letter for r in ra)
    sb = "".join(r.one_letter for r in rb)

    matcher = difflib.SequenceMatcher(None, sa, sb, autojunk=False)
    pairs = []
    for block in matcher.get_matching_blocks():
        for k in range(block.size):
            x, y = ra[block.a + k], rb[block.b + k]
            if x.atom("CA") is not None and y.atom("CA") is not None:
                pairs.append((x, y))
    return pairs


# ------------------------------------------------------------ superposition

def kabsch(mobile: np.ndarray, target: np.ndarray):
    """Optimal rigid-body fit of mobile onto target. Returns (rotation, translation)."""
    mc, tc = mobile.mean(axis=0), target.mean(axis=0)
    P, Q = mobile - mc, target - tc
    V, _, Wt = np.linalg.svd(P.T @ Q)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    R = V @ D @ Wt
    return R, tc - mc @ R


def apply_transform(coords: np.ndarray, R: np.ndarray, t: np.ndarray) -> np.ndarray:
    return coords @ R + t


def tm_score(pred: np.ndarray, ref: np.ndarray, l_target: int | None = None) -> float:
    """
    TM-score for already-superposed, already-paired CA coordinates.

    Not a full TM-align (no search over alternative superpositions), which is
    what you want here: the superposition is chosen deliberately, on the core.
    """
    if len(pred) == 0:
        return 0.0
    L = l_target or len(ref)
    d0 = 1.24 * (max(L, 19) - 15) ** (1.0 / 3.0) - 1.8
    d0 = max(d0, 0.5)
    d = np.linalg.norm(pred - ref, axis=1)
    return float(np.sum(1.0 / (1.0 + (d / d0) ** 2)) / L)


# ------------------------------------------------------------------ results

@dataclass
class StateComparison:
    ref_name: str
    core_rmsd: float
    region_rmsd: float | None
    tm: float
    n_core: int
    n_region: int


@dataclass
class StateRecoveryResult:
    prediction: str
    verdict: str                  # name of the closest reference state
    margin: float                 # region RMSD gap between the two states, A
    comparisons: list = field(default_factory=list)
    ligand: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)

    @property
    def decisive(self) -> bool:
        """A margin under 1 A means the two states are not really distinguished."""
        return self.margin >= 1.0

    def summary(self) -> str:
        lines = [f"prediction: {self.prediction}"]
        for c in self.comparisons:
            reg = f"{c.region_rmsd:6.2f}" if c.region_rmsd is not None else "   n/a"
            lines.append(
                f"  vs {c.ref_name:<12s} core_rmsd={c.core_rmsd:5.2f} A  "
                f"region_rmsd={reg} A  TM={c.tm:.3f}  (n_core={c.n_core}, n_region={c.n_region})"
            )
        flag = "" if self.decisive else "   <-- MARGIN TOO SMALL, not a real call"
        lines.append(f"  verdict: {self.verdict}  (margin {self.margin:.2f} A){flag}")
        for k, v in self.ligand.items():
            lines.append(f"  ligand {k}: {v}")
        for n in self.notes:
            lines.append(f"  ! {n}")
        return "\n".join(lines)


# ---------------------------------------------------------------- the metric

def _paired_coords(pairs, seqid_filter=None, use_ref_seqid=True):
    pred_xyz, ref_xyz = [], []
    for ref_res, pred_res in pairs:
        if seqid_filter is not None:
            key = ref_res.seqid if use_ref_seqid else pred_res.seqid
            if key not in seqid_filter:
                continue
        pred_xyz.append(pred_res.atom("CA"))
        ref_xyz.append(ref_res.atom("CA"))
    if not pred_xyz:
        return np.zeros((0, 3)), np.zeros((0, 3))
    return np.array(pred_xyz), np.array(ref_xyz)


def compare_to_reference(
    prediction: Structure,
    reference: Structure,
    ref_name: str,
    region_seqids: set | None = None,
    core_seqids: set | None = None,
    pred_chain: str | None = None,
    ref_chain: str | None = None,
) -> StateComparison:
    """
    Superpose prediction onto reference (on core_seqids, or everything) and
    report RMSD over the core and, separately, over the state-defining region.

    region_seqids / core_seqids are numbered in the REFERENCE structure's
    numbering, since that is what you can look up in PyMOL.
    """
    pairs = residue_correspondence(reference, prediction, ref_chain, pred_chain)
    if not pairs:
        return StateComparison(ref_name, float("nan"), None, 0.0, 0, 0)

    # Superpose on the core. Default core = everything except the state region,
    # which is the whole point: the region must not drive the alignment.
    if core_seqids is None and region_seqids is not None:
        core_seqids = {r.seqid for r, _ in pairs} - set(region_seqids)

    core_pred, core_ref = _paired_coords(pairs, core_seqids)
    if len(core_pred) < 3:
        core_pred, core_ref = _paired_coords(pairs, None)

    R, t = kabsch(core_pred, core_ref)

    core_fit = apply_transform(core_pred, R, t)
    core_rmsd = float(np.sqrt(((core_fit - core_ref) ** 2).sum(axis=1).mean()))

    all_pred, all_ref = _paired_coords(pairs, None)
    tm = tm_score(apply_transform(all_pred, R, t), all_ref)

    region_rmsd = None
    n_region = 0
    if region_seqids:
        reg_pred, reg_ref = _paired_coords(pairs, set(region_seqids))
        n_region = len(reg_pred)
        if n_region:
            reg_fit = apply_transform(reg_pred, R, t)
            region_rmsd = float(np.sqrt(((reg_fit - reg_ref) ** 2).sum(axis=1).mean()))

    return StateComparison(ref_name, core_rmsd, region_rmsd, tm,
                           len(core_pred), n_region)


def state_recovery(
    prediction: Structure,
    references: dict,
    region_seqids: dict | None = None,
    core_seqids: dict | None = None,
    pred_chain: str | None = None,
    ref_chains: dict | None = None,
) -> StateRecoveryResult:
    """
    The headline call. references maps {state_name: Structure}, e.g.
    {"DFG-in": ..., "DFG-out": ...}.

    region_seqids maps {state_name: set of residue numbers} defining the
    state-sensitive region in each reference's own numbering. If you give the
    same region for both, that's fine and usual.

    Verdict is the state with the lowest region RMSD (falling back to core
    RMSD if no region was given). Margin is the gap to the runner-up -- if it
    is under 1 A, treat the call as undecided and say so out loud.
    """
    region_seqids = region_seqids or {}
    core_seqids = core_seqids or {}
    ref_chains = ref_chains or {}

    comparisons = []
    for name, ref in references.items():
        comparisons.append(
            compare_to_reference(
                prediction, ref, name,
                region_seqids=region_seqids.get(name),
                core_seqids=core_seqids.get(name),
                pred_chain=pred_chain,
                ref_chain=ref_chains.get(name),
            )
        )

    def key(c):
        return c.region_rmsd if c.region_rmsd is not None else c.core_rmsd

    ranked = sorted(comparisons, key=key)
    verdict = ranked[0].ref_name
    margin = (key(ranked[1]) - key(ranked[0])) if len(ranked) > 1 else float("inf")

    notes = []
    if not any(c.region_rmsd is not None for c in comparisons):
        notes.append(
            "no state-defining region given, so this is a global-RMSD call and "
            "is much less sensitive -- define the region in targets.yaml"
        )

    return StateRecoveryResult(
        prediction=prediction.path,
        verdict=verdict,
        margin=margin,
        comparisons=comparisons,
        notes=notes,
    )


# ------------------------------------------------------------ ligand metrics

def ligand_centroid_distance(prediction: Structure, reference: Structure,
                             pred_chain=None, ref_chain=None) -> float | None:
    """
    Distance between predicted and experimental ligand centroids, after
    superposing the protein on its core. Robust where atom-matched RMSD is
    fragile (atom naming differs, symmetry, partial occupancy).
    """
    pred_ligs, ref_ligs = prediction.ligands(), reference.ligands()
    if not pred_ligs or not ref_ligs:
        return None

    pairs = residue_correspondence(reference, prediction, ref_chain, pred_chain)
    if len(pairs) < 3:
        return None
    pred_ca, ref_ca = _paired_coords(pairs, None)
    R, t = kabsch(pred_ca, ref_ca)

    pl = max(pred_ligs, key=lambda r: len(r.heavy_atoms()))
    rl = max(ref_ligs, key=lambda r: len(r.heavy_atoms()))
    pc = apply_transform(pl.heavy_atoms(), R, t).mean(axis=0)
    return float(np.linalg.norm(pc - rl.heavy_atoms().mean(axis=0)))


def pocket_contacts(structure: Structure, cutoff: float = 4.5,
                    chain: str | None = None) -> set:
    """Residue numbers with any heavy atom within `cutoff` of the largest ligand."""
    ligs = structure.ligands()
    if not ligs:
        return set()
    lig = max(ligs, key=lambda r: len(r.heavy_atoms()))
    lig_xyz = lig.heavy_atoms()

    contacts = set()
    for res in structure.chain(chain):
        xyz = res.heavy_atoms()
        if len(xyz) == 0:
            continue
        d = np.linalg.norm(xyz[:, None, :] - lig_xyz[None, :, :], axis=-1)
        if d.min() <= cutoff:
            contacts.add(res.seqid)
    return contacts


def contact_recovery(
    prediction: Structure,
    reference: Structure,
    restrict_to: set | None = None,
    cutoff: float = 4.5,
    pred_chain=None,
    ref_chain=None,
) -> dict:
    """
    How much of the experimental binding-site contact set did the prediction
    reproduce? Independent of ligand RMSD, and the right question to ask when
    testing whether a ligand reached a state-specific sub-pocket.

    Pass restrict_to = the back-pocket residues (in reference numbering) to
    ask specifically: did the ligand occupy the pocket that only exists in the
    target state?
    """
    ref_contacts = pocket_contacts(reference, cutoff, ref_chain)
    pred_contacts_local = pocket_contacts(prediction, cutoff, pred_chain)

    # translate the prediction's contacts into reference numbering
    pairs = residue_correspondence(reference, prediction, ref_chain, pred_chain)
    pred_to_ref = {p.seqid: r.seqid for r, p in pairs}
    pred_contacts = {pred_to_ref[s] for s in pred_contacts_local if s in pred_to_ref}

    if restrict_to is not None:
        ref_contacts &= set(restrict_to)
        pred_contacts &= set(restrict_to)

    if not ref_contacts:
        return {"recall": None, "jaccard": None,
                "n_ref": 0, "n_pred": len(pred_contacts)}

    tp = len(ref_contacts & pred_contacts)
    union = len(ref_contacts | pred_contacts)
    return {
        "recall": tp / len(ref_contacts),
        "jaccard": tp / union if union else 0.0,
        "n_ref": len(ref_contacts),
        "n_pred": len(pred_contacts),
    }
