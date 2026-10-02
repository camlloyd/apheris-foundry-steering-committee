# Integrity note — H12 numbering, and what the references actually resolve

> **⚠ Superseded in part by `results/OFFICIAL_RESCORE_FINDINGS.md`.** This
> note correctly predicted the 501–507 window would be unmeasurable against
> 4ZJW — but the actual fix wasn't "use our 479–486 window instead," it was
> "use the organizers' own reference pair (5VB7/6T4I)," which *does* resolve
> 484–507 in full. Read that file before concluding 479–486 is the right
> window to report on a slide.

Raised as a risk before the final slide: the challenge deck appears to
reference an H12 window around residues 501–507, while every H1/H2/H3 result
in this project uses an auto-derived state region of 479–486 (construct/3KYT
numbering). Checked both halves of this before writing it up.

## 1. Numbering itself is correct — no construct/UniProt offset bug

The construct used throughout (full RORγ LBD, 243 residues) is numbered
identically to 3KYT's deposited author numbering, residues 265–507. Verified
directly against the construct sequence at the three residues HANDOFF.md
names as the agonist-lock triad:

```
residue 479 -> H  (His479)
residue 502 -> Y  (Tyr502)
residue 506 -> F  (Phe506)
```

All three match the literature description exactly (`H479-Y502-F506`). There
is no offset error between the construct numbering used in this project and
the canonical/UniProt numbering the deck's slide is drawn from — both halves
of the project are talking about the same residues by the same numbers.

## 2. The real issue: 4ZJW does not resolve residues 502–507 at all

Checked what the H12-out reference structure (4ZJW, used throughout as the
inactive/inverse-agonist reference) actually covers:

```
4ZJW chain A resolved range: 268–501
Residues 501–507 resolved?  only 501.  502, 503, 504, 505, 506, 507: absent.
```

If "the H12 window" is read as residues 501–507 (the C-terminal portion
containing Tyr502/Phe506), **4ZJW has real coordinates for exactly one of
those seven residues.** Any RMSD-to-4ZJW computed over that window would
either fail outright (missing atoms) or silently degrade to a near-single-
residue comparison — not a meaningful 7-residue segment RMSD. This is not a
coding bug; it is a genuine, unavoidable fact about what this crystal
structure resolved (residues 487–500 are a documented gap, plausibly the
actual physical signature of H12 becoming disordered in the inverse-agonist-
bound state rather than adopting some other fixed pose — see
`H2_PREREGISTRATION.md`'s discussion of the same gap in the template×depth
follow-up).

## 3. Why 479–486 is the defensible choice, not a workaround

The 479–486 region used throughout H1–H3 was not picked to dodge the 4ZJW
gap — it was auto-derived *before* any prediction ran, from the CA
displacement profile between 3KYT and 4ZJW at a 2 Å threshold (see
`targets_rorgamma.json`'s `_state_region_note`). It happens to also be the
only part of H12 that both reference structures resolve in full, which is
why it is usable as a scoring region at all. The kit's own pre-registered
sensitivity check (testing 479–486 against two alternatives — the region
plus loop 286–290, and the full 470–507 window) found all three gave the
same crystal state calls, so the choice is not fragile — but any claim about
"H12" on the slide should say explicitly which window is meant, since 479–486
and "501–507" are genuinely different, non-overlapping segments of the same
helix, and only one of them has usable reference coordinates.

## Recommendation for the slide / report

State explicitly: *"H12 state calls use the auto-derived, pre-registered
segment 479–486 — the only part of helix 12 resolved in both reference
structures. The deck's broader H12 window (through residue 507) includes
residues 4ZJW does not resolve at all (502–507), so state-calling against
that wider window is not possible with this reference pair."* This is a
finding worth stating plainly, not a gap to paper over — it is itself
evidence for the project's recurring theme (H12 destabilization in the
inverse-agonist state shows up as literal absence of coordinates, not a
cleanly alternative pose).
