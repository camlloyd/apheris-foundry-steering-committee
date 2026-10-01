# apheris-foundry-steering-committee

RORγ Two-State Steering — Hackathon Day 1 Plan (v2, replaces kinase plan)
Date: 2026-10-01, etc.venues Chancery Lane. Team: 2 (Lily + 1 teammate). Hard constraint: Lily leaves at 17:30 (assumed from "I will have to at 5:30:pm") → figure done by ~17:15, handoff package by 17:30, teammate owns the evening.

1. Target and hypothesis
Target: RORγ ligand-binding domain (RORC / NR1F3, UniProt P51449, confirmed, no artifact in the config), ~230-residue LBD construct.
The state switch is helix 12 (H12 / AF-2). Agonist 25-hydroxycholesterol (25-HC, CCD HC2) packs H12 against the core (agonist lock H479–Y502–F506), forming the LXXLL coactivator groove [71, 75]. Inverse agonists destabilize or displace H12 → coactivator (SRC2 / RIP140 / NCOA1 / NCOA2 LXXLL peptide) cannot bind → inflammatory cytokine / Th17 transcription drops (challenge slide).
Critical caveat — state–function decoupling: some inverse agonists co-crystallize with H12 still in the agonist position (diosgenin 9VZQ [70]; "water-trapping" class in 5NTK [62, 79]). Therefore the two reference structures are chosen empirically from CA displacement profiles, not from ligand labels.
Hypothesis: state-specific pocket conditioning (which residues the ligand is asked to contact) steers Boltz-2 / OpenFold3 between H12-in and H12-out, while naive / shared-pocket conditioning does not. Controls are the differentiator: KinConfBench showed mode collapse and pose≠state [9, 22]; Boltz-2 defaults to the dominant state on kinases [5]; apo RORγ LBD crystallizes active [71, 78], so the unconditioned baseline is expected to sit in H12-in.
**Pre-registered success criterion (write before any results):** A4 median H12 displacement is closer to the inactive reference than A3, A5 and A8 by more than the seed spread (headline = recovery rate over 5 seeds). A null result is reported as a finding. Written into the config and runbook before any prediction ran.
**Threats to validity:** (a) contact-set leakage: pocket sets derive from the same crystals as the references, so they may encode ligand identity rather than H12 state; mitigated by the A2 / A4 controls; (b) template memory: 3KYT / 4ZJW are pre-2025 and may sit in training data or templates; mitigated by the wrong-template (A6) and MSA-free arms, plus any post-cutoff pharma-drop data; (c) the auto-derived H12 region is unvalidated until checked on crystals (Gate 1, closed, see §4).
2. Reference pair (final choice made empirically in execution step 1)
ACTIVE reference: 3KYT — RORγ LBD + 25-HC (HC2) + SRC2-2 LXXLL peptide (chain B), canonical agonist conformation [72, 74, 78].
INACTIVE candidates: 5IXK (BIO399, AF-2 destabilized) [63, 64]; 4ZJW (4P1) and 4ZJR (4P3) biaryl carboxylamides, H479 H-bond key [68, 69]; a steric-clash-class member of the 5NTK set [62, 79]; MRL-871 allosteric structure (unprecedented H12 position) [67] as backup.
Selection procedure (kit does this): download all candidates → state_recovery2 CA displacement profile of each vs 3KYT → pick the candidate with the largest, cleanest H12-localized displacement; the state region is auto-derived (threshold max(2 Å, median), segment merging). If the top two candidates disagree by <1 Å on the region, flag the call as non-decisive and fall back (§7).
3. Experimental grid (adapted from A0–A8 / B0–B3)
Primary scorer: state_recovery_generic(pred, {"active": 3KYT, "inactive": chosen_ref}) on the auto-derived H12 region. Secondary readouts: ligand pose RMSD, pocket-contact recovery.

Arm	Protein	Ligand	Conditioning	Purpose
A0	LBD	25-HC	none	baseline — expect H12-in (apo RORγ crystallizes active)
A1	LBD	25-HC	active pocket set	steering → active
A2	LBD	25-HC	inactive pocket set	cross / conflict control
A3	LBD	inverse agonist	none	baseline — where does the model default?
A4	LBD	inverse agonist	inactive pocket set	steering → inactive (the money arm)
A5	LBD	inverse agonist	active pocket set	cross / conflict control
A6	LBD	25-HC	wrong template (kinase ref)	negative control
A7	LBD	– (apo)	none	apo control
A8	LBD	25-HC	shared-pocket-only	control: pocket constraints alone ≠ state steering
B3 (optional)	LBD + SRC2-2 12-mer	25-HC	active pocket set	LXXLL peptide co-fold; peptide-contact recovery as orthogonal state readout (slide features the coactivator explicitly — good for the pitch)
Inverse agonist ligand: SR2211 (SMILES in hand [73]) or 4P1 / 4P3 once the CCD codes are confirmed from the downloaded refs.
Pocket sets: ligand-contact residues (≤4 Å) computed from each reference with gemmi; the inactive set is expected to include back-pocket residues near S404 and H12-facing contacts.
Core arms (priority): A0, A1, A2, A3, A4, A5, A8. A2 (shared-pocket-style cross control) is the cheapest arm and the only clean test that pocket conditioning per se does not flip the state, so it stays core. A6, A7 only with spare quota.
Seeds: 5 on every inverse-agonist core arm (the headline is a recovery rate; n=2 cannot support one), 2 on the 25-HC direction. Full plan 49 predictions; quota-tight floor 27. Submit the entire grid as ONE batch.
Model: single-model-first. Run the whole grid on Boltz-2 (or whichever model Foundry exposes pocket constraints for); add OpenFold3 mirrors only after the grid is done.
4. Execution mechanics
**Go/no-go gates.** (1) Scorer validation (sandbox, before 11:30) — CLOSED: the scorer called all five crystals correctly and decisively: 3KYT→H12-in, 4ZJW→H12-out, 5IXK→H12-out, 5NTK→H12-out, and 4ZJR→H12-in (margin 2.33 Å), the inverse-agonist decoy whose H12 is still packed. The scorer handles the state–function decoupling case honestly, which is itself a pitch slide. (2) Schema gate (by 12:15, teammate, via foundry_emit.py probe on the laptop with Hub access): if Foundry exposes no pocket constraints, pivot the conditioning lever (templates, or a ligand with known H12 preference) instead of patching; the HPC Boltz-2 fallback must also be checked for constraint support.
Foundry: apheris-foundry CLI v0.5.0 — workflows run --workflow predict --input request.json --model-params @file, jobs logs / download [52]. Emit every arm with foundry_emit.py (probe / emit / collect subcommands). Exact request/params schema to be captured at Q&A (11:05–11:15) or from organizers at lunch. Fallback: classic compute-spec Docker path [44–47].
Optional de-risk (ask organizers it's OK): 2–3 smoke predictions via sandbox HPC Boltz-2 to confirm the state region and pocket sets behave before committing Foundry quota.
All scoring is client-side with the kit: state_recovery2 + score_run.py + figure.py. No kinase motif machinery needed (DFG logic drops out; motif_overrides stays unused).
5. Timeline (minute-by-minute, Day 1)
Now–11:30 (I work in the sandbox in parallel with the intros): download RORγ refs, run displacement profiles, pick the pair, derive region + pocket sets, build targets_rorgamma.json, smoke-test the emitters, package rorgamma_addon.zip + runbook. — Lily: capture Foundry/OpenFold schema details, GPU quota, and any pharma-data drop format.
11:30–12:15 team formation: register the team, confirm roles.
12:15–13:00 (lunch): validate request.json against the captured schema; submit one smoke arm; fix issues.
13:00–13:45 (Build I): submit the full grid as one batch; monitor.
13:45–15:30: collect finished jobs as they land; score; interim figure.
15:45–17:00 (Build II): complete scoring; B3 peptide arm if time; sensitivity check (state region ±5 residues).
17:00–17:15: final figure_day1.png + summary_table.md.
17:15–17:30: handoff package to teammate (configs, results.csv, figure, evening runbook).
Evening (teammate): extend arms/seeds, polish figure, draft the 5-min pitch. Lily remote-optional.
6. Task division — Science Lead (Lily) vs Foundry Ops (teammate)
6.1 Science Lead (Lily, until 17:30)
S1. Now–11:30 — references (delegated to the sandbox): trigger the reference pipeline: download 3KYT / 5IXK / 4ZJW / 4ZJR / one 5NTK member (+ MRL-871 backup), run state_recovery2 displacement profiles, pick the pair, auto-derive the H12 region, compute active/inactive pocket sets, build targets_rorgamma.json, smoke-test the emitters, package rorgamma_addon.zip + runbook. Lily reviews the displacement-profile output for 5 min when it lands.
S2. Now–11:30 — challenge intel: from the challenge talk: exact judging criteria, expected deliverable, whether a pharma/ligand data drop exists (format, deadline).
S3. 11:05–11:15 Q&A (backup if the teammate is mid-setup): ask the schema questions F3.
S4. 11:30–12:15: team registration; pick the inverse-agonist ligand (SR2211 vs 4P1/4P3) once CCD codes are confirmed; approve the reference pair.
S5. 12:15–13:00 (lunch): review the smoke-arm request.json; freeze the grid; set scoring thresholds (margin <1 Å = non-decisive call).
S6. 13:00–15:30 (Build I): score predictions as the teammate downloads them (state_recovery2 → results.csv); flag failing arms early; decide extensions.
S7. 15:45–17:00 (Build II): sensitivity check (state region ±5 residues); score the B3 peptide arm if submitted; interpret the control arms (A2 / A5 / A6 / A8).
S8. 17:00–17:15: final figure_day1.png + summary_table.md; write the pitch outline.
S9. 17:15–17:30: handoff (§6.3).
6.2 Foundry Ops (teammate, owns the evening)
F1. Now–11:05: install/verify apheris-foundry CLI v0.5.0, authenticate, run foundry_emit.py probe to capture the live request schema.
F2. 10:35–11:05 OpenFold intro: capture the model-params schema — especially the pocket-constraint syntax — plus per-job runtime and whether Boltz-2, OpenFold3, or both are exposed.
F3. 11:05–11:15 Q&A — the four questions: (a) exact request.json / model-params schema incl. pocket constraints; (b) GPU quota and max concurrent jobs per team; (c) pharma/ligand data drop — where and in what format; (d) when submissions close for Day-2 judging.
F4. 11:30–12:15: team registration logistics; make sure BOTH members have Foundry workspace access.
F5. 12:15–13:00 (lunch): submit one smoke arm (A0); verify end-to-end (submit → run → download → feeds score_run.py).
F6. 13:00–13:45 (Build I): submit the full grid as ONE batch; monitor; pull logs on failures; resubmit.
F7. 13:45–17:00: download finished predictions into preds/<arm>/<seed>/; keep the job tracker current.
F8. Evening (owns it): run the extension queue Lily left; rescore + regenerate the figure; polish slides; rehearse the 5-min pitch; Day-2 demo prep.
6.3 Handoff (17:15–17:30) and shared decisions
Handoff package: rorgamma_addon.zip (targets_rorgamma.json, emitter grids, scoring scripts), results.csv so far, figure_day1.png/.svg + summary_table.md, evening runbook (prioritized extension queue, rescore + figure commands), pitch outline, open risks.
Shared: grid freeze at ~12:45 (both); pitch storyline agreed before 17:30; check-in channel for evening questions.
Swap roles if the teammate is the stronger structural biologist — the checklists transfer 1:1.
7. Fallbacks
Foundry schema unknown or blocked → classic Docker compute-spec path; if Foundry is entirely unavailable → run the grid via sandbox HPC Boltz-2 and state so honestly in the pitch.
No clean inactive reference (all candidates keep H12-in) → use the MRL-871 allosteric structure [67] or an apo-vs-bound contrast; report the decoupling itself as a finding — it is insight, not failure.
GPU quota tight → cut to the 27-prediction floor (core arms only, drop A6 / A7 / B3, reduce seeds on the 25-HC direction first).
Peptide arm fails or crowds the schedule → drop B arms; peptide recovery is a bonus readout only.
8. Deliverables
figure_day1.png / .svg + summary_table.md; results.csv; targets_rorgamma.json; rorgamma_addon.zip (new filename — the results mount forbids overwriting); evening handoff runbook; pitch outline.

9. Pitch line
"State-specific pocket conditioning steers co-folding models between RORγ agonist and inverse-agonist H12 states — and the controls show naive conditioning can't." Novelty vs prior art: a controlled head-to-head of ligand-induced state recovery with state-specific (not shared) pocket conditioning on a nuclear receptor, extending the kinase DFG-era benchmarks [5, 9, 12, 22].

10. Compute estimate
Full plan: 49 predictions (quota-tight floor: 27); a ~230-residue LBD plus a small ligand runs in minutes per job on the Lyceum GPUs — fits comfortably inside Build I–II. Client-side scoring is <1 min total.
