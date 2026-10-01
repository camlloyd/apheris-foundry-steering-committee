# Handoff message — Lily → Foundry Ops teammate (Day 1 evening)

*Ready to paste; edit names/times as needed.*

---

Status going into the evening session — science side is in good shape, one scoring gap to close tonight.

**Done (my S1–S9 list):**
- S1 references: 3KYT/4ZJW pair picked and in use; H12 region derived; pocket sets computed. Can you confirm `rorgamma_addon.zip` + runbook are actually packaged on the VM? That's the one piece I can't see from here.
- S2 challenge intel: judging criteria confirmed from the deck — **60% prediction accuracy / 30% scientific approach / 10% reproducibility**. Presentations 10:30 tomorrow.
- S4/S5: 4P1 chosen as the inverse agonist; grid frozen; 1 Å margin threshold in force.
- S6: H1 scored — 0/122 structurally sound predictions reached H12-out across nine conditioning levers.
- Since then: H2 titration (fold-stability cliff between 32 and 128 MSA sequences; 0/200 H12-out; pre-registered negative), H3 (four arms — apo, 25-HC, 4P1, BIO592 — H12-region confidence partitions by pharmacological class, p = 0.0079, confound arm clean), and all three follow-ups (cliff bisection, template×depth, coactivator co-fold) are complete. Please pull the follow-up numbers into `H3_RESULTS_REPORT.md` tonight.

**One integrity fix you should know about before the slides:** the kit's "H12 window" (479–486) is the H11–H12 junction, not the challenge's H12 (501–507). I verified against the PDB files: 4ZJW models 479–486 but not 501–507; 3KYT models both. Our results stand (the junction is state-informative), but tonight we re-extract per-residue pLDDT over 501–507 from the existing H3 outputs — zero GPU, ~10 min — and report both windows. Details in `report_additions.md`, section 1.

**Tonight (GPU):** full-pool generation on the deposited set (~1,400 compounds, ~3h per the deck) — approved to run overnight. Please check the deposited set first: does it include H12-state labels? That decides whether we can validate the classifier or submit blind. Also confirm the GPU is free so we don't collide.

**Tomorrow 09:00–10:15:** score the pool with compute-metrics → run `h12_state_classifier.py` (calibrate on the H3 arms first: `--calibrate`) → state calls per compound → final figure (`make_figure_day1.py --pool-csv ...`) → slides. I'll take the slides and the report; can you own the reproducibility package (pre-registrations, scripts, CSVs, logs, compute usage in one runnable folder)?

---

## Slide outline — "Steering Committee — RORγ H12" (7 slides, mapped to the rubric)

1. **The challenge, restated + our one-line answer.** "322 structures, three independent approaches, zero H12-out: the model cannot be steered into the therapeutic state — but its confidence reads it. We turned a pre-registered negative into a state classifier."
2. **The cliff (H2).** Figure panel a. Fold or agonist, nothing between; declared negative, pre-registered. *(30% — failures analysed.)*
3. **Nine levers, zero exits (H1 + follow-ups).** 0/122 + bisection + template×depth; the prior locks the agonist state whenever the domain folds. *(30%.)*
4. **H3: confidence partitions by pharmacological class.** Figure panel b. Contrast passed (p = 0.0079, complete separation); C3 kills the rigidity confound; ordering prediction failed and we say so; confidence–accuracy decoupling (best geometry, lowest confidence). *(30% + £50 novelty.)*
5. **The classifier on the pool.** Figure panel c + accuracy vs deposited labels (or honest blind-call statement). State the ceiling: state calls come from the confidence signal; coordinates never adopt H12-out. *(60%.)*
6. **The coactivator co-fold** (if the result is clean): SRC2 placed in the groove with an inverse agonist bound — a structure that cannot exist biologically. *(£50 novelty visual; drop if ambiguous.)*
7. **Reproducibility + integrity.** Runnable workflow folder (pre-registrations, scripts, data, compute logs); the 479–486 vs 501–507 numbering note, stated before a judge asks. *(10%.)*
