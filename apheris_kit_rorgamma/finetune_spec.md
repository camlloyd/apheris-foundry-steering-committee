# Stretch arm: federated fine-tuning spec (execute ONLY if all inference arms are scored)

**Owner:** ML teammate. **Earliest start:** Day 1, 17:00, and only if A0–A8 are scored.
**Why it's the stretch:** federated fine-tuning on proprietary structures is Apheris's
flagship capability (their demonstrated pattern is federated LoRA — FRA-LoRA on ESM-2
across the AISB pharma consortium). Showing even a toy version of it on their own
platform is worth more than another inference arm. Showing a *failed* half-run is
worth nothing — hence the gate.

## Question

Does fine-tuning on state-annotated proprietary structures recover the target state
*without* any inference-time steering — and does it beat the best inference lever?

## Data prep (from triage.py output)

1. From the pharma drop, take every structure of the target protein(s) in the
   **target state** (triage labels: DFG-out / rare cluster). Expect single digits —
   that is fine for LoRA, fatal for full fine-tuning. If fewer than 5 examples,
   do not start; write it up as "future work" instead.
2. Format as training samples in whatever the Hub's fine-tune workflow advertises
   (probe: `apheris-foundry workflows list` — look for a `finetune`/`train` workflow).
   If only the classic Gateway path exists (custom Docker image + Compute Spec),
   skip: image build + governance approval will not fit in the evening.
3. Hold out ONE target-state structure as the eval case. Never eval on a training
   example — that is the memorisation trap in miniature.

## Configuration sketch

- **Method:** LoRA on the diffusion/transformer attention projections (query/value
  first), rank 8–16, alpha = 2×rank, dropout 0.05.
- **Optimiser:** AdamW, lr 1e-4, cosine decay, 500–1000 steps, batch by 512-token
  crops. On a single Lyceum GPU this is 1–3 hours — overnight, not afternoon.
- **Federation:** if the structures sit behind two different custodians, use the
  platform's federated averaging (FRA-LoRA-style full-rank aggregation of adapters).
  If one custodian, plain single-site LoRA — say so honestly.
- **Freeze everything else.** The trunk keeps general folding; the adapter learns
  the state bias.

## Evaluation protocol (identical to the inference grid)

Run the fine-tuned model on the held-out target with **default settings** (full MSA,
no template, no pocket constraint) and score with `score_run.py score`. Success =
target-state recovery at default settings where the baseline failed. Anything less
(e.g. needs a template anyway) is reported as-is.

## Failure modes to state up front

- **n is tiny:** LoRA on <10 structures learns the specific complex, not the state.
  Mitigation: include the public two-state pairs (ABL1, p38α) in training and hold
  out the pharma structure for eval — then the claim is "adapts to the state class".
- **Catastrophic forgetting:** check one unrelated fold (any small protein) after
  training; if it degrades, say so.
- **Governance latency:** if the Compute Spec isn't approved by 21:00, stop and put
  the spec on the "what we'd do next" slide. A credible plan is a deliverable.
