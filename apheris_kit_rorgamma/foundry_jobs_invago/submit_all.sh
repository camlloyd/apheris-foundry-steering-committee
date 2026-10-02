#!/bin/bash
# Submit the whole grid as ONE batch. Do not trickle jobs through a gateway.
set -e

# A0_baseline: default inference, full MSA, no steering
apheris-foundry workflows run --workflow predict --input foundry_jobs_invago/A0_baseline.request.json --model-params @foundry_jobs_invago/A0_baseline.params.json

# A1_pocket_specific: pocket conditioned on state-SPECIFIC contacts (the hypothesis)
apheris-foundry workflows run --workflow predict --input foundry_jobs_invago/A1_pocket_specific.request.json --model-params @foundry_jobs_invago/A1_pocket_specific.params.json

# A2_pocket_shared: NEGATIVE CONTROL: pocket conditioned on contacts shared by both states -- should not flip the state
apheris-foundry workflows run --workflow predict --input foundry_jobs_invago/A2_pocket_shared.request.json --model-params @foundry_jobs_invago/A2_pocket_shared.params.json

# A3_template_target: template of the target state
apheris-foundry workflows run --workflow predict --input foundry_jobs_invago/A3_template_target.request.json --model-params @foundry_jobs_invago/A3_template_target.params.json

# A4_template_wrong: NEGATIVE CONTROL: template of the WRONG state -- if this also lands on the target state, the model is ignoring the template and you are measuring memorisation
apheris-foundry workflows run --workflow predict --input foundry_jobs_invago/A4_template_wrong.request.json --model-params @foundry_jobs_invago/A4_template_wrong.params.json

# A5_msa_free: MSA-free: removes the evolutionary prior that pulls the model to the dominant state
apheris-foundry workflows run --workflow predict --input foundry_jobs_invago/A5_msa_free.request.json --model-params @foundry_jobs_invago/A5_msa_free.params.json

# A6_pocket_plus_template: state-specific pocket + target template
apheris-foundry workflows run --workflow predict --input foundry_jobs_invago/A6_pocket_plus_template.request.json --model-params @foundry_jobs_invago/A6_pocket_plus_template.params.json

# A7_pocket_msa_free: state-specific pocket, MSA-free
apheris-foundry workflows run --workflow predict --input foundry_jobs_invago/A7_pocket_msa_free.request.json --model-params @foundry_jobs_invago/A7_pocket_msa_free.params.json

# A8_conflict: DIAGNOSTIC: target-state pocket vs WRONG-state template. Which lever wins tells you their relative strength
apheris-foundry workflows run --workflow predict --input foundry_jobs_invago/A8_conflict.request.json --model-params @foundry_jobs_invago/A8_conflict.params.json

# B0_of3_baseline: OpenFold3 default
apheris-foundry workflows run --workflow predict --input foundry_jobs_invago/B0_of3_baseline.request.json --model-params @foundry_jobs_invago/B0_of3_baseline.params.json

# B1_of3_template: OpenFold3 with target-state template (CIF direct mode)
apheris-foundry workflows run --workflow predict --input foundry_jobs_invago/B1_of3_template.request.json --model-params @foundry_jobs_invago/B1_of3_template.params.json

# B2_of3_msa_free: OpenFold3 MSA-free
apheris-foundry workflows run --workflow predict --input foundry_jobs_invago/B2_of3_msa_free.request.json --model-params @foundry_jobs_invago/B2_of3_msa_free.params.json

# B3_of3_template_msa_free: OpenFold3 target template, MSA-free
apheris-foundry workflows run --workflow predict --input foundry_jobs_invago/B3_of3_template_msa_free.request.json --model-params @foundry_jobs_invago/B3_of3_template_msa_free.params.json
