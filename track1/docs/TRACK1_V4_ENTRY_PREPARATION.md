# Prepare the completed Track 1 V4 entry

The preparation command verifies the completed V4 run and creates a separate `submission_ready/` folder. It preserves the original ranking, prediction bytes, V3 evidence and V4 manifest. It does not repeat genome annotation, ranking or Exomiser.

```bash
.venv/bin/python scripts/prepare_track1_v4_entry.py
```

Outputs stay under `work/private_analysis.nosync/track1_v4/submission_ready/`. The command captures diagnostics in a private log and returns generic status only. Patient findings must not be pasted into a hosted assistant.

## What is checked

- Completed V3/V4 manifests, recorded source and evidence hashes, effective V4 configuration and current input/reference provenance.
- The CSV matches the existing V4 automatic selection byte-for-byte, including both members of every exported pair and the ordinal EPCR values.
- Each exported allele is checked again against the original VCF genotype and normalized GRCh38 reference.
- The completed Exomiser comparison matches the input VCF, exact positive phenotype payload and recorded result hash.
- Source files are copied into a new staging directory. The ZIP uses an explicit inventory and cannot absorb leftovers from an older generated repository.

The report adds deterministic explanations of the leading hypotheses using their recorded evidence. These are automated descriptions. No clinical reviewer, confirmed phase, calibrated probability or independent accuracy result is implied. Neither the ranked hypotheses nor their scores change during preparation.

## Participant details

Edit `submission_ready/entry_identity.json` locally. Required facts are:

| Field | What to supply |
|---|---|
| `team` | The team name or participant identifier for this entry |
| `github_url` | The actual `https://github.com/OWNER/REPOSITORY` URL |
| `ai_details.plan_tier` | The real plan/tier used for Codex |
| `ai_details.data_controls` | Accurate relevant data-handling settings; the existing text records the participant's opt-out statement and the lack of independent verification |
| `ai_details.other_ai_tools` | Other tools used, with their providers/plans/controls; use `None` only if accurate |

Provider `OpenAI` and product `Codex` are prefilled from the known workflow. Correct them or add other providers if necessary. The preparation command creates the disclosure from these fields when they are complete. Nonempty text does not verify account settings or contractual compliance.

Run the same command after editing. To read details from a different local file, use `--identity /path/to/entry_identity.json`. The previous prepared package is archived locally, and a new one replaces it only after the preparation checks pass. Both upload files receive team/model filenames once a team is supplied.

## Files and submission

`predictions.csv` and `methods_submission.md` are the canonical local files. Team-named copies are created when possible. `repository/` and `repository.zip` contain source only. `submission_bundle.zip` contains the CSV, report and source ZIP for local convenience. `START_HERE.md` names the exact upload files. `readiness.json` distinguishes technical checks from unresolved participant and remote checks.

The [official Track 1 form](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026) expects separate prediction and methods uploads plus a GitHub URL. Log in, check the remaining attempts and applicable rules, and use the team-named files. The source repository can remain private during the hackathon and must become public afterward according to the organizer instructions. Publish only the generated source folder; the rest of the local package contains patient-derived findings.

Preparation does not authenticate, publish code, upload patient findings, accept rules or consume a submission attempt. A well-formed GitHub URL does not establish repository accessibility. “V4” is the software version and does not imply four submissions were made.

The original V4 runner's package remains a record of its run. Use this preparation command for final entry packaging and identity changes. No new model run is necessary for those changes.
