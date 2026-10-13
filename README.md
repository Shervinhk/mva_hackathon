# MVA Hackathon

**Team: Kinetochore Compass**

Research code for the MVA Hackathon 2026. Each track has its own directory and instructions.

## Tracks

- [Track 1: variant prediction](track1/README.md) — the V4 approach, built on reusable V3 evidence, with reproducible source, configuration and fabricated tests.
- [Track 2](track2/README.md) — reserved for future work. No Track 2 model or results are included yet.

## Start with Track 1

Run Track 1 commands from its own directory:

```bash
cd track1
python -m venv .venv
source .venv/bin/activate
pip install -r requirements-track1.txt
python -m unittest discover -s track1_v4/tests -v
```

Read the [V4 method](track1/track1_v4/README.md) and [entry preparation guide](track1/docs/TRACK1_V4_ENTRY_PREPARATION.md) before running an analysis. V4 reuses a completed V3 analysis; authorized inputs and public references must be obtained separately. Moving the project into `track1/` does not change its internal relative directory layout.

## Data and interpretation

This repository contains source code, public methodology and fabricated test examples. Patient inputs, variant predictions, case-specific reports, logs and reference databases are excluded. The prediction CSV and case methods report are submitted separately through the authorized challenge form.

The output is a research shortlist. Software validation does not establish a clinical diagnosis, phase confirmation or calibrated probability. See [AI use and privacy](track1/docs/AI_USE_AND_PRIVACY.md) for the workflow's disclosure history and limitations. Participant-specific AI plan/tier details belong in the separately prepared submission report.

The repository can host both tracks. Their code, documentation and entry files remain separate so reviewers can identify the method for each submission.
