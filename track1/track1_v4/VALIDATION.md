# V4 validation record

The development workspace passed all 54 V4 fabricated regression tests and all 153 existing V3/workflow regression tests. Tests use synthetic variants, toy ontology data, fabricated Word documents and temporary manifests; patient findings are not returned to the assistant.

V4 checks cover negative and family observations, ambiguous and contradictory scope, HPO aliases and information-content bounds, ClinVar review precedence, X-linked context, transcript compatibility, stale-cache refusal, original-validation requirements, pair export, duplicate allele sets, score ties, participant placeholders, local logging and official pair/F-max mechanics.

The official scorer regression uses its public implementation and fabricated truth only. That check skips in a standalone source archive where the separately acquired public scorer is absent. This does not constitute an independent clinical cohort, calibration, or a leaderboard evaluation.

Actual execution completion is recorded locally in `work/private_analysis.nosync/track1_v4/stage_status.json`. The V4 manifest binds the reused V3 evidence, relevant inputs/references, source and output hashes. Detailed findings and V3/V4 differences are in that private directory. A generated file or a passing software test does not establish biological causality.

The full local V4 run completed, including a matching Exomiser comparison and a validated prediction CSV. Final checks confirmed V4 source/output hashes, correct bundle contents, original-call/reference validation, and unchanged V3 manifest, evidence and prediction CSV. The source archive passed its standalone V4 checks; the three official-scorer tests were skipped as one class because the public scorer snapshot is deliberately acquired separately. These tests passed in the full workspace.

No patient inputs or results were included in the source archive. A local check confirmed case-specific input-path metadata was absent from its contents. No upload or challenge submission occurred. Team, GitHub destination and actual AI plan/disclosure fields remain participant-controlled and unresolved until supplied.

## Entry preparation

The separate entry-preparation command has 18 passing fabricated checks for provenance, exact CSV/model binding, fresh source staging, archive paths and missing participant details. The local rationale renderer has five passing fabricated checks for grounded evidence reporting, uncertain phase, missing population evidence, complete-hypothesis sensitivity matching and markup handling. These additions do not change the recorded V4 ranking implementation or predictions.

The prepared entry is in `work/private_analysis.nosync/track1_v4/submission_ready/`. Its separate preparation record binds the original analysis manifest, unchanged prediction bytes, revised methods and current source archive. Participant details and remote submission state remain separate from technical preparation.
