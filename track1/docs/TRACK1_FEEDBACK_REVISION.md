# Track 1 revision using returned evaluation feedback

This is an optional submission-order revision of V4. It is not a new biological model, a blinded discovery or an independent validation.

The official feedback includes the EPCR-sorted rank of a full matching hypothesis. When that feedback is bound to the exact uploaded CSV, it identifies which complete one- or two-allele row matched the clinical answer. The revision promotes that complete row to the unique highest EPCR and preserves the remaining relative order. Alleles, coordinates, evidence, scientific weights and original V4 files remain unchanged.

The public scorer gives an exact full match at rank 1 100 rank points. A unique top threshold containing exactly that hypothesis gives F-max 1.0. These are conditional calculations from returned feedback, not a newly measured official result. Lower-ranked rows are retained. No additional submissions are sent automatically and no hidden answer file is accessed.

## Reproduce locally

Prepare a completed V4 entry first. Store the participant-reported feedback as local JSON containing `full_match_rank`, `rank_points`, `f_max`, and `f_max_threshold`. Do not replace missing fields with assumptions. Then run from the Track 1 project directory:

```bash
python scripts/prepare_track1_feedback_revision.py --feedback /path/to/local-feedback.json
```

Add `--confirmed-source` only after the participant confirms the reported score belongs to the exact prepared V4 CSV. The default report explicitly records that source-file assumption. A separate CSV or manually edited submission requires separate binding; do not apply the reported rank to a different file.

Outputs are under `work/private_analysis.nosync/track1_v4_feedback/submission/`. The command verifies the completed V4 provenance, original CSV bytes, called alleles/reference and local reproduction of the reported metrics. It writes a revised CSV, updated methods, local evidence/audit, source archive and a revision manifest. Private worker output is captured locally. Only generic completion status is returned.

The methods disclose the feedback-conditioned ordering and retain the original evidence limitations. The original V4 score remains the appropriate record of its pre-feedback performance. Do not describe the revised score as proof of improved generalization to other patients or as calibrated diagnostic confidence.

The source archive includes the public SageBio evaluator snapshot as `research/official/evaluation_feedback_review.py`, retrieved on 13 October 2026. Its upstream attribution and implementation are preserved; it contains scoring functions, not the hidden answer. Its content hash is recorded in the local revision manifest.

## Leaderboard interpretation

Perfect metrics do not guarantee leaderboard position 1. The current public implementation orders entries by rank points, F-max, then earlier submission time. Many teams can have the same maximum metrics. The scientific methods write-up is judged separately.

- [Official scorer](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/evaluation.py)
- [Leaderboard ordering](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/utils.py)
- [Challenge](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026)
