# Track 1 V4

V4 is a separate extension of the completed V3 variant prioritizer for Track 1. It reuses V3's normalized candidate alleles, annotations, population evidence and clinical assertion records. No genome download, read alignment, database indexing or model training is repeated.

## What changes

- Clinical observations retain subject, present/absent/unknown state, source row and available onset text. Ambiguous or contradictory observations are excluded from the positive proband query. Ordinary unqualified clinical table rows remain provisional positive proband observations; this assumption and the source text are retained locally. This is a conservative parser, not a claim of clinical curation.
- Phenotype similarity excludes terms without usable information content from both sums and reports unscored terms. Unsupported descendants no longer suppress a usable parent observation. Alternate HPO identifiers are resolved before handling contradictions. Scores must lie within zero and one.
- X-linked dominant, recessive and unspecified inheritance remain distinct. Genotype, ploidy and GRCh38 pseudoautosomal context determine whether the hypothesis has sufficient context for the automated shortlist. Unknown context is preserved for review; no sex or parental origin is inferred.
- Matching ClinVar expert-panel and guideline assertions take precedence over ordinary submissions. Lower-review dissent remains visible, and conflicts at the active level remain unresolved.
- A consequence outside the supplied disease-specific mechanism rules is deferred when it has no matching condition-specific clinical support. This means insufficient compatibility evidence, not proof that the allele is benign. No new mechanism bonus is added.
- Recessive pairs require a shared plausibly affected transcript annotation. Different selected transcripts remain flagged even if a shared relevant annotation exists. This does not establish phase, functional effect or nonsense-mediated decay.
- Automatic shortlist eligibility, evidence review and clinical confirmation are distinct states. Full compact sensitivity rankings and a V3–V4 comparison make ordering changes inspectable.

Numeric evidence weights, quality criteria and population cutoffs are inherited from the recorded V3 configuration. The changes correct evidence interpretation rather than tune this case to a known answer. Negative observations are retained and excluded from the positive query; V4 does not yet use a negative-phenotype likelihood model. NMD, penetrance, complete structural-variant analysis and transcript function are not newly established by these corrections.

## Folder layout

This folder contains public source and fabricated tests. V3 source and outputs remain separate.

| Location | Contents |
|---|---|
| `track1_v4/` | V4 code, configuration, description and tests |
| `scripts/`, `config/` | Reused V3 implementation and public configuration |
| `work/private_analysis.nosync/track1_v3/` | Completed V3 evidence used as the baseline |
| `work/private_analysis.nosync/track1_v4/` | New V4 analysis, structured phenotypes, sensitivity and comparison |
| `work/private_analysis.nosync/track1_v4/submission/` | New CSV, methods, source archive, local bundle and readiness report |

Patient data, derived findings and detailed logs stay in the excluded private tree. Do not return them to a hosted assistant. The public wrapper redirects Python, native-library and subprocess output into its private execution log and emits generic completion status only. Sync exclusion was verified for the parent folder on this computer; verify it before using another computer or destination.

## Run using the completed V3 analysis

From the project root, using its existing environment:

```bash
.venv/bin/python track1_v4/run_local.py
```

V4 checks the completed V3 manifest, all recorded evidence hashes, current input VCF and clinical-document hashes, GRCh38 FASTA, used HPO references and shared V3 source before reusing data. It refuses a mismatched cache. The original genotype and reference allele of each exported prediction are checked again.

Exomiser is reused only if completion, recorded VCF provenance, reference verification, result integrity and the exact positive HPO payload match. If the corrected input changes that payload, the existing local Exomiser wrapper runs in a separate V4 directory using already downloaded references. No new data preparation is requested. A comparator failure is disclosed and does not masquerade as a completed comparison.

For a deliberately skipped comparison, pass `--skip-comparator`. To rerun it explicitly, pass `--refresh-comparator`. These options do not change the evidence ranking algorithm.

Before each rerun, prior V4 submission files are moved into `previous_submissions/`. Participant details are retained, but an older CSV does not remain at the active submission path after an interrupted run. V3 files are never overwritten.

## Structured phenotype corrections

The generated `phenotypes.json` is a local audit record. Its observations include source text and must remain private. To provide curated corrections, save a separate JSON file containing `schema_version: 1`, the exact clinical-document `source_sha256`, an explicit `reviewed_locally` boolean, and a nonempty `observations` list. Each observation requires `hpo_id`, `subject` (`proband`, `family`, or `unknown`), `present` (`true`, `false`, or `null`), and a positive integer `source_row`. Optional fields include `onset`, `label` and `source_text`. Retain uncertainty when the source does not settle the observation.

```bash
.venv/bin/python track1_v4/run_local.py \
  --phenotype-override work/private_analysis.nosync/track1_v4/curated_phenotypes.json
```

No automatic run claims that this human review occurred. A zero usable-positive-phenotype result stops the run with a local record for correction rather than inventing symptoms.

## Track 1 output and interpretation

The CSV contains up to ten distinct eligible allele hypotheses in model order. A recessive pair occupies one row. EPCR values are strictly ordinal `1/(entry rank+1)`, not calibrated causal probabilities; upstream utility ties are disclosed in the comparison. The report explains that the formal EPCR label is probabilistic but this implementation supplies ordering information only.

The official scorer gives complete pair matches precedence over partial ones, so adding rows is not universally harmless for rank points. Fabricated tests cover pair, overlap and tie semantics. No hidden answer or leaderboard probing is used to choose predictions. Clinical accuracy and improvement over V3 require external validation or official evaluation; neither is inferred from software tests or tool agreement.

The original run package includes `predictions.csv`, `methods_submission.md`, `comparison_v3_v4.md`, `entry_identity.json`, `readiness.json`, `repository.zip` and `submission_bundle.zip`. For final packaging and participant-detail changes, use the separate [entry preparation command](../docs/TRACK1_V4_ENTRY_PREPARATION.md):

```bash
.venv/bin/python scripts/prepare_track1_v4_entry.py
```

It creates `work/private_analysis.nosync/track1_v4/submission_ready/` using fresh source staging, adds local explanations of the leading hypotheses, and preserves the prediction CSV byte-for-byte. Edit that folder's `entry_identity.json` and repeat the preparation command to fill in the team, real repository URL and actual AI-use details. Both upload files receive team/model filenames when available. No ranking/comparator rerun, publication, account-setting change or challenge submission occurs. The ZIP is a local convenience bundle; the challenge form expects separate CSV/report/GitHub fields.

For a fresh checkout of the source archive, install `requirements-track1.txt`, obtain authorized inputs and public references separately, and reproduce V3 first. V4 needs its completed local evidence files; those patient-derived artifacts are deliberately absent from the public source archive.

## Validation and references

Run fabricated V4 checks from the project root:

```bash
.venv/bin/python -m unittest discover -s track1_v4/tests -v
```

Public scorer checks skip when its separately obtained snapshot is unavailable. They do not require patient files. See `VALIDATION.md` for the final checks performed on this build.

- [ClinVar classification precedence](https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/)
- [ClinVar review levels](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/)
- [GRCh38 pseudoautosomal regions](https://www.ncbi.nlm.nih.gov/grc/human)
- [BCFtools consequence annotation](https://samtools.github.io/bcftools/howtos/csq-calling.html)
- [MANE transcripts](https://www.ncbi.nlm.nih.gov/refseq/MANE/)
- [Official Track 1 evaluator](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/evaluation.py)

Full privacy history and AI-disclosure limitations remain in `docs/AI_USE_AND_PRIVACY.md`. Public participant methods had been seen previously; this is not described as a blinded discovery.
