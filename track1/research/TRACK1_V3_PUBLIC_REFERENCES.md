# Track 1 v3 public reference preparation

Prepared 7 October 2026; documentation updated 8 October 2026. All preparation requests use fixed complete-release URLs. The preparation does not open challenge data or patient-derived results. A separate offline function accepts locally supplied allele keys and writes its results only to the caller's local directory.

## Coverage and storage choice

The previous 11-gene public catalogue is useful for detailed source checks but cannot give equal population and missense coverage to genes outside that panel. Version 3 therefore reuses Exomiser's complete production GRCh38 variant MVStore for population and pretrained predictor annotations. This avoids duplicating its large dataset. The comparator and its public reference downloads are managed separately; the adapter alone does not establish that the production database has finished downloading or that Exomiser has run successfully.

Current Exomiser releases support population annotations and multiple effect predictors. Its stored allele schema contains repeated frequency and pathogenicity records, with source enums including AlphaMissense and SpliceAI. Actual availability is checked per allele; supported sources are not guaranteed to be present for every record. [Exomiser source schema](https://github.com/exomiser/Exomiser/blob/15.1.0/exomiser-core/src/main/proto/allele.proto), [installation requirements](https://exomiser.readthedocs.io/en/latest/installation.html).

The complete public AlphaMissense hg38 archive was size-checked at 642,961,469 compressed bytes. It was **not downloaded** for this build because the chosen production database already provides genome-wide predictor coverage. Similarly, the 6.295 GB compressed ClinVar RCV XML was inspected only to understand its schema; it was not downloaded in full. Preparation code retains optional import paths for these resources, but their presence in code is not a claim that they were used.

Instead, the public ClinVar variant summary and individual submission summary are downloaded in full. Together they occupy 858,560,130 compressed bytes. These support a transparent condition-specific **SCV** join without claiming that individual submissions are aggregate RCV assertions.

## ClinVar disease-specific evidence

The variant summary maps ClinVar Variation IDs to GRCh38 alleles. The submission summary attaches each submitted interpretation, its condition, review status, submitter and date to a Variation ID and SCV accession. The parser uses only the condition identifiers on that individual SCV; it never assigns all diseases associated with a variant to every submission. Somatic classification columns are not imported as germline evidence. [ClinVar tabular schema](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/README).

Alleles are restricted to ordinary sequence representations on primary GRCh38 chromosomes. REF is checked against the local public hg38 FASTA, common padding is trimmed and ordinary indels are left aligned. Variations with multiple source AlleleIDs are excluded before assembly or representation filtering, including composites whose other component has no supported GRCh38 representation. Variations mapping to multiple normalized primary alleles, or with an unresolved primary GRCh38 mapping, are also excluded. These exclusions and REF mismatches are counted; they must not be treated as evidence of absence.

Each retained assertion stores its SCV accession, exact submitted-condition IDs, classification, review status/stars, evaluation date, submitter, collection method, allele origin, Variation ID and whether it contributes to the current aggregate classification. An explicit noncontributing flag should exclude the assertion from active classification evidence while leaving it available for provenance. Differences across distinct conditions are not automatically conflicts for the disease being tested.

SCV index version 4 preserves distinct condition or provenance payloads when the same SCV accession occurs more than once for a normalized allele. Only rows with an identical normalized allele, accession and serialized assertion payload are deduplicated. Thus an accession alone does not identify a unique retained assertion. Metadata distinguishes input assertions accepted for indexing (`submission_assertions_indexed`), distinct records retained (`distinct_assertion_records_retained`) and identical assertion rows deduplicated (`identical_assertion_rows_deduplicated`). The tab-separated parser preserves literal quotation marks in public free text, permits long fields and rejects inconsistent column counts; CSV quotation rules must not merge separate ClinVar records.

Many submitted conditions use MedGen identifiers, while phenotype profiles use OMIM. A separate complete public NCBI `mim2gene_medgen` snapshot provides aliases. Only rows labeled **phenotype** are accepted; gene entries are excluded. This establishes disease-ID correspondence independently of the patient's variants and independently of a ClinVar variant's list of conditions. Unmapped conditions remain unmatched. [NCBI mapping](https://ftp.ncbi.nlm.nih.gov/gene/DATA/mim2gene_medgen).

## SQLite interface

The local file is `references/track1_v3/public_annotations.sqlite`. Consumers must verify source completeness in `metadata` before analysis.

| Object | Purpose |
|---|---|
| `cv_variant` | Normalized primary-GRCh38 allele mappings by ClinVar Variation ID |
| `scv_records` | Distinct submitted assertion payloads, keyed by normalized allele, SCV accession and payload |
| `clinvar_scv` view | Allele-keyed JSON arrays of submitted assertions |
| `clinvar_rcv` view | Compatibility alias of `clinvar_scv`; payloads explicitly say `assertion_type: SCV` |
| `metadata` | Completeness, source hashes, counts, normalization and condition aliases |
| `population` | Empty compatibility table; population evidence is supplied by the separate Exomiser adapter |

The source view contract is `chrom, pos, ref, alt, payload`, where `chrom` omits `chr`, mitochondrial `MT` maps to `M`, and coordinates are one-based. The compatibility view name must never be presented in reports as proof that RCV assertions were loaded.

## Offline Exomiser adapter

`scripts/prepare_track1_v3_references.py` exposes:

```python
annotate_keys(keys, mvstore_path, output_dir) -> dict
```

`keys` contains `(chrom, pos, ref, alt)` tuples. The result is keyed by `chrom:pos:ref:alt`, with `population`, `other_predictors`, `database_match` and `alphamissense_annotations` fields. The wrapper opens the MVStore read-only. It performs no downloads or network requests. Inputs, outputs and diagnostic logs stay in `output_dir`; no allele values are printed by the function.

The public-only preparation step compiles a small Java adapter against Exomiser 15.1.0 libraries using Eclipse ECJ 3.39.0 and Java 21. It records the wrapper source, class and compiler hashes. Each offline call verifies the cached class receipt against the current embedded Java source. Calling `annotate_keys` before compilation, or with a stale cached class, fails rather than fetching dependencies during private analysis.

Frequency units require care: when counts exist, the adapter uses `AC / AN`; older frequency-only records store percentages, so those values are divided by 100. An absent entry is unknown, not zero frequency. Exomiser's protobuf does not retain per-allele gnomAD FILTER or callable-region flags; these are returned as `null`, with a stated limitation. Population groups are kept separate rather than averaged or added. [Frequency representation](https://github.com/exomiser/Exomiser/blob/15.1.0/exomiser-core/src/main/java/org/monarchinitiative/exomiser/core/model/frequency/Frequency.java), [official conversion code](https://github.com/exomiser/Exomiser/blob/15.1.0/exomiser-core/src/main/java/org/monarchinitiative/exomiser/core/model/AlleleProtoAdaptor.java).

The adapter includes `af_basis` (`AC/AN`, `stored_percentage/100`, or `unknown_or_invalid`). Without positive AN, a stored frequency of zero cannot be distinguished from an unset protobuf field, so it returns `af: null` and `observed: null`. Likewise, a zero homozygote field is emitted as `homozygotes: null` with `homozygotes_zero_or_missing: true` when counts exist; positive homozygote counts are retained. A known AN with AC zero remains an explicit zero frequency, while callability and per-allele quality remain unknown.

Stored predictor scores do not include a transcript identifier or protein substitution. Therefore AlphaMissense and SpliceAI from this source are emitted under `other_predictors` with `transcript_unresolved: true`; the adapter emits no invented transcript-specific AlphaMissense annotation. Existing independently checked transcript-specific annotations may be retained separately. These predictors provide computational evidence, not functional assays or patient-level diagnostic probabilities.

## Reproduction and checks

```bash
# Public downloads and ClinVar indexing only
.venv/bin/python scripts/prepare_track1_v3_references.py

# Verify cached source hashes and build/reuse the index without network access
.venv/bin/python scripts/prepare_track1_v3_references.py --offline

# One-time public adapter compilation
.venv/bin/python scripts/prepare_track1_v3_references.py --prepare-mvstore-adapter

# Synthetic tests, never challenge data
.venv/bin/python -m unittest discover -s tests -p test_track1_v3_references.py -v
```

Fourteen synthetic tests passed, including an actual Java MVStore round trip using the cached adapter. They check count-based versus percentage frequencies, protobuf missingness, mitochondrial chromosome conversion, padding/coordinate shifts and haplotype-preserving repeat normalization, REF mismatches, source receipt validation, ClinVar comment/header handling, long public fields, literal quotes, malformed rows, condition isolation, composite exclusion and phenotype-only aliases. The repeated-SCV test verifies that two distinct condition/classification payloads survive under the same allele and accession while an identical duplicate is removed. A fresh-process regression also imports the offline adapter after the socket guard is installed and verifies that networking remains disabled. These validate software behavior, not diagnostic performance.

Downloads have URL, retrieval timestamp, source byte count and SHA-256 receipts. Cached and offline sources must match all of those receipt fields. NCBI compressed sources are checked against published MD5 values on acquisition. ClinVar index metadata records both tabular source hashes, the actual uncompressed reference FASTA hash, source URLs and index version; reuse requires each to match, including `SCV_INDEX_VERSION = 4`. Preparation enforces at least 6 GiB free space. Partial downloads and incomplete SQLite sources remain distinguishable from completed references. Final run counts and completion status are stored in SQLite metadata rather than inferred from file existence.
