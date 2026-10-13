# Track 1 v3 independent local comparator

Exomiser 15.1.0 is installed with portable Java 21 under `references/exomiser`. It is a separate implementation from the project ranker. Because both pipelines use some of the same public annotations and the case has already been studied, this is not an independent clinical validation or a blinded analysis.

The production local comparison completed successfully on October 9. Its reference integrity and unchanged-input checks passed, and the resulting review package accepted it as matching the current main analysis. This operational status discloses no patient finding and does not establish clinical accuracy or a leaderboard score. Detailed comparison results remain in the excluded local review package.

## Reference data and integrity

The production comparator uses Exomiser's 2602 GRCh38 and phenotype releases. The official [2602 release announcement](https://github.com/exomiser/Exomiser/discussions/640) documents the source versions and links to both fixed archives. Downloads use complete public archives and byte ranges selected only by ZIP structure; no patient input controls any network request. Optional `phenix/` assets are omitted because the selected preset uses hiPhive.

The wrapper downloads one compressed member at a time, retains hashed chunk receipts for interrupted transfers, verifies ZIP member size and CRC after decompression, and records SHA-256 hashes. It enforces a 6 GiB minimum disk reserve. Before a local analysis, it verifies that every required database and Ensembl transcript file is listed in the manifests, then checks the exact size and SHA-256 of all retained members. Incomplete `.part` files cannot satisfy this check. These locally computed hashes establish consistency with the acquired archives, not independent publisher-signed authentication.

Public receipts are stored in `references/exomiser/data/2602_hg38.zip.manifest.json`, `2602_phenotype.zip.manifest.json`, and `production_verification.json`. The shared variant database is `references/exomiser/data/2602_hg38/2602_hg38_variants.mv.db`.

### Recovering an evicted public variant database

macOS File Provider evicted the large public variant database after its original acquisition. Repeated cloud hydration attempts timed out. Its generated, zero-allocation cloud placeholder was removed only after checking its dataless flag and exact size against the acquisition manifest; `public_placeholder_cleanup.json` records that operation. The original acquisition manifest remains unchanged. A previous verification receipt does not establish current availability.

`scripts/restore_exomiser_public_cache.py` restores this one fixed member directly from the official archive. It checks absolute HTTP byte ranges, keeps resumable 4 MiB compressed chunks with SHA-256 receipts, and deletes each consumed chunk only after decompressed output is durably written. Final size, ZIP CRC and the original acquisition SHA-256 must all match before the script publishes an internal relative symlink. The target is `references/exomiser/data/local_cache.nosync/`, verified as excluded from cloud sync on this Mac. It refuses to replace an existing materialized canonical database.

Allow approximately 36 GiB of free space for this restoration, including the enforced 6 GiB reserve, and more for subsequent analysis. Download interruption preserves verified chunks; inflation interruption restarts decompression from the beginning and redownloads consumed chunks. The script does not read patient files. A server response with an incorrect byte range is rejected rather than accepted as a cache block.

```bash
.venv/bin/python scripts/restore_exomiser_public_cache.py --workers 32
# After a successful restore, run the analysis, comparator and review packaging.
.venv/bin/python scripts/run_track1_v3_local.py
```

## Analysis scope

The wrapper uses the vendor's `preset-exome-analysis.yml` and evaluates all input variants before its vendor filters; there is no fixed candidate-gene whitelist. The preset removes many noncoding variant classes, so this is a coding-focused analysis of a whole-genome input rather than comprehensive noncoding genome prioritization. The comparator uses Ensembl transcripts, supplied positive proband HPO terms, the selected single sample, and the vendor inheritance modes. ClinVar variant whitelisting is explicitly disabled.

The export groups vendor-designated contributing PASS variants by gene and inheritance mode. It preserves compound pairs as one hypothesis, records groups with more than two alleles as incompatible with a one-or-two-variant challenge row, and never treats a comparator pair as confirmed phase. Allele joins use a `chr` prefix and minimal trimming, not repeat-aware FASTA left alignment. Exomiser scores are not calibrated diagnostic probabilities.

## Privacy and output

Patient-derived analysis YAML, JVM logs, native-library output, result tables and hypothesis JSON are written under `work/private_analysis.nosync/track1_exomiser`. The private parent folder was verified as excluded from macOS cloud sync on this computer. The wrapper redirects Python, native stdout/stderr and subprocess output to private files. Its terminal response is only a generic completion or failure message. The main local artifacts are `comparator_status.json` and `comparator_hypotheses.json`; they must not be printed into model context or transmitted to remote services. A custom output directory must have its sync exclusion checked separately.

The wrapper never submits challenge predictions. Completion status records this explicitly and also records that official challenge scoring and clinical accuracy have not been established.

Each comparator rerun first marks its status incomplete and writes raw results to a new run directory. It records hashes of the VCF, phenotype input, HPO payload, preset, generated configuration, wrapper, extractor, Java runtime, Exomiser libraries, verified reference manifests and result exports. Input hashes are checked again after execution. The review package excludes a comparator whose VCF or phenotype differs from the current analysis, whose verification is incomplete, or whose exported hypotheses changed after the run.

## Public synthetic demonstration

The official supplied Pfeiffer example ran successfully with the vendor test hg19 database and human-only phenotype preset. The expected public example gene, FGFR2, ranked first. The test genomic database includes only demonstration variants and the phenotype data are restricted; one synthetic positive example does not measure clinical sensitivity, specificity or generalization. This is a software integration demonstration, not an accuracy benchmark.

Public reproducible outputs are under `research/track1_exomiser_benchmark`, including `benchmark_status.json`, analysis configuration and logs. Wrapper regression tests use synthetic data and cover resume integrity, local-only inputs, path traversal, reference size/hash failure, OS-descriptor privacy capture, grouped hypotheses, phase caveats and disabled ClinVar whitelisting.

After the successful demonstration, the two large public test phenotype database files were removed to preserve disk space for the production release. The benchmark report and source manifests remain. Run `--prepare-benchmark-data` again before repeating the demonstration; a cleanup receipt records the removed files and their hashes.

## Reproduction

```bash
# Public downloads only; no patient file reads.
.venv/bin/python scripts/run_track1_exomiser.py --prepare-data
.venv/bin/python scripts/run_track1_exomiser.py --verify-data

# Public vendor demonstration only.
.venv/bin/python scripts/run_track1_exomiser.py --prepare-benchmark-data
.venv/bin/python scripts/run_track1_exomiser.py --benchmark

# Fully local comparator; detailed output stays in private local files.
.venv/bin/python scripts/run_track1_exomiser.py

# Synthetic regression tests only.
.venv/bin/python -m unittest discover -s tests -p test_track1_exomiser.py -v
```
