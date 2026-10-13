# Track 1 v3: biology, evidence and the complete local workflow

This is a research workflow for ranking possible genetic explanations of a child's condition. It does not train a model on the child, establish a clinical diagnosis or use the hidden challenge answer. The newest code is `scripts/track1_v3.py`. Actual run completion, reference coverage and comparator status are recorded in the private output manifests; installing a capability is not evidence that it ran successfully.

## The biological question in plain language

DNA contains instructions that cells use to make proteins. A gene is a stretch of those instructions. Most genes are present in two copies, one inherited from each parent. A variant is a difference from a reference DNA sequence. Everyone has many variants; most do not cause disease.

The challenge is therefore not to find a DNA difference. It is to identify a difference, or a combination of differences, that can explain the clinical features through a credible biological mechanism.

In mosaic variegated aneuploidy, different cells can have abnormal numbers of chromosomes. Chromosomes must be copied and distributed accurately when a cell divides. The spindle checkpoint is part of the machinery that delays division until chromosomes are appropriately attached. Other proteins help organize the spindle, maintain chromosome structure or repair DNA. Disruption of these processes can lead to chromosome instability. Different genes and different changes within a gene can have different effects; membership in a chromosome-instability gene list is not proof of pathogenicity.

“Mosaic” describes differences between populations of cells. The proportion of sequencing reads carrying one small variant does not, on its own, diagnose chromosome mosaicism or measure chromosome numbers.

Three distinctions matter throughout the analysis:

| Question | What it means | Evidence that helps |
|---|---|---|
| Is the DNA call real? | The reported change is supported by the sequencing data. | Call quality, depth, allele balance, reference matching and, where available, read review. |
| Can it disrupt the relevant biological function? | The change affects the right transcript/protein in a way known to cause the disease. | Variant consequence, condition-specific clinical assertions, functional studies and computational predictions. |
| Does it explain this case? | The genotype, inheritance and clinical features fit the proposed condition. | Disease-specific phenotype matching, both alleles for recessive disease, phase, segregation and clinical review. |

A strong answer to one question does not automatically answer the others.

## Why transcripts and inheritance matter

Cells copy DNA into RNA and remove intervening segments called introns. The remaining exons form a transcript, which may be translated into a protein. A gene can have several transcripts. The same DNA change can have different consequences in different transcripts, so v3 evaluates transcript-specific evidence before selecting the annotation used for a hypothesis.

A missense change replaces one amino acid, a building block of a protein. It may be harmless or may impair a particular function. A stop-gained change introduces an early stop; a frameshift changes the reading frame; a splice change can alter how RNA is assembled. These categories are useful clues, but an early stop is not automatically a disease-causing loss of function. Relevant exon usage, the affected protein region, and nonsense-mediated decay—the cell's removal of some abnormal RNAs—matter. The workflow flags these review needs rather than assigning an automatic clinical classification. [ClinGen loss-of-function recommendations](https://pubmed.ncbi.nlm.nih.gov/30192042/)

For an autosomal recessive condition, the usual hypothesis requires both copies of the gene to be affected. A homozygous call carries the same alternate allele on both copies. A compound-heterozygous hypothesis involves two different variants. Those variants must normally be in trans, on opposite copies; a cis pair affects the same copy. Unresolved phase remains unresolved. Existing genotype phase fields can provide evidence, but software must not invent parental origin or phase.

For a dominant condition, one changed copy can be sufficient, but the relevant mechanism still matters. Reduced function, a new or increased activity, and interference with the normal protein are different mechanisms. X-linked conditions also require attention to chromosome context, ploidy, inheritance subtype and penetrance. A novel dominant variant need not already appear in ClinVar to deserve review. Absence of parental data does not establish that it arose de novo.

## What the input data can and cannot tell us

The supplied VCF contains small-variant calls and sample-level evidence such as genotype and call quality. FASTQ files contain sequencing reads. The main v3 workflow starts from the supplied VCF; it does not silently claim that the reads were realigned or that a new variant caller was run.

The clinical document is converted locally into Human Phenotype Ontology identifiers. The ontology lets the model recognize relationships between specific and more general clinical descriptions. Proband and family observations remain distinct. A feature that was not mentioned is not treated as an observed absence. Technical sequencing lanes are not independent patients or valid train/test splits.

No structural-variant, copy-number, repeat-expansion or genome-wide regulatory analysis is established by this small-variant pipeline. The separate ClinVar screen and Exomiser comparator are useful additional views, but their actual filtering and completion must be inspected in their local records. Missing or uncalled variants cannot be recovered simply by changing ranking weights.

## What changed from v2

| Area | v2 limitation | v3 behavior |
|---|---|---|
| Population evidence | Modern reference catalogues covered eleven genes. | An offline adapter supports the genome-wide public Exomiser variant database; retained frequency evidence and missingness are recorded explicitly. The run manifest states whether the database was actually used. |
| Clinical assertions | A disease identifier could be combined with an aggregate ClinVar classification spanning several conditions. | Only a matching condition-specific SCV or RCV assertion supplies condition-specific support. Conflicting assertions are retained and flagged. |
| Allele representation | Minimal trimming could miss equivalent indels in repeats. | Reference-checked normalization left-aligns sequence indels and retains original source alleles for genotype verification. |
| Transcripts | The annotation with the largest consequence weight was chosen before all evidence was considered. | Each available transcript is evaluated against the disease and compatible predictor evidence. |
| Inheritance | Existing ClinVar support was required for dominant candidates; heterozygous X-linked hypotheses were excluded. | Novel dominant and X-linked hypotheses can enter review, with explicit mechanism and context requirements. Incomplete recessive and unknown-inheritance findings are kept in a separate review output. |
| Phase | Every unknown pair received one fixed penalty. | Phase is an explicit evidence state; an alternative penalty is tested in sensitivity analysis rather than described as a probability of trans phase. |
| Submission | The first hypothesis was exported with provisional EPCR 0.60. | The automated export preserves the top ten distinct complete allele hypotheses and uses disclosed ordinal EPCR weights, not causal probabilities. An optional human-reviewed export remains available. Both routes are bound to the analysis hash and validated. |
| Comparison | No established independent prioritizer had been run. | A separate local Exomiser workflow and public demonstration setup are provided. Completion and limitations are recorded, not assumed. |

## How the ranking works

1. Verify the input/reference provenance and annotate primary human chromosome calls with BCFtools csq and the MANE transcript reference. Reuse an annotation only when its input hashes match.
2. Retain observed sequence alleles, check quality, and normalize them against GRCh38. Preserve the original allele and genotype evidence.
3. Join public annotations locally. Population data, condition-specific clinical assertions and computational predictions are separate evidence types. Missing annotation never means “known to be rare.”
4. Match clinical features to individual gene–disease profiles rather than pooling every disease associated with a gene. Construct inheritance-compatible single-allele or two-allele hypotheses. Reject established cis recessive pairs.
5. Evaluate transcript and mechanism evidence. Rank complete hypotheses with a transparent utility score; keep incomplete genotype and unresolved-inheritance findings in a separate review file.
6. Compare alternative scoring choices, check against Exomiser when its run completes, and export a local ranked submission draft. The automated export preserves the model order; an optional human-reviewed alternative can supply independently justified confidence.

The default score is an explicitly uncalibrated weighted combination of genetic evidence and disease-specific phenotype similarity. Exact weights and effect-category values are pinned in `config/track1_v3.json` and copied into every run manifest. For a two-allele hypothesis, the weaker allele receives 85% of the genetic component and the stronger receives 15%. This makes a poorly supported second allele difficult to hide behind one strong allele. The arithmetic is a research design choice, not a biological law or a learned diagnostic probability.

Condition-specific clinical evidence can strengthen allele support, but a pathogenic variant for one disease need not explain this child's symptoms. Review tiers describe the strength and gaps in evidence; they do not turn a clinical assertion into proof of a case-level explanation.

The default common-evidence ranking does not simply add every predictor score. AlphaMissense, REVEL and other predictors can share biological signals and are not independent votes. Genomic database scores may also lack the transcript context needed to transfer a prediction safely. Compatible predictor evidence is reported, and the effect of adding AlphaMissense support is examined separately. SpliceAI evidence, when available, is computational support and does not establish the RNA outcome. [AlphaMissense study](https://doi.org/10.1126/science.adg7492)

Population filtering uses observed public frequencies with explicit source and availability information. A browser catalogue can contain sample counts and quality fields that a compact compiled database omits. Those missing quality and callability details must remain visible; a database miss cannot be interpreted as a measured frequency of zero. Generic frequency thresholds are screening rules, not gene-specific maximum credible frequencies or penetrance estimates.

The phenotype model is semantic similarity, not a full clinical likelihood model. It does not comprehensively model age of onset, feature frequency, penetrance or explicitly absent features. A likelihood-based method such as LIRICAL is a relevant future comparator; it is not claimed to have run here. [LIRICAL primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC7477017/)

## What validation means here

There are several different kinds of validation, and they answer different questions:

- Fabricated regression tests check normalization, assertion matching, inheritance, missing evidence, phase handling and safe export. They do not estimate diagnostic accuracy.
- Original-call and reference checks confirm that review alleles occur in the supplied called genotype and match GRCh38. They do not prove that the call is biologically causal.
- Sensitivity analyses show whether a rank depends heavily on one scoring choice. Stability is not correctness.
- Exomiser is an independently implemented comparator. Agreement is useful, but the tools can share ClinVar, population and phenotype references, so agreement is not independent experimental evidence. [Exomiser evaluation](https://www.nature.com/articles/s41525-024-00456-2)
- A released toy demonstration verifies installation and expected behavior on that example. It is not a held-out clinical benchmark, and no percentage accuracy should be inferred from it.

True performance assessment needs separate, appropriately licensed cases with known diagnoses, a fixed evaluation protocol, and separate tuning and evaluation sets. It should report top-1/top-5/top-10 recovery, recovery of both recessive alleles, annotation coverage and failure categories. The current challenge case and its technical lanes cannot supply that cohort. No leaderboard feedback or hidden clinical answer is used to set weights.

## Local execution and review

Public reference preparation and private analysis are separate operations. The reference preparers download fixed public releases, record hashes and never build an internet request from a patient's genotype or phenotype. Large reference downloads are resumable and must preserve disk headroom. Reference-specific commands and provenance are in `research/TRACK1_V3_PUBLIC_REFERENCES.md` and comparator setup is in `research/TRACK1_V3_COMPARATOR.md`.

Once public references are ready:

```bash
.venv/bin/python scripts/run_track1_v3_local.py
```

The wrapper captures detailed execution output in local logs and returns only generic status. New private files use `work/private_analysis.nosync/`, whose exclusion from macOS cloud sync was verified on this computer. The analysis is in its `track1_v3/` subfolder, the comparator in `track1_exomiser/`, and the review package in `submission_track1_v3/`. If you use a different output location, verify its sync exclusion first; the folder name alone is not a guarantee on every service or computer.

The automated `predictions.csv` contains up to ten distinct complete hypotheses in the existing model order. Its `epcr` field is `1/(submission_rank+1)`: an ordinal encoding used solely to preserve strict ordering in the challenge scorer. These numbers are not estimated or calibrated causal probabilities. The formal field name is a probability label, so the methods explicitly disclose this limitation. The organizer states that relative order determines scoring and its methods template prefers automated output; human review is not a challenge requirement. The current public scorer gives identical rank points and F-max for strictly monotonic encodings of the same ranking. The absolute threshold value changes and has no clinical interpretation.

`methods_submission.md` is the concise report; `submission_bundle.zip` contains that report, the CSV and the source archive. Fill `entry_identity.json` with team, GitHub and actual AI-use details, then rebuild the export to apply them. The existing analysis can be packaged without rerunning the genome using `scripts/build_track1_v3_submission.py --automated`. Detailed local findings remain in `methods.md` and `comparison.md`.

Open `work/private_analysis.nosync/submission_track1_v3/index.html` locally for optional human review. It makes no network requests. Review the original calls, both members of recessive pairs, disease mechanism, assertions, population evidence, transcript context and phase. The page can save a reviewed selection with your own confidence and rationale; save it in the excluded folder. No automatic export claims that this review happened.

For an independent export check, supply the saved manifest to the Python finalizer:

```bash
.venv/bin/python scripts/build_track1_v3_submission.py \
  --selection /path/to/selection.reviewed.json \
  > work/private_analysis.nosync/track1_v3/finalize.log 2>&1
```

The finalizer rejects stale evidence, incomplete genotype hypotheses and malformed CSV records. In reviewed mode it also requires the supplied review and confidence rationale. Team identity, a real GitHub repository URL, accurate AI-use disclosure and a publication decision remain participant-controlled. Follow `HOW_TO_SUBMIT.md` for the actual entry steps. A local file does not consume a challenge submission slot. No upload is performed by these scripts.

The source archive includes code, configurations and fabricated tests. It excludes patient inputs, detailed findings, private logs and downloaded reference databases. The ranked review table is not an official submission file. A methods draft containing unresolved participant fields is not yet ready for publication.

## Privacy and interpretation

Patient data, local findings and detailed logs must not be returned to a hosted assistant, uploaded for annotation or placed in a public repository. Local reports are for the participant's review. The AI-use disclosure must accurately describe the product/tier and actual controls, the participant's reported training opt-out, and any earlier conversation exposure. The pipeline does not establish or alter account-level data settings. See `docs/AI_USE_AND_PRIVACY.md`.

This approach is designed to make evidence and uncertainty easier to assess. It does not promise a higher leaderboard score or claim a clinically confirmed result. Additional read review, family sequencing, RNA experiments or functional assays can resolve uncertainties that another numerical model cannot.

## Primary sources and technical references

- [Official challenge](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026).
- [BCFtools consequence calling](https://samtools.github.io/bcftools/howtos/csq-calling.html). The workflow uses local per-record effects; it does not claim haplotype-aware frame rescue or compound-codon annotation.
- [MANE reference transcripts](https://www.ncbi.nlm.nih.gov/refseq/MANE/).
- [ClinVar variant–condition records](https://www.ncbi.nlm.nih.gov/clinvar/docs/details/).
- [ClinVar review status](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/).
- [Exomiser source and documentation](https://github.com/exomiser/Exomiser).
- [MAD2L1BP chromosome-instability study](https://insight.jci.org/articles/view/170079).
- [CEP192 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10716027/).
- [SLF2 and SMC5 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9636423/).

Public participant methods had been consulted during earlier project work. This project is therefore not described as a blinded discovery. Published results from other teams are not our validation results.
