"""Local Track 1 files; no patient information is read on import or transmitted."""
from pathlib import Path
import csv
import json
from collections import Counter
import zipfile

from build_track1_v3_submission import csv_row, hypothesis_id, cell, source_archive, archive_previous_exports
from build_track1_submission import ACK
from track1_automated_export import automated_choices
from track1_entry_identity import participant_lines, validate_identity
from validate_track1_submission import FIELDS, validate_csv, validate_package
from track1_baseline import digest, write_json

ROOT = Path(__file__).resolve().parents[1]


def source_bundle(output):
    source_archive(output)
    repo = output / 'repository'
    import shutil
    for path in sorted((ROOT / 'track1_v4').rglob('*')):
        if path.is_file() and path.suffix in {'.py', '.md', '.json'} and '__pycache__' not in path.parts:
            target = repo / path.relative_to(ROOT)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)
    (repo / 'README.md').write_text(
        '# MVA Track 1 V4\n\nSee [the V4 approach and reproduction steps](track1_v4/README.md). '
        'This source-only archive includes V3 dependencies and fabricated examples/tests. '
        'Patient inputs, results, identity details, local logs and public reference databases are excluded.\n')
    with zipfile.ZipFile(output / 'repository.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(repo.rglob('*')):
            if path.is_file() and '__pycache__' not in path.parts:
                archive.write(path, Path('mva-track1-v4') / path.relative_to(repo))


def build_package(output, ranked, evidence, manifest, source_v3_rows, identity, comparator, validation):
    output.mkdir(parents=True, exist_ok=True)
    write_json(output / 'package_status.json', {'completed':False, 'csv_created':False, 'submitted':False})
    archive_previous_exports(output)
    stored_manifest = json.loads((output.parent/'run_manifest.json').read_text())
    if (stored_manifest != manifest or manifest.get('completed') is not True
            or manifest.get('model') != 'mva_track1_v4'):
        raise ValueError('Completed current V4 manifest required')
    hashes = manifest.get('output_hashes', {})
    if not isinstance(hashes, dict) or 'ranked_hypotheses.json' not in hashes:
        raise ValueError('V4 ranked output provenance required')
    for name, expected in hashes.items():
        if not isinstance(name,str) or Path(name).name != name or digest(output.parent/name) != expected:
            raise ValueError('V4 output evidence changed or path invalid')
    if ranked != json.loads((output.parent/'ranked_hypotheses.json').read_text()):
        raise ValueError('Supplied ranking differs from recorded V4 output')
    if any(row.get('automatic_shortlist_eligible') is not True for row in ranked):
        raise ValueError('Ranking includes a deferred V4 hypothesis')
    identity, missing = validate_identity(identity)
    write_json(output / 'entry_identity.json', identity)
    if validation.get('original_genotypes_checked') is not True or validation.get('normalized_reference_checked') is not True:
        raise ValueError('Original calls and reference must be verified before V4 export')
    selected, selection = automated_choices(ranked, digest(output.parent / 'run_manifest.json'), hypothesis_id,
                                            manifest.get('v4_config',{}).get('maximum_submission_hypotheses',10))
    selection['identity'] = identity
    selection['model'] = 'mva_track1_v4'
    selection['case_level_human_review'] = False
    selection['confidence_policy'] = selection['confidence_policy'].replace('v3','V4')
    for choice in selection['choices']:
        choice['reason'] = choice['reason'].replace('v3','V4')
    write_json(output / 'selection.automated.json', selection)
    by_id = evidence if isinstance(evidence, dict) else {v['id']: v for v in evidence}
    with (output / 'predictions.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for row, choice in selected:
            choice = dict(choice, reason='Automated V4 order; EPCR is an ordinal weight, not a causal probability; no human review claimed.')
            writer.writerow(csv_row(row, choice, by_id))
    if not validate_csv(output / 'predictions.csv')['valid']:
        (output / 'predictions.csv').unlink()
        raise ValueError('V4 export failed schema validation')
    unique_v3, seen = {}, set()
    for row in source_v3_rows:
        key = tuple(sorted(row['members']))
        if key not in seen:
            seen.add(key)
            unique_v3[key] = len(seen)
    ties = Counter()
    seen = set()
    for row in ranked:
        key = tuple(sorted(row['members']))
        if key not in seen:
            seen.add(key)
            ties[row['score']] += 1
    comparison = ['# V3 and V4 Track 1 comparison', '',
        'Local patient-derived findings. Ranks are over distinct complete allele sets; repeated disease interpretations are collapsed. '
        'Rank changes do not establish improved accuracy. Upstream exact utility ties are disclosed; unique ordinal EPCR values do not establish confidence differences.', '',
        '| V4 entry rank | Gene | Disease | Alleles | V3 distinct rank | Utility tie group size | Phase |',
        '|---|---|---|---|---|---|---|']
    for n, (row, _) in enumerate(selected, 1):
        comparison.append('| ' + ' | '.join(map(cell, [n, row['gene'], row['disease_id'], ' + '.join(row['members']),
            unique_v3.get(tuple(sorted(row['members'])), 'not in complete V3 ranking'), ties[row['score']], row['phase']])) + ' |')
    (output / 'comparison_v3_v4.md').write_text('\n'.join(comparison)+'\n')
    lines = ['# MVA Track 1 V4 methods', '', *participant_lines(identity), '',
        '## Abstract', '',
        'V4 extends the completed V3 small-variant prioritizer using the same locally stored annotated alleles and public evidence. '
        'It corrects clinical-observation extraction and phenotype-score normalization, preserves X-linked inheritance subtypes, '
        'applies explicit genotype and mechanism compatibility checks, and respects expert-panel or guideline precedence in '
        'condition-specific ClinVar evidence. Numeric ranking weights remain those recorded by V3, permitting an interpretable '
        'comparison. The output is an automated shortlist of complete one- or two-allele hypotheses for Track 1. '
        'There is no model training, independent clinical accuracy estimate, or calibrated causal probability.', '',
        '## Input reuse and provenance', '',
        'The authorized challenge VCF and clinical document stay local. V3 candidate evidence, its manifest and recorded output '
        'hashes are verified before reuse. Current VCF, clinical document, GRCh38 FASTA and used HPO references are checked against '
        'their recorded provenance. Genome annotation, population-database joins and ClinVar indexing are reused; FASTQs were not '
        'realigned. Source code, effective configuration, reference and output hashes are recorded in the V4 manifest.', '',
        '## Track 1 method', '',
        'A structured observation contains the HPO identifier, subject, present/absent/unknown state, source row and any onset text. '
        'Conservative rules exclude ambiguous observations from the positive proband query. A source-bound local JSON override '
        'supports explicit curation. No human review is presumed. Absent observations are retained but are not scored as negative '
        'likelihood evidence in this version. Unmentioned symptoms are not interpreted as absent.', '',
        'Phenotype similarity excludes unsupported information-content terms from both numerator and denominator, reports them '
        'as unscored, and enforces a bounded score. Disease profiles retain recessive and dominant X-linked subtypes. '
        'Genotype, ploidy, pseudoautosomal context and documented disease-specific compatibility determine shortlist eligibility; '
        'unresolved cases remain in the exploratory review output. Sex and de novo status are not inferred.', '',
        'Clinical support remains condition-specific. Expert/guideline evidence takes precedence over ordinary submissions, '
        'while lower-level disagreement is retained. Equal-precedence conflicts remain explicit. Unsupported consequences '
        'against a sourced mechanism are deferred when no matching clinical support resolves the uncertainty. Recessive pair '
        'transcript coherence is checked, and unresolved phase is not called confirmed compound heterozygosity.', '',
        'Default utility retains 60% genetic and 40% phenotype weighting. A pair uses 85% weaker-allele and 15% stronger-allele '
        'genetic utility. Population thresholds, quality criteria and numeric consequence weights are unchanged from the base '
        'configuration; they remain heuristic research choices. Missing population evidence remains unknown. The inherited '
        'MANE/per-record annotation scope does not constitute comprehensive CNV, structural, repeat or deep-intronic analysis.', '',
        '## Automated shortlist and EPCR', '',
        'Up to ten distinct eligible complete allele sets are exported in model order, with both recessive members in one row. '
        'All rows are alternative primary hypotheses; no incidental diagnosis is asserted. EPCR=1/(entry rank+1) encodes ordering '
        'only and is not an estimated or calibrated probability. Underlying utility ties are reported, and tie-breaking does not '
        'imply differential biological confidence. Official pair and F-max semantics are tested with fabricated examples. '
        'No hidden answer or leaderboard feedback is used to choose these predictions.', '',
        '## Results and comparison', '',
        '| Entry rank | Gene and disease | GRCh38 allele hypothesis | Inheritance | Phase |',
        '|---|---|---|---|---|']
    for n, (row, _) in enumerate(selected, 1):
        lines.append('| ' + ' | '.join(map(cell, [n, row['gene']+' / '+row['disease_id'], ' + '.join(row['members']), row['inheritance'], row['phase']])) + ' |')
    lines += ['', 'Comparator status: ' + str(comparator.get('mode', 'unavailable')) + '.', '',
        'Exomiser is a comparison implementation, not statistically independent biological confirmation. Reuse requires matching '
        'VCF provenance and the exact positive HPO payload. A changed payload requires a separate local run or explicit unavailable '
        'status. The V3/V4 and comparator reports compare allele combinations as well as genes. Sensitivity outputs retain full '
        'compact rankings rather than silently treating absence from a truncated table as exclusion.', '',
        '## Validation and limitations', '',
        'Each exported allele was verified against its original called genotype and normalized GRCh38 reference. CSV format and '
        'analysis binding were checked locally. Fabricated regressions cover the corrected failures; they do not estimate clinical '
        'accuracy. Results remain limited by phenotype interpretation, transcript coverage, population callability, phase and '
        'incomplete mechanism evidence. No better leaderboard score is asserted before official evaluation.', '',
        f"The local V4 ranking stage recorded {manifest.get('elapsed_seconds', 0)/60:.2f} minutes, excluding original V3 preparation, any new comparator execution and packaging. No total financial cost is estimated.", '',
        '## AI use and acknowledgement', '',
        'OpenAI Codex assisted with public-code research, implementation and documentation. The participant reports training '
        'opt-out was already enabled; plan/tier and account controls have not been independently verified. The participant-supplied '
        'disclosure must state the actual provider, tier and controls. Earlier conversation contained patient-specific material '
        'and source review exposed case-specific filenames, later removed from shared code. Current V4 data processing and results '
        'stay in local programs and excluded files. No retrospective no-disclosure or no-training guarantee is claimed. '
        'Prior exposure to public participant methods means this work is not described as blinded discovery.', '',
        'Dataset: https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data . Organizer project: '
        'https://www.synapse.org/Synapse:syn76251147/wiki/642892 . Verify the organizer-current formal citation at publication.', '',
        'Methods references: https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/ ; https://www.ncbi.nlm.nih.gov/grc/human ; '
        'https://samtools.github.io/bcftools/howtos/csq-calling.html ; https://www.ncbi.nlm.nih.gov/refseq/MANE/ ; '
        'https://github.com/exomiser/Exomiser .', '', ACK, '']
    (output / 'methods_submission.md').write_text('\n'.join(lines))
    write_json(output / 'readiness.json', validate_package(output/'predictions.csv', output/'methods_submission.md', identity.get('github_url') or None, True))
    (output/'START_HERE.md').write_text(
        '# Track 1 V4 local package\n\nOpen predictions.csv and methods_submission.md locally. '
        'comparison_v3_v4.md explains rank changes; ../phenotypes.json records input interpretation; '
        '../review_hypotheses.json preserves deferred hypotheses. ../comparator_comparison.md gives the separate comparison.\n\n'
        'Complete entry_identity.json with team, GitHub URL and accurate AI disclosure, then run the V4 wrapper again. '
        'Publish only repository/ source after account/destination details are supplied. Follow the official Submit Track 1 form '
        'for separate CSV, Markdown report and GitHub URL; rename a copy of the CSV with your team and model. '
        'submission_bundle.zip is a local convenience bundle. No upload or challenge submission is performed. '
        'Patient findings and private logs must not be pasted into an AI assistant.\n')
    source_bundle(output)
    with zipfile.ZipFile(output/'submission_bundle.zip', 'w', zipfile.ZIP_DEFLATED) as bundle:
        for name in ('predictions.csv', 'methods_submission.md', 'repository.zip'):
            bundle.write(output/name, name)
    result = {'completed': True, 'csv_created': True, 'submitted': False,
        'source_archive_sha256': digest(output/'repository.zip'), 'participant_fields_missing': missing,
        'clinical_accuracy_established': False, 'prediction_csv_sha256': digest(output/'predictions.csv')}
    write_json(output/'package_status.json', result)
    return result
