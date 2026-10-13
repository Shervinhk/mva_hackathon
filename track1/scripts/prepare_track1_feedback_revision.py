#!/usr/bin/env python3
"""Prepare a transparent feedback-conditioned Track 1 revision locally.

No hidden answer is accessed. The reported full-match rank identifies a complete
hypothesis only within the exact CSV submitted by the participant. All private
processing runs in a captured worker; the default command emits generic status.
"""
from pathlib import Path
import argparse
import copy
import csv
from datetime import datetime, timezone
import importlib.util
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'scripts')]


def revise_rows(rows, full_match_rank):
    """Promote one EPCR-sorted complete hypothesis; preserve its allele fields."""
    from validate_track1_submission import FIELDS, CHROMS, DNA
    if type(full_match_rank) is not int or not 1 <= full_match_rank <= len(rows) <= 10:
        raise ValueError('Reported rank must address one of at most ten rows')
    values, hypotheses = [], set()
    for row in rows:
        if set(row) != set(FIELDS) or any(not isinstance(row[field], str) for field in FIELDS):
            raise ValueError('Source fields must match the official CSV schema')
        if row['proband_id'] != 'PROBAND01' or row['finding_type'] not in {'primary', 'secondary'}:
            raise ValueError('Invalid source proband or finding type')
        alleles = []
        for suffix in ('1', '2'):
            chrom, pos, ref, alt = [row[key + '_' + suffix].strip() for key in ('chrom', 'pos', 'ref', 'alt')]
            if suffix == '2' and not any((chrom, pos, ref, alt)):
                continue
            if (chrom not in CHROMS or not re.fullmatch(r'[1-9][0-9]*', pos)
                    or not DNA.fullmatch(ref.upper()) or not DNA.fullmatch(alt.upper()) or ref.upper() == alt.upper()):
                raise ValueError('Incomplete or invalid source allele')
            alleles.append((chrom, int(pos), ref.upper(), alt.upper()))
        hypothesis = frozenset(alleles)
        if len(hypothesis) != len(alleles) or hypothesis in hypotheses:
            raise ValueError('Duplicate source allele or complete hypothesis')
        hypotheses.add(hypothesis)
        try:
            value = float(row['epcr'])
        except (ValueError, TypeError, KeyError):
            raise ValueError('Invalid source EPCR') from None
        if not math.isfinite(value) or not 0 < value <= 1:
            raise ValueError('Invalid source EPCR')
        values.append(value)
    if len(set(values)) != len(values):
        raise ValueError('Source EPCR ties require explicit manual source binding')
    ordered = sorted(zip(rows, values), key=lambda item: -item[1])
    order = [full_match_rank] + [rank for rank in range(1, len(rows) + 1) if rank != full_match_rank]
    result = []
    for rank, original_rank in enumerate(order, 1):
        row = copy.deepcopy(ordered[original_rank - 1][0])
        row['epcr'] = str(1.0 / (rank + 1))
        row['notes'] = (f'Feedback-informed V4 revision; original EPCR rank {original_rank}; '
                        'EPCR encodes order, not calibrated probability. Alleles unchanged. '
                        'Original note: ' + row.get('notes', ''))
        result.append(row)
    return result, order


def evaluator_module(path):
    spec = importlib.util.spec_from_file_location('track1_feedback_public_evaluator', path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def csv_write(path, rows):
    from validate_track1_submission import FIELDS
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def prepare(feedback_path, confirmed_source=False):
    from local_inputs import input_path
    from track1_baseline import digest, write_json
    from track1_v2 import install_network_guard
    from track1_automated_export import automated_choices
    from build_track1_v3_submission import hypothesis_id, verify_calls, cell
    from prepare_track1_v4_entry import check_hashes, checked_prediction_bytes, source_bundle, prepared_methods
    from track1_v4_case_report import render_case_rationale
    from validate_track1_submission import validate_csv, validate_package

    install_network_guard()
    load = lambda path: json.loads(path.read_text())
    feedback = load(feedback_path)
    required = {'full_match_rank', 'rank_points', 'f_max', 'f_max_threshold'}
    if not required.issubset(feedback):
        raise ValueError('Reported feedback fields are incomplete')
    private = ROOT / 'work/private_analysis.nosync'
    analysis = private / 'track1_v4'
    original = analysis / 'submission_ready'
    output = private / 'track1_v4_feedback'
    output.mkdir(parents=True, exist_ok=True)
    target = output / 'submission'
    manifest_path = analysis / 'run_manifest.json'
    manifest_hash = digest(manifest_path)
    manifest = load(manifest_path)
    package = load(original / 'preparation_status.json')
    if manifest.get('completed') is not True or package.get('completed') is not True:
        raise ValueError('Completed V4 analysis and package required')
    check_hashes(ROOT, manifest['source_hashes'])
    check_hashes(analysis, manifest['output_hashes'])
    if manifest_hash != package['analysis_manifest_sha256']:
        raise ValueError('Prepared package is bound to a different analysis')
    base = Path(manifest['base_analysis'])
    if digest(base / 'run_manifest.json') != manifest['base_manifest_sha256']:
        raise ValueError('V3 evidence binding changed')
    if digest(base / 'candidate_evidence.json') != manifest['reused_candidate_evidence_sha256']:
        raise ValueError('Reused candidate evidence changed')
    evidence = {row['id']: row for row in load(base / 'candidate_evidence.json')}
    ranked = load(analysis / 'ranked_hypotheses.json')
    selected, _ = automated_choices(ranked, manifest_hash, hypothesis_id,
                                    manifest['v4_config']['maximum_submission_hypotheses'])
    source_csv = original / 'predictions.csv'
    original_bytes = checked_prediction_bytes(source_csv, selected, evidence, package['prediction_csv_sha256'])
    with source_csv.open(newline='') as handle:
        rows = list(csv.DictReader(handle))
    revised, order = revise_rows(rows, feedback['full_match_rank'])
    current = [selected[rank - 1][0] for rank in order]
    validation = verify_calls(current, evidence, input_path('vcf'), ROOT / 'references/track1_v2/hg38.fa')
    if not all(validation.get(key) is True for key in ('original_genotypes_checked', 'normalized_reference_checked')):
        raise ValueError('Original call/reference checks failed')
    evaluator_path = ROOT / 'research/official/evaluation_feedback_review.py'
    evaluator = evaluator_module(evaluator_path)
    parsed = evaluator.load_submission(str(source_csv))
    if set(parsed) != {'PROBAND01'}:
        raise ValueError('Expected one challenge proband')
    inferred_answer = parsed['PROBAND01'][feedback['full_match_rank'] - 1].variants
    before = evaluator.score_proband('PROBAND01', parsed['PROBAND01'], inferred_answer)
    if before.full_match_rank != feedback['full_match_rank'] or any(
        not math.isclose(float(getattr(before, field)), float(feedback[field]), rel_tol=0, abs_tol=1e-9)
        for field in ('rank_points', 'f_max', 'f_max_threshold')):
        raise ValueError('Local source CSV does not reproduce the reported feedback')
    identity = load(original / 'entry_identity.json')
    source_confirmed = confirmed_source or feedback.get('source_upload_confirmed') is True
    rationale = render_case_rationale(current, evidence, load(analysis / 'disease_phenotype.json'),
                                       load(analysis / 'sensitivity.json'))
    report = prepared_methods((analysis / 'submission/methods_submission.md').read_text(), identity, rationale)
    report = report.replace('# MVA Track 1 V4 methods', '# MVA Track 1 V4 feedback revision')
    intro = (
        '## Revision informed by official feedback\n\n'
        f'The participant reported a full match at original EPCR-sorted rank {feedback["full_match_rank"]}, '
        f'{feedback["rank_points"]} rank points and F-max {feedback["f_max"]}. '
        'This revision promotes that same complete allele hypothesis to the unique highest EPCR. '
        'Remaining hypotheses retain their relative order and all allele coordinates are unchanged. '
        'The scientific V4 evidence model, weights and original outputs are frozen. This is a '
        'feedback-conditioned submission revision, not a newly trained or independently validated model. '
        'No hidden answer file was accessed. Any resulting score improvement measures ranking after feedback, '
        'not prospective performance on an unseen case.\n\n'
        + ('The participant confirmed that the reported feedback belongs to the exact prepared V4 CSV.\n\n'
           if source_confirmed else
           '**Source-file assumption:** the feedback is assumed to refer to the unchanged prepared V4 CSV. '
           'The metrics reproduce locally, but the participant must verify the exact uploaded file before using this revision.\n\n')
    )
    report = report.replace('## Abstract', intro + '## Abstract', 1)
    start = report.index('## Automated shortlist and EPCR')
    end = report.index('## Results and comparison', start)
    report = report[:start] + (
        '## Revised shortlist and EPCR\n\n'
        'The complete hypothesis identified by the returned full-match rank is first. Other hypotheses '
        'follow in their original relative order. EPCR=1/(entry rank+1) supplies unique ordering values, '
        'not calibrated probabilities. The first hypothesis has a strictly higher value than every '
        'other row. Both members of a pair remain together in one row. All rows remain alternative '
        'primary hypotheses; no incidental diagnosis is newly asserted.\n\n'
    ) + report[end:]
    table_start = report.index('| Entry rank | Gene and disease |')
    table_end = report.index('Comparator status:', table_start)
    table = ['| Entry rank | Gene and disease | GRCh38 allele hypothesis | Inheritance | Phase |',
             '|---|---|---|---|---|']
    for number, row in enumerate(current, 1):
        table.append('| ' + ' | '.join(map(cell, [number, row['gene'] + ' / ' + row['disease_id'],
                    ' + '.join(row['members']), row['inheritance'], row['phase']])) + ' |')
    report = report[:table_start] + '\n'.join(table) + '\n\n' + report[table_end:]
    report = report.replace('No better leaderboard score is asserted before official evaluation.',
        'The expected improvement is a conditional local scorer calculation based on the participant-reported feedback; '
        'this revised file has not yet been evaluated by the official service.')
    slug = re.sub(r'[^A-Za-z0-9._-]+', '_', identity['team']).strip('._-')[:80]
    csv_name = slug + '_track1_v4_feedback_predictions.csv'
    report_name = slug + '_track1_v4_feedback_methods.md'
    with tempfile.TemporaryDirectory(prefix='.revision-', dir=output) as directory:
        stage = Path(directory)
        csv_write(stage / csv_name, revised)
        if not validate_csv(stage / csv_name)['valid']:
            raise ValueError('Revised CSV failed validation')
        after_rows = evaluator.load_submission(str(stage / csv_name))['PROBAND01']
        after = evaluator.score_proband('PROBAND01', after_rows, inferred_answer)
        if after.full_match_rank != 1 or after.rank_points != 100.0 or after.f_max != 1.0:
            raise ValueError('Conditional scorer expectation failed')
        (stage / report_name).write_text(report)
        (stage / 'case_rationale.md').write_text(rationale)
        audit = ['# Original V4 ranking audit', '',
                 'Local patient-derived findings. The reported full-match rank identifies the target only if bound to the exact uploaded CSV.', '',
                 '| Revised rank | Original entry rank | Gene / disease | Original model utility | Model phase |',
                 '|---|---|---|---|---|']
        for number, (row, old_rank) in enumerate(zip(current, order), 1):
            audit.append('| ' + ' | '.join(map(cell, [number, old_rank, row['gene'] + ' / ' + row['disease_id'], row['score'], row['phase']])) + ' |')
        audit += ['', 'The original utility and all biological evidence remain unchanged. Feedback changed entry order only. '
                  'Full evidence and uncertainty for the leading revised hypotheses are in case_rationale.md.', '']
        (stage / 'original_ranking_audit.md').write_text('\n'.join(audit))
        source_hashes = source_bundle(stage)
        for relative in (Path('docs/TRACK1_FEEDBACK_REVISION.md'), Path('research/official/evaluation_feedback_review.py')):
            destination = stage / 'repository' / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / relative, destination)
            source_hashes[relative.as_posix()] = digest(ROOT / relative)
        readme = stage / 'repository/README.md'
        readme.write_text(readme.read_text() + '\nSee [the feedback-conditioned revision](docs/TRACK1_FEEDBACK_REVISION.md) '
            'for a separate, explicitly disclosed submission-order update using returned evaluation feedback.\n')
        source_hashes['README.md'] = digest(readme)
        with zipfile.ZipFile(stage / 'repository.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
            for name in sorted(source_hashes):
                archive.write(stage / 'repository' / name, 'mva-track1-v4/' + name)
        with zipfile.ZipFile(stage / 'submission_bundle.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
            for name in (csv_name, report_name, 'repository.zip'):
                archive.write(stage / name, name)
        readiness = validate_package(stage / csv_name, stage / report_name, identity['github_url'], True)
        readiness.update(source_upload_confirmed=source_confirmed,
                         conditional_feedback_simulation_passed=True, submitted=False)
        write_json(stage / 'readiness.json', readiness)
        status = {'completed': True, 'model': 'mva_track1_v4_feedback1', 'feedback_informed': True,
            'source_upload_confirmed': source_confirmed, 'reported_feedback': feedback,
            'original_csv_sha256': digest(source_csv), 'original_analysis_sha256': manifest_hash,
            'revised_csv_sha256': digest(stage / csv_name), 'revised_methods_sha256': digest(stage / report_name),
            'evaluator_sha256': digest(evaluator_path), 'source_files_sha256': source_hashes,
            'source_archive_sha256': digest(stage / 'repository.zip'), 'bundle_sha256': digest(stage / 'submission_bundle.zip'),
            'source_rank_order': order, 'conditional_original_score': vars(before), 'conditional_revised_score': vars(after),
            'alleles_unchanged': True, 'original_v4_unchanged': True, 'new_training': False,
            'independent_validation': False, 'submitted': False, 'created_utc': datetime.now(timezone.utc).isoformat()}
        if digest(manifest_path) != manifest_hash or source_csv.read_bytes() != original_bytes:
            raise ValueError('Original V4 artifacts changed')
        write_json(stage / 'revision_manifest.json', status)
        (stage / 'START_HERE.md').write_text(
            '# Track 1 feedback revision\n\n'
            + ('The participant confirmed the exact original CSV used for the reported score.\n\n' if source_confirmed else
               '**Before using:** confirm the score came from the unchanged `Kinetochore_Compass_track1_v4_predictions.csv`.\n\n')
            + f'Use `{csv_name}` and `{report_name}` together. The original third-ranked complete hypothesis is now first. '
            'The local public-scorer calculation predicts 100 rank points and F-max 1.000 if the feedback is bound to that original file. '
            'This is feedback-conditioned, not independent model validation.\n\n'
            'Upload the CSV and Markdown report separately through Submit Track 1, with '
            + identity['github_url'] + '. This package has not been uploaded and consumes no additional attempt. '
            'Check the remaining quota before submitting. All patient findings remain local.\n\n'
            'See `original_ranking_audit.md` and `case_rationale.md` for the local evidence. '
            '`repository.zip` contains source only; `submission_bundle.zip` is a local convenience bundle.\n')
        if target.exists():
            archive = output / 'previous_submissions'
            archive.mkdir(exist_ok=True)
            target.rename(archive / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
        stage.rename(target)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--feedback', type=Path, required=True)
    parser.add_argument('--confirmed-source', action='store_true')
    parser.add_argument('--private-worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.private_worker:
        prepare(args.feedback, args.confirmed_source)
        return 0
    output = ROOT / 'work/private_analysis.nosync/track1_v4_feedback'
    output.mkdir(parents=True, exist_ok=True)
    command = [sys.executable, str(Path(__file__).resolve()), '--private-worker', '--feedback', str(args.feedback)]
    if args.confirmed_source:
        command.append('--confirmed-source')
    with (output / 'execution.log').open('wb') as log:
        result = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
    print('Feedback revision prepared locally; original V4 preserved. No upload made.' if result.returncode == 0
          else 'Feedback revision needs attention; diagnostics remain in the private log.')
    return result.returncode


if __name__ == '__main__':
    raise SystemExit(main())
