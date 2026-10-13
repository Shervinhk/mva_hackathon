#!/usr/bin/env python3
"""Prepare an existing V4 entry locally without reranking or uploading.

The public entry point captures all worker output in the excluded private tree.
Only generic completion status is returned. This script does not read data on
import and never changes the completed scientific run or original submission.
"""
from pathlib import Path
import argparse
import csv
from datetime import datetime, timezone
import io
import json
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'scripts')]

V4_FILES = (
    '__init__.py', 'config.json', 'phenotypes.py', 'ranking.py', 'run.py',
    'run_local.py', 'export.py', 'README.md', 'VALIDATION.md',
    'tests/test_phenotypes.py', 'tests/test_ranking.py', 'tests/test_pipeline.py',
)


def read_json(path):
    return json.loads(path.read_text())


def check_hashes(root, hashes):
    from track1_baseline import digest
    if not isinstance(hashes, dict) or not hashes:
        raise ValueError('Missing provenance inventory')
    for name, expected in hashes.items():
        path = root / name
        if path.is_symlink() or not path.resolve().is_relative_to(root.resolve()) or not path.is_file():
            raise ValueError('Invalid provenance path')
        if digest(path) != expected:
            raise ValueError('Recorded artifact changed')


def entry_identity(supplied):
    """Require explicit plan/tools facts; never derive a tier from opt-out."""
    from track1_entry_identity import validate_identity
    identity, missing = validate_identity(supplied)
    details = supplied.get('ai_details', {})
    if not isinstance(details, dict):
        raise ValueError('AI details must be an object')
    details = {key: details.get(key, default) for key, default in {
        'provider': 'OpenAI', 'product': 'Codex', 'plan_tier': '',
        'data_controls': 'On 13 October 2026, the participant confirmed personally checking the relevant data-sharing/model-training settings and that these were off before any patient-data work began. Historical account settings and provider-side processing were not independently inspected by the assistant.',
        'other_ai_tools': '',
    }.items()}
    if any(not isinstance(value, str) for value in details.values()):
        raise ValueError('AI details must be text')
    from track1_entry_identity import PLACEHOLDER
    unknown = {'unknown', 'not sure', 'unsure', 'not provided', 'to confirm', 'n/a'}
    absent = [key for key, value in details.items()
              if not value.strip() or PLACEHOLDER.search(value) or value.strip().lower() in unknown]
    if not absent:
        identity['ai_disclosure'] = (
            f"Provider: {details['provider']}. Product: {details['product']}. "
            f"Plan/tier: {details['plan_tier']}. Data controls: {details['data_controls']} "
            f"Other AI tools: {details['other_ai_tools']}. "
            'Software and documentation assistance; no model was trained or fine-tuned by this pipeline. '
            'See the report for the prior-disclosure history and limits of account verification.'
        )
    identity, missing = validate_identity(identity)
    identity['ai_details'] = details
    return identity, sorted(set(missing + ['ai_details.' + key for key in absent]))


def source_bundle(stage):
    """Use a new directory and explicit source inventory, never old repo files."""
    from build_track1_v3_submission import source_archive
    if (stage / 'repository').exists():
        raise ValueError('Source staging must be fresh')
    for directory in ('scripts', 'tests', 'config'):
        if (ROOT / directory).is_symlink() or any(path.is_symlink() for path in (ROOT / directory).iterdir()):
            raise ValueError('Source inputs must not be symbolic links')
    source_archive(stage)
    repo = stage / 'repository'
    # V3's archive enumerates exactly the copied public allowlist.
    with zipfile.ZipFile(stage / 'repository.zip') as original:
        names = [Path(name) for name in original.namelist()]
        if any(path.is_absolute() or '..' in path.parts or len(path.parts) < 2
               or path.parts[0] != 'mva-track1-v3' for path in names):
            raise ValueError('Invalid original source archive path')
        allowed = {Path(*path.parts[1:]) for path in names}
    for name in V4_FILES:
        relative = Path('track1_v4') / name
        source = ROOT / relative
        if source.is_symlink() or not source.is_file():
            raise ValueError('Expected regular V4 public source')
        target = repo / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        allowed.add(relative)
    (repo / 'README.md').write_text(
        '# MVA Track 1 V4\n\nSee [the V4 method](track1_v4/README.md) and '
        '[entry preparation](docs/TRACK1_V4_ENTRY_PREPARATION.md).\n\n'
        'Source, configurations and fabricated tests only. Patient inputs, findings, '
        'identity details and local logs are excluded. Obtain authorized inputs and '
        'public references separately to reproduce V3 and then V4.\n')
    relative = Path('docs/TRACK1_V4_ENTRY_PREPARATION.md')
    shutil.copy2(ROOT / relative, repo / relative)
    allowed.add(relative)
    with zipfile.ZipFile(stage / 'repository.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
        for relative in sorted(allowed):
            path = repo / relative
            if path.is_symlink() or not path.is_file():
                raise ValueError('Invalid staged source')
            archive.write(path, Path('mva-track1-v4') / relative)
    actual = {path.relative_to(repo) for path in repo.rglob('*') if path.is_file()}
    if actual != allowed:
        raise ValueError('Source staging inventory mismatch')
    return {path.as_posix(): __import__('hashlib').sha256((repo / path).read_bytes()).hexdigest()
            for path in sorted(allowed)}


def checked_prediction_bytes(csv_path, selected, evidence, recorded_hash):
    from build_track1_v3_submission import csv_row
    from track1_baseline import digest
    from validate_track1_submission import FIELDS, validate_csv
    if digest(csv_path) != recorded_hash or not validate_csv(csv_path)['valid']:
        raise ValueError('Original prediction file failed validation')
    expected = io.StringIO(newline='')
    writer = csv.DictWriter(expected, fieldnames=FIELDS)
    writer.writeheader()
    for row, choice in selected:
        choice = dict(choice, reason='Automated V4 order; EPCR is an ordinal weight, not a causal probability; no human review claimed.')
        writer.writerow(csv_row(row, choice, evidence))
    contents = csv_path.read_bytes()
    if contents != expected.getvalue().encode():
        raise ValueError('CSV differs from completed V4 order')
    return contents


def prepared_methods(original, identity, rationale):
    from track1_entry_identity import participant_lines
    title, separator, body = original.partition('## Abstract')
    if not separator:
        raise ValueError('Expected V4 methods structure')
    report = '# MVA Track 1 V4 methods\n\n' + '\n'.join(participant_lines(identity)) + '\n\n' + separator + body
    marker = '## Validation and limitations'
    if marker not in report:
        raise ValueError('Expected methods validation section')
    report = report.replace(
        'The participant reports training opt-out was already enabled; plan/tier and account controls have not been independently verified.',
        'On 13 October 2026, the participant confirmed personally checking the relevant data-sharing/model-training settings and that these were off before any patient-data work began. This confirmation applies from the start of the work. The assistant did not independently inspect historical account settings or provider-side processing; the actual plan/tier remains to be supplied.')
    report = report.replace(
        'No retrospective no-disclosure or no-training guarantee is claimed.',
        'Earlier conversation input does not establish training use or imply that opt-out was enabled afterward. No claim of entirely offline operation or deletion of earlier conversation content is made.')
    ai_start, ai_end = report.find('## AI use and acknowledgement'), report.find('Dataset: ')
    if ai_start >= 0 and ai_end > ai_start:
        report = report[:ai_start] + (
            '## AI use and acknowledgement\n\n'
            'OpenAI Codex assisted with code development, public-source research and documentation. '
            'The participant confirms personally checking the relevant data-sharing/model-training settings '
            'and that these were off before any work with the patient data began. The final submission must '
            'specify the actual provider, plan/tier and any other AI tools in the participant disclosure above. '
            'V3 and V4 genomic processing and detailed results were handled locally. Some patient-related '
            'content appeared in earlier conversation and source review; the overall workflow was not entirely '
            'offline. This does not establish model-training use. Public participant methods were encountered '
            'during research, so this work is not described as blinded discovery.\n\n'
            '## Data-deletion commitment\n\n'
            'I will delete the hackathon data within 30 days of hackathon close from all environments '
            'under my control, including local machines, cloud instances, notebooks, private repositories, '
            'and intermediate or derived datasets. Once deletion is complete, I will email '
            'RarediseaserealkidMVAhackathon2026@synapse.org to confirm, as required by the organizers. '
            '[Official data-deletion rules](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/rules.py).\n\n'
            '## References and acknowledgement\n\n'
        ) + report[ai_end:]
    return report.replace(marker, rationale.rstrip() + '\n\n' + marker, 1)


def prepare(identity_path=None):
    from local_inputs import input_path
    from track1_baseline import digest, write_json
    from track1_v2 import install_network_guard
    from track1_v4.run import checked_base
    from track1_v4_case_report import render_case_rationale
    from track1_automated_export import automated_choices
    from build_track1_v3_submission import hypothesis_id, verify_calls
    from validate_track1_submission import validate_package

    install_network_guard()
    output = ROOT / 'work/private_analysis.nosync/track1_v4'
    original = output / 'submission'
    target = output / 'submission_ready'
    manifest_path = output / 'run_manifest.json'
    manifest = read_json(manifest_path)
    manifest_hash = digest(manifest_path)
    if (manifest.get('completed') is not True or manifest.get('model') != 'mva_track1_v4'
            or read_json(output / 'stage_status.json').get('completed') is not True):
        raise ValueError('Completed V4 run required')
    check_hashes(ROOT, manifest['source_hashes'])
    check_hashes(output, manifest['output_hashes'])
    if digest(ROOT / 'track1_v4/config.json') != manifest['v4_config_sha256']:
        raise ValueError('V4 configuration changed')
    base = Path(manifest['base_analysis'])
    if digest(base / 'run_manifest.json') != manifest['base_manifest_sha256']:
        raise ValueError('V3 baseline manifest changed')
    vcf, phenotype = input_path('vcf'), input_path('phenotype')
    fasta = ROOT / 'references/track1_v2/hg38.fa'
    _, current = checked_base(base, vcf, phenotype, fasta, ROOT / 'references/track1')
    if current != manifest['input_hashes']:
        raise ValueError('Inputs changed since V4')
    evidence = {row['id']: row for row in read_json(base / 'candidate_evidence.json')}
    ranked = read_json(output / 'ranked_hypotheses.json')
    if any(row.get('automatic_shortlist_eligible') is not True for row in ranked):
        raise ValueError('Deferred hypothesis in automatic ranking')
    selected, _ = automated_choices(ranked, manifest_hash, hypothesis_id,
        manifest['v4_config']['maximum_submission_hypotheses'])
    previous_status = read_json(original / 'package_status.json')
    prediction = checked_prediction_bytes(original / 'predictions.csv', selected, evidence,
        previous_status['prediction_csv_sha256'])
    validation = verify_calls([row for row, _ in selected], evidence, vcf, fasta)
    if not all(validation.get(key) is True for key in ('original_genotypes_checked', 'normalized_reference_checked')):
        raise ValueError('Original genotype/reference verification failed')
    # Compare the completed comparator's exact phenotype payload, VCF and output.
    comparator_matching = False
    for directory in (output / 'exomiser', ROOT / 'work/private_analysis.nosync/track1_exomiser'):
        status_path = directory / 'comparator_status.json'
        hypotheses = directory / 'comparator_hypotheses.json'
        if not status_path.is_file() or not hypotheses.is_file():
            continue
        status = read_json(status_path)
        comparator_matching |= (status.get('completed') is True
            and status.get('provenance_schema_version') == 1
            and status.get('inputs_unchanged_during_run') is True
            and status.get('reference_size_and_hash_verified') is True
            and status.get('analysis_input_hashes', {}).get('vcf') == current['vcf']
            and status.get('hpo_payload_sha256') == manifest['positive_hpo_payload_sha256']
            and status.get('output_hashes', {}).get('hypotheses') == digest(hypotheses))
    if not comparator_matching:
        raise ValueError('Previously completed matching comparator no longer verifies')
    identity_file = identity_path or (target / 'entry_identity.json' if target.is_dir() else original / 'entry_identity.json')
    identity, missing = entry_identity(read_json(identity_file))
    rationale = render_case_rationale([row for row, _ in selected], evidence,
        read_json(output / 'disease_phenotype.json'), read_json(output / 'sensitivity.json'))
    methods = prepared_methods((original / 'methods_submission.md').read_text(), identity, rationale)
    with tempfile.TemporaryDirectory(prefix='.entry-preparation-', dir=output) as temporary:
        stage = Path(temporary)
        (stage / 'predictions.csv').write_bytes(prediction)
        (stage / 'methods_submission.md').write_text(methods)
        (stage / 'case_rationale.md').write_text(rationale)
        shutil.copy2(original / 'comparison_v3_v4.md', stage / 'comparison_v3_v4.md')
        write_json(stage / 'entry_identity.json', identity)
        readiness = validate_package(stage / 'predictions.csv', stage / 'methods_submission.md', identity['github_url'] or None, True)
        readiness.update(participant_fields_missing=missing, technical_checks_passed=True,
            complete_local_entry=readiness['local_validation_passed'] and not missing,
            github_access_verified=False, ai_account_controls_verified=False, submitted=False)
        write_json(stage / 'readiness.json', readiness)
        source_hashes = source_bundle(stage)
        upload_files = ['predictions.csv', 'methods_submission.md']
        if identity['team'] and 'team' not in missing:
            slug = re.sub(r'[^A-Za-z0-9._-]+', '_', identity['team']).strip('._-')[:80]
            if slug:
                upload_files = [slug + '_track1_v4_predictions.csv', slug + '_track1_v4_methods.md']
                shutil.copy2(stage / 'predictions.csv', stage / upload_files[0])
                shutil.copy2(stage / 'methods_submission.md', stage / upload_files[1])
        with zipfile.ZipFile(stage / 'submission_bundle.zip', 'w', zipfile.ZIP_DEFLATED) as bundle:
            for name in (*upload_files, 'repository.zip'):
                bundle.write(stage / name, name)
        pending = ', '.join(missing) if missing else 'None in local metadata; remote checks remain.'
        (stage / 'START_HERE.md').write_text(
            '# Track 1 V4 entry\n\nThe completed model and prediction CSV are unchanged. '
            'Original calls, reference alleles, evidence provenance, matching Exomiser comparison and CSV format were checked again. '
            'The methods report now includes evidence-based explanations for the leading hypotheses generated locally. '
            'Software checks do not establish clinical correctness.\n\n'
            '## Files\n\n- Predictions: `' + upload_files[0] + '`\n- Methods: `' + upload_files[1] + '`\n'
            '- GitHub source: `repository/` or `repository.zip`\n'
            '- Local convenience bundle: `submission_bundle.zip`\n'
            '- Case explanation: `case_rationale.md`\n\n'
            '## Complete the entry details\n\nPending: ' + pending + '\n\n'
            'Edit `entry_identity.json`: team, real GitHub repository URL, actual AI plan/tier and other AI tools '
            '(write "None" only if that is accurate). The existing opt-out statement is participant-reported. '
            'Run from the project root:\n\n```bash\n.venv/bin/python scripts/prepare_track1_v4_entry.py\n```\n\n'
            'This updates the entry without rerunning annotation, ranking or Exomiser. Both upload files receive '
            'team/model names when a team is supplied. V4 is a software version, not a claim that four submissions were used.\n\n'
            '## Submit\n\nOpen the official [Submit Track 1 form](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026). '
            'Sign in with Hugging Face, check remaining attempts and rules, attach the CSV and Markdown report separately, '
            'and supply the real repository URL. Publish source files only from `repository/`; this folder excludes patient data/results. '
            'The repo may remain private during the hackathon but must be public afterward under the organizer rules. '
            'The bundle is for local convenience, not the form upload format. No upload, account-setting verification '
            'or challenge submission was performed here.\n')
        if (digest(manifest_path) != manifest_hash
                or digest(base / 'run_manifest.json') != manifest['base_manifest_sha256']
                or (original / 'predictions.csv').read_bytes() != prediction):
            raise ValueError('Original analysis changed during preparation')
        status = {'completed': True, 'technical_checks_passed': True,
            'predictions_unchanged': True, 'v3_unchanged': True, 'v4_analysis_unchanged': True,
            'matching_comparator_available': comparator_matching, 'csv_valid': True,
            'case_rationale_added': True, 'source_staged_fresh': True,
            'source_files_sha256': source_hashes, 'analysis_manifest_sha256': manifest_hash,
            'prediction_csv_sha256': digest(stage / 'predictions.csv'),
            'methods_sha256': digest(stage / 'methods_submission.md'),
            'source_archive_sha256': digest(stage / 'repository.zip'),
            'bundle_sha256': digest(stage / 'submission_bundle.zip'),
            'participant_fields_missing': missing, 'complete_local_entry': readiness['complete_local_entry'],
            'clinical_accuracy_established': False, 'submitted': False,
            'created_utc': datetime.now(timezone.utc).isoformat()}
        write_json(stage / 'preparation_status.json', status)
        if target.exists():
            archive = output / 'previous_ready_packages'
            archive.mkdir(exist_ok=True)
            target.rename(archive / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
        stage.rename(target)
    return status


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--identity', type=Path)
    parser.add_argument('--private-worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.private_worker:
        prepare(args.identity)
        return 0
    output = ROOT / 'work/private_analysis.nosync/track1_v4'
    output.mkdir(parents=True, exist_ok=True)
    command = [sys.executable, str(Path(__file__).resolve()), '--private-worker']
    if args.identity:
        command += ['--identity', str(args.identity)]
    with (output / 'entry_preparation.log').open('wb') as log:
        result = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
    print('V4 entry prepared locally. Predictions unchanged; no upload made.' if result.returncode == 0
          else 'V4 preparation needs attention. Diagnostics remain in the private log.')
    return result.returncode


if __name__ == '__main__':
    raise SystemExit(main())
