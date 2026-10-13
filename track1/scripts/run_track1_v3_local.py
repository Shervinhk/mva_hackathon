"""Run v3 locally while keeping all detailed logs and findings private."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys
import time
import traceback

from track1_baseline import write_json
from run_track1_exomiser import require_local_reference

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fresh-annotation', action='store_true')
    parser.add_argument('--skip-comparator', action='store_true',
                        help='Run and package v3 without claiming a new Exomiser execution.')
    parser.add_argument('--review-only', action='store_true',
                        help='Prepare the manual review page without automated CSV export.')
    parser.add_argument('--output-root', type=Path,
                        default=ROOT / 'work/private_analysis.nosync',
                        help='Local private output root; verify exclusion from any cloud-sync service.')
    parser.add_argument('--reference-index', type=Path,
                        default=ROOT / 'references/track1_v3/public_annotations.sqlite')
    parser.add_argument('--exomiser-mvstore', type=Path,
                        default=ROOT / 'references/exomiser/data/2602_hg38/2602_hg38_variants.mv.db')
    args = parser.parse_args()
    out = args.output_root.resolve() / 'track1_v3'
    out.mkdir(parents=True, exist_ok=True)
    status_path = out / 'execution_status.json'
    status = {'started_utc': datetime.now(timezone.utc).isoformat(), 'completed': False,
              'steps': [], 'private_outputs_returned_to_assistant': False,
              'new_comparator_run_requested': not args.skip_comparator}
    started = time.perf_counter()
    # Replace any previous run's status before the first potentially long step.
    write_json(status_path, status)
    for path in (args.reference_index, args.exomiser_mvstore):
        if not path.is_file():
            write_json(status_path, {**status, 'failure': 'Required public reference missing; finish reference preparation.'})
            print('Public references are incomplete; private analysis was not started.', flush=True)
            return 2
        try:
            require_local_reference(path)
        except RuntimeError:
            write_json(status_path, {**status, 'failure': 'A required public reference is a cloud placeholder; restore its local contents.'})
            print('A public reference needs local restoration; private analysis was not started.', flush=True)
            return 2
    command = [sys.executable, 'scripts/track1_v3.py', '--gff',
               'references/track1_v2/MANE.v1.5.bcfcsq.gff3', '--reference-index',
               str(args.reference_index), '--exomiser-mvstore', str(args.exomiser_mvstore),
               '--output', str(out), '--annotated-vcf', str(out / 'annotated.bcf')]
    if not args.fresh_annotation:
        cached = ROOT / 'outputs/track1_v2/annotated.bcf'
        if cached.is_file():
            command += ['--annotated-vcf', str(cached), '--reuse-annotation']
    steps = [('analysis', command)]
    if not args.skip_comparator:
        steps.append(('independent_comparator', [sys.executable, 'scripts/run_track1_exomiser.py',
                      '--output', str(args.output_root.resolve() / 'track1_exomiser')]))
    steps.append(('submission_review_package', [sys.executable, 'scripts/build_track1_v3_submission.py',
                  '--analysis', str(out), '--output', str(args.output_root.resolve() / 'submission_track1_v3'),
                  '--comparator-output', str(args.output_root.resolve() / 'track1_exomiser')]
                  + ([] if args.review_only else ['--automated'])))
    essential_failed = False
    for name, command in steps:
        print('Running local step: ' + name, flush=True)
        with (out / (name + '_run.log')).open('w') as log:
            result = subprocess.run(command, cwd=ROOT, stdout=log, stderr=log, stdin=subprocess.DEVNULL)
        status['steps'].append({'step': name, 'success': result.returncode == 0,
                                'return_code': result.returncode})
        status['elapsed_seconds'] = time.perf_counter() - started
        write_json(status_path, status)
        if result.returncode:
            print('Local step requires attention; detailed diagnostics remain in its private log.', flush=True)
            if name != 'independent_comparator':
                essential_failed = True
                break
            # A comparator failure is disclosed in the package; do not discard
            # the main analysis or imply independent validation succeeded.
    status['completed'] = not essential_failed and all(step['success'] for step in status['steps'])
    status['main_analysis_and_package_completed'] = not essential_failed
    status['elapsed_seconds'] = time.perf_counter() - started
    write_json(status_path, status)
    if essential_failed:
        return 1
    print('Local analysis and review package completed. No upload was made.', flush=True)
    return 0 if status['completed'] else 3


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except Exception:
        directory = ROOT / 'work/private_analysis.nosync/track1_v3'
        directory.mkdir(parents=True, exist_ok=True)
        (directory / 'wrapper_error.log').write_text(traceback.format_exc())
        print('Local wrapper requires repair; private diagnostics remain on this computer.')
        raise SystemExit(1)
