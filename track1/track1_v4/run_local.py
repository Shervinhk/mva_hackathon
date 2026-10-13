"""Public entry point: capture all private Python/native/subprocess output."""
from pathlib import Path
import argparse
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path)
    parser.add_argument('--phenotype-override', type=Path)
    parser.add_argument('--identity', type=Path)
    parser.add_argument('--refresh-comparator', action='store_true')
    parser.add_argument('--skip-comparator', action='store_true')
    args = parser.parse_args()
    output = ROOT/'work/private_analysis.nosync/track1_v4'
    output.mkdir(parents=True,exist_ok=True)
    command = [sys.executable,str(ROOT/'track1_v4/run.py')]
    for key in ('base','phenotype_override','identity'):
        if getattr(args,key) is not None:
            command += ['--'+key.replace('_','-'), str(getattr(args,key))]
    for key in ('refresh_comparator','skip_comparator'):
        if getattr(args,key):
            command.append('--'+key.replace('_','-'))
    print('Running Track 1 V4 locally using existing V3 evidence.',flush=True)
    with (output/'execution.log').open('wb') as log:
        result = subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL)
    print('Track 1 V4 completed; private files saved locally. No upload was made.' if result.returncode == 0
          else 'Track 1 V4 needs attention; diagnostics remain in the private execution log.',flush=True)
    return result.returncode


if __name__ == '__main__':
    raise SystemExit(main())
