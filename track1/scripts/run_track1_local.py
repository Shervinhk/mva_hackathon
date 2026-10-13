"""Run private analysis with all detailed output confined to local log files."""
from __future__ import annotations

if __package__:
    from .local_inputs import input_path
else:
    from local_inputs import input_path
import argparse, json, subprocess, sys, time
from pathlib import Path
from track1_baseline import extract_phenotypes,write_json

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fresh-annotation',action='store_true')
    args=parser.parse_args();out=ROOT/'outputs/track1_v2';out.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter()
    # Private extraction stays inside this process, with no values printed.
    write_json(out/'phenotypes.json',extract_phenotypes(input_path('phenotype')))
    model=[sys.executable,'scripts/track1_v2.py','--gff','references/track1_v2/MANE.v1.5.bcfcsq.gff3']
    if not args.fresh_annotation:model.append('--reuse-annotation')
    steps=[('analysis',model),('clinvar_screen',[sys.executable,'scripts/track1_clinvar_screen.py','--phenotypes',str(out/'phenotypes.json')]),
           ('submission_package',[sys.executable,'scripts/build_track1_submission.py'])]
    statuses=[]
    for name,command in steps:
        print('Running local step: '+name,flush=True)
        with (out/(name+'_run.log')).open('w') as log:
            result=subprocess.run(command,cwd=ROOT,stdout=log,stderr=log,stdin=subprocess.DEVNULL)
        statuses.append({'step':name,'success':result.returncode==0})
        write_json(out/'execution_status.json',{'steps':statuses,'elapsed_seconds':time.perf_counter()-start})
        if result.returncode:
            print('Local step requires repair; private diagnostics retained in its log.',flush=True)
            return result.returncode
    print('Local analysis and submission package completed. No data were uploaded.',flush=True)
    return 0

if __name__=='__main__':
    try:raise SystemExit(main())
    except Exception:
        # Avoid automatically exposing file content or private exceptions in a tool response.
        import traceback
        p=ROOT/'outputs/track1_v2';p.mkdir(parents=True,exist_ok=True)
        (p/'wrapper_error.log').write_text(traceback.format_exc())
        print('Local wrapper requires repair; private diagnostics retained locally.')
        raise SystemExit(1)
