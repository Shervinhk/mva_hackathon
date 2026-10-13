"""Private V4 execution. Invoke through run_local.py to capture native logs."""
from pathlib import Path
import argparse
import copy
import json
import subprocess
import sys
import time
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'scripts')]
import pysam
from local_inputs import input_path
from track1_baseline import digest, write_json
from track1_v2 import install_network_guard
from track1_v3 import augment_mechanism_profiles
from build_track1_v3_submission import verify_calls, allele_identity
from track1_automated_export import automated_choices
from build_track1_v3_submission import hypothesis_id
from run_track1_exomiser import payload_sha256
from track1_v4.phenotypes import load_observations, positive_terms, DiseaseMatcherV4
from track1_v4.ranking import rank, distinct_rows
from track1_v4.export import build_package


def checked_base(base, vcf, phenotype, fasta, references):
    manifest = json.loads((base / 'run_manifest.json').read_text())
    if manifest.get('completed') is not True or manifest.get('provenance_schema_version') != 1:
        raise ValueError('Completed provenance-bound V3 analysis required')
    expected = {name+'.json' for name in ['candidate_evidence', 'ranked_hypotheses', 'review_hypotheses',
        'disease_phenotype', 'sensitivity', 'review_counts', 'annotation_counts']}
    if set(manifest.get('output_hashes', {})) != expected:
        raise ValueError('V3 evidence inventory incomplete')
    for name, recorded in manifest['output_hashes'].items():
        if digest(base/name) != recorded:
            raise ValueError('V3 evidence changed after its recorded run')
    current = {'vcf': digest(vcf), 'fasta': digest(fasta), 'phenotype': digest(phenotype)}
    if (current['vcf'] != manifest['input_hashes']['vcf'] or current['fasta'] != manifest['input_hashes']['fasta']
            or current['phenotype'] != manifest['phenotype_sha256']):
        raise ValueError('Current inputs differ from V3 provenance; reuse refused')
    for name in ('hp.obo', 'genes_to_phenotype.txt'):
        if digest(references/name) != manifest['reference_hashes'][name]:
            raise ValueError('Used HPO reference differs from V3 provenance')
    for name, field in [('track1_v3.py','code_sha256'), ('track1_v2.py','shared_v2_sha256'), ('track1_baseline.py','shared_baseline_sha256')]:
        if digest(ROOT/'scripts'/name) != manifest[field]:
            raise ValueError('Shared V3 implementation differs from completed baseline')
    return manifest, current


def eligible_comparator(directory, current, query):
    status_path, result_path = directory/'comparator_status.json', directory/'comparator_hypotheses.json'
    if not status_path.exists() or not result_path.exists():
        return None
    status = json.loads(status_path.read_text())
    valid = (status.get('completed') is True and status.get('provenance_schema_version') == 1
        and status.get('inputs_unchanged_during_run') is True and status.get('reference_size_and_hash_verified') is True
        and status.get('analysis_input_hashes', {}).get('vcf') == current['vcf']
        and status.get('hpo_payload_sha256') == payload_sha256(sorted(set(query)))
        and status.get('output_hashes', {}).get('hypotheses') == digest(result_path))
    return {'status': status, 'hypotheses': json.loads(result_path.read_text()), 'directory': str(directory)} if valid else None


def compare_comparator(path, ranked, comparator, fasta_path):
    lines = ['# V4 Exomiser comparison', '', 'Local patient-derived findings.', '',
        'Mode: ' + comparator['mode'] + '.', '',
        'Comparison requires matching VCF provenance and the exact positive HPO payload. '
        'Shared evidence and software inputs mean agreement is not independent clinical validation.', '']
    if comparator.get('record'):
        with pysam.FastaFile(str(fasta_path)) as fasta:
            comp = [(r, allele_identity(r['variant_ids'], fasta)) for r in comparator['record']['hypotheses']['hypotheses']]
            lines += ['| V4 distinct rank | Gene | Exomiser exact allele rank | Exomiser gene rank |', '|---|---|---|---|']
            for n, row in enumerate(distinct_rows(ranked, 20), 1):
                key = allele_identity(row['members'], fasta)
                exact = [r['rank'] for r,k in comp if k == key]
                gene = [r['rank'] for r,_ in comp if r['gene'] == row['gene']]
                lines.append(f"| {n} | {row['gene']} | {min(exact) if exact else 'no exact exported match'} | {min(gene) if gene else 'not exported'} |")
    else:
        lines.append('No current matching completed comparison is claimed. Main V4 predictions remain separately available.')
    path.write_text('\n'.join(lines)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path, default=ROOT/'work/private_analysis.nosync/track1_v3')
    parser.add_argument('--output', type=Path, default=ROOT/'work/private_analysis.nosync/track1_v4')
    parser.add_argument('--phenotype-override', type=Path)
    parser.add_argument('--identity', type=Path)
    parser.add_argument('--refresh-comparator', action='store_true')
    parser.add_argument('--skip-comparator', action='store_true')
    args = parser.parse_args()
    install_network_guard()
    started = time.perf_counter()
    output = args.output.resolve()
    base = args.base.resolve()
    if output == base or base in output.parents or output in base.parents:
        raise ValueError('V4 output must be separate from the V3 baseline')
    output.mkdir(parents=True, exist_ok=True)
    # Preserve identity and all previous entry files while making a new run's
    # incomplete state explicit. Never leave an older CSV at the active path.
    prior_identity = output/'submission/entry_identity.json'
    saved_identity = json.loads(prior_identity.read_text()) if prior_identity.exists() else None
    if (output/'submission').exists():
        previous = output/'previous_submissions'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        previous.parent.mkdir(exist_ok=True)
        (output/'submission').rename(previous)
    (output/'submission').mkdir()
    if saved_identity is not None:
        write_json(output/'submission/entry_identity.json', saved_identity)
    write_json(output/'submission/package_status.json', {'completed':False,'csv_created':False,'submitted':False})
    write_json(output/'run_manifest.json', {'completed': False, 'model': 'mva_track1_v4'})
    write_json(output/'stage_status.json', {'stage': 'verifying_reusable_v3_evidence', 'completed': False})
    vcf, phenotype = input_path('vcf'), input_path('phenotype')
    fasta, references = ROOT/'references/track1_v2/hg38.fa', ROOT/'references/track1'
    base_manifest, current = checked_base(base, vcf, phenotype, fasta, references)
    base_hash = digest(base/'run_manifest.json')
    variants = json.loads((base/'candidate_evidence.json').read_text())
    v3_rows = json.loads((base/'ranked_hypotheses.json').read_text())
    config = copy.deepcopy(base_manifest['config'])
    config['model'] = 'mva_track1_v4'
    cfg = json.loads((ROOT/'track1_v4/config.json').read_text())
    observations = load_observations(phenotype, args.phenotype_override)
    write_json(output/'phenotypes.json', observations)
    matcher = DiseaseMatcherV4(references/'hp.obo', references/'genes_to_phenotype.txt', config['aliases'])
    query = positive_terms(observations, matcher.base.alts)
    if not query:
        raise ValueError('No unambiguous positive proband terms; review local phenotype record')
    genes = {a['gene'] for v in variants for a in v['annotations']}
    profiles = {g: matcher.gene(g, query, []) for g in sorted(genes)}
    augment_mechanism_profiles(profiles, matcher, query, config)
    write_json(output/'stage_status.json', {'stage': 'reranking_cached_evidence', 'completed': False})
    ranked, review, counts = rank(variants, profiles, config)
    if not ranked:
        write_json(output/'review_hypotheses.json', review)
        raise ValueError('No complete eligible V4 hypotheses; retain exploratory review')
    scenarios = {}
    for name, options in [('without_phenotype', {'pheno_weight': 0}), ('phenotype_0.2', {'pheno_weight': .2}),
        ('phenotype_0.5', {'pheno_weight': .5}), ('without_condition_assertions', {'no_clinvar': True}),
        ('AlphaMissense_weight_0.5', {'predictor_weight': .5}), ('unknown_phase_factor_0.8', {'unknown_phase_factor': .8})]:
        result, _, _ = rank(variants, profiles, config, **options)
        scenarios[name] = [{k:r[k] for k in ('rank','gene','disease_id','members','inheritance','phase','score')} for r in result]
    selected, _ = automated_choices(ranked, base_hash, hypothesis_id, cfg['maximum_submission_hypotheses'])
    validation = verify_calls([r for r,_ in selected], {v['id']:v for v in variants}, vcf, fasta)
    write_json(output/'validation_evidence.json', validation)
    hashes = {}
    for name, payload in [('ranked_hypotheses',ranked),('review_hypotheses',review),('disease_phenotype',profiles),
                           ('sensitivity',scenarios),('review_counts',counts)]:
        write_json(output/(name+'.json'), payload)
        hashes[name+'.json'] = digest(output/(name+'.json'))
    hashes['phenotypes.json'] = digest(output/'phenotypes.json')
    code = {p.relative_to(ROOT).as_posix():digest(p) for p in (ROOT/'track1_v4').rglob('*.py') if '__pycache__' not in p.parts}
    code.update({f'scripts/{name}':digest(ROOT/'scripts'/name) for name in
                 ('track1_v3.py','track1_v2.py','track1_baseline.py','build_track1_v3_submission.py',
                  'track1_automated_export.py','track1_entry_identity.py','validate_track1_submission.py','run_track1_exomiser.py')})
    manifest = {'completed': True, 'provenance_schema_version': 1, 'model': 'mva_track1_v4',
        'base_manifest_sha256': base_hash, 'reused_candidate_evidence_sha256': base_manifest['output_hashes']['candidate_evidence.json'],
        'base_analysis': str(base), 'input_hashes': current, 'config':config, 'v4_config':cfg,
        'v4_config_sha256':digest(ROOT/'track1_v4/config.json'), 'output_hashes':hashes, 'source_hashes':code,
        'phenotype_override_sha256':digest(args.phenotype_override) if args.phenotype_override else None,
        'positive_hpo_payload_sha256':payload_sha256(sorted(set(query))),
        'reference_hashes':base_manifest['reference_hashes'], 'annotation_reused':True,
        'public_annotation_reused':True, 'training':'none', 'network_used':False,
        'clinical_accuracy_established':False, 'elapsed_seconds':time.perf_counter()-started,
        'created_utc':datetime.now(timezone.utc).isoformat()}
    write_json(output/'run_manifest.json', manifest)
    write_json(output/'stage_status.json', {'stage':'comparator_and_packaging', 'completed':False})
    comparator = {'mode':'unavailable', 'record':None}
    if not args.skip_comparator:
        if not args.refresh_comparator:
            for path, mode in [(output/'exomiser','reused_matching_v4_comparator'),
                               (ROOT/'work/private_analysis.nosync/track1_exomiser','reused_matching_v3_comparator')]:
                record = eligible_comparator(path, current, query)
                if record:
                    comparator = {'mode':mode, 'record':record}
                    break
        if comparator['record'] is None:
            # Supply only validated positive terms to the existing local wrapper.
            termfile = output/'exomiser_positive_observations.json'
            write_json(termfile, [{'hpo_id':t, 'subject':'proband','present':True} for t in query])
            command = [sys.executable, str(ROOT/'scripts/run_track1_exomiser.py'), '--phenotypes',str(termfile), '--output',str(output/'exomiser')]
            with (output/'comparator_execution.log').open('w') as log:
                result = subprocess.run(command,cwd=ROOT,stdout=log,stderr=log)
            record = eligible_comparator(output/'exomiser',current,query) if result.returncode == 0 else None
            comparator = {'mode':'new_v4_comparator' if record else 'unavailable_after_local_run', 'record':record}
    write_json(output/'comparator_reuse.json', {'mode':comparator['mode'], 'matching_comparator_available':comparator['record'] is not None,
        'independent_clinical_validation':False})
    compare_comparator(output/'comparator_comparison.md',ranked,comparator,fasta)
    identity_path = args.identity or output/'submission/entry_identity.json'
    if args.identity and not args.identity.is_file():
        raise ValueError('Supplied identity file does not exist')
    if not identity_path.exists():
        identity_path = ROOT/'work/private_analysis.nosync/submission_track1_v3/entry_identity.json'
    identity = json.loads(identity_path.read_text()) if identity_path.exists() else {}
    build_package(output/'submission',ranked,variants,manifest,v3_rows,identity,comparator,validation)
    if digest(base/'run_manifest.json') != base_hash:
        raise ValueError('V3 baseline changed during V4 execution')
    write_json(output/'stage_status.json', {'stage':'complete','completed':True,'v3_baseline_unchanged':True,
        'prediction_file_created':True,'submitted':False})


if __name__ == '__main__':
    main()
