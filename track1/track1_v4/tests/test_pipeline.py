"""V4 provenance, packaging and scorer regressions using fabricated files only.

No patient inputs, private configuration, cached analysis or real findings are
opened. The official evaluator is optional public source, never an answer key.
"""
import copy
import csv
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock
import zipfile

from track1_v4 import run as pipeline
from track1_v4 import export as package
from track1_v4 import run_local


PUBLIC_ROOT = Path(__file__).resolve().parents[2]


def save(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload))
    return pipeline.digest(path)


def base_fixture(root):
    base, refs = root / 'base', root / 'references'
    base.mkdir(); refs.mkdir(); (root / 'scripts').mkdir()
    paths = {}
    for name in ('vcf', 'phenotype', 'fasta'):
        paths[name] = root / (name + '.synthetic')
        paths[name].write_text('Fabricated ' + name)
    outputs = {name + '.json': save(base / (name + '.json'), {'synthetic': name})
               for name in ('candidate_evidence', 'ranked_hypotheses', 'review_hypotheses',
                            'disease_phenotype', 'sensitivity', 'review_counts', 'annotation_counts')}
    ref_hashes = {}
    for name in ('hp.obo', 'genes_to_phenotype.txt'):
        (refs / name).write_text('Fabricated public reference ' + name)
        ref_hashes[name] = pipeline.digest(refs / name)
    manifest = {'completed': True, 'provenance_schema_version': 1,
                'output_hashes': outputs,
                'input_hashes': {k: pipeline.digest(paths[k]) for k in ('vcf', 'fasta')},
                'phenotype_sha256': pipeline.digest(paths['phenotype']),
                'reference_hashes': ref_hashes}
    for name, field in (('track1_v3.py', 'code_sha256'), ('track1_v2.py', 'shared_v2_sha256'),
                        ('track1_baseline.py', 'shared_baseline_sha256')):
        path = root / 'scripts' / name
        path.write_text('# Fabricated source ' + name)
        manifest[field] = pipeline.digest(path)
    save(base / 'run_manifest.json', manifest)
    return base, refs, paths, manifest


def toy_hypotheses():
    evidence = {}
    for pos in (10, 20, 30):
        key = f'chr1:{pos}:A:C'
        evidence[key] = {'id': key, 'chrom': 'chr1', 'pos': pos, 'ref': 'A', 'alt': 'C'}
    pair = ['chr1:10:A:C', 'chr1:20:A:C']
    first = {'rank': 1, 'gene': 'SYNTHETIC1', 'disease_id': 'OMIM:000001',
             'inheritance': 'AR', 'phase': 'unknown', 'score': .7, 'members': pair,
             'automatic_shortlist_eligible': True, 'eligible_for_manual_selection': True}
    rows = [first, {**first, 'rank': 2, 'disease_id': 'OMIM:000002', 'members': pair[::-1]},
            {**first, 'rank': 3, 'gene': 'SYNTHETIC2', 'inheritance': 'AD',
             'phase': 'not_applicable', 'members': ['chr1:30:A:C']}]
    return rows, evidence


def fake_source_bundle(output):
    with zipfile.ZipFile(output / 'repository.zip', 'w') as archive:
        archive.writestr('README.md', 'Fabricated public source only.')


def export_fixture(root):
    rows, evidence = toy_hypotheses()
    output = root / 'analysis' / 'submission'
    output.parent.mkdir()
    validation = {'original_genotypes_checked': True, 'normalized_reference_checked': True,
                  'clinical_causality_verified': False, 'review_alleles_checked': 3}
    payloads = {'ranked_hypotheses.json': rows, 'review_hypotheses.json': [],
                'disease_phenotype.json': {}, 'sensitivity.json': {}, 'review_counts.json': {},
                'phenotypes.json': {'observations': []}, 'validation_evidence.json': validation}
    hashes = {name: save(output.parent / name, value) for name, value in payloads.items()}
    manifest = {'completed': True, 'provenance_schema_version': 1, 'model': 'mva_track1_v4',
                'elapsed_seconds': 1.5, 'output_hashes': hashes, 'source_hashes': {},
                'clinical_accuracy_established': False}
    save(output.parent / 'run_manifest.json', manifest)
    identity = {'team': 'Synthetic team', 'github_url': 'https://github.com/example/synthetic-v4',
                'ai_disclosure': 'Synthetic provider and plan; fabricated opt-out statement for testing only.'}
    return output, rows, evidence, manifest, identity, {'mode': 'unavailable', 'record': None}, validation


class CacheProvenanceTests(unittest.TestCase):
    def check(self, root, fixture):
        base, refs, paths, _ = fixture
        with mock.patch.object(pipeline, 'ROOT', root), \
                mock.patch.object(pipeline, 'input_path', side_effect=AssertionError('No private config access')):
            return pipeline.checked_base(base, paths['vcf'], paths['phenotype'], paths['fasta'], refs)

    def test_completed_unchanged_cache_is_reusable(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); fixture = base_fixture(root)
            manifest, current = self.check(root, fixture)
            self.assertEqual(manifest, fixture[3])
            self.assertEqual(set(current), {'vcf', 'fasta', 'phenotype'})

    def test_incomplete_or_wrong_schema_and_inventory_fail_closed(self):
        for mutation in ('incomplete', 'wrong_schema', 'missing_output', 'extra_output'):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as temp:
                root = Path(temp); fixture = base_fixture(root)
                base, _, _, manifest = fixture
                if mutation == 'incomplete':
                    manifest['completed'] = False
                elif mutation == 'wrong_schema':
                    manifest['provenance_schema_version'] = 0
                elif mutation == 'missing_output':
                    manifest['output_hashes'].pop('candidate_evidence.json')
                else:
                    manifest['output_hashes']['unexpected.json'] = 'untrusted'
                save(base / 'run_manifest.json', manifest)
                with self.assertRaises(ValueError):
                    self.check(root, fixture)

    def test_changed_evidence_inputs_references_or_shared_source_refuse_reuse(self):
        for target in ('candidate_evidence.json', 'vcf', 'fasta', 'phenotype', 'hp.obo',
                       'genes_to_phenotype.txt', 'track1_v3.py', 'track1_v2.py', 'track1_baseline.py'):
            with self.subTest(target=target), tempfile.TemporaryDirectory() as temp:
                root = Path(temp); fixture = base_fixture(root)
                base, refs, paths, _ = fixture
                path = (paths[target] if target in paths else root / 'scripts' / target
                        if target.endswith('.py') else base / target if target.endswith('.json') else refs / target)
                path.write_text('Fabricated changed bytes')
                with self.assertRaises(ValueError):
                    self.check(root, fixture)

    def test_comparator_requires_matching_query_provenance_and_output_hash(self):
        for mutation in (None, 'completed', 'schema', 'inputs_changed', 'reference_unverified',
                         'vcf', 'query', 'output'):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as temp:
                directory = Path(temp)
                result_path = directory / 'comparator_hypotheses.json'
                result_hash = save(result_path, {'hypotheses': []})
                current = {'vcf': 'fabricated-input-hash'}
                query = ['HP:0000001', 'HP:0000002']
                status = {'completed': True, 'provenance_schema_version': 1,
                          'inputs_unchanged_during_run': True, 'reference_size_and_hash_verified': True,
                          'analysis_input_hashes': current.copy(),
                          'hpo_payload_sha256': pipeline.payload_sha256(query),
                          'output_hashes': {'hypotheses': result_hash}}
                field = {'completed': 'completed', 'schema': 'provenance_schema_version',
                         'inputs_changed': 'inputs_unchanged_during_run',
                         'reference_unverified': 'reference_size_and_hash_verified'}.get(mutation)
                if field:
                    status[field] = False
                elif mutation == 'vcf':
                    status['analysis_input_hashes']['vcf'] = 'changed'
                elif mutation == 'query':
                    status['hpo_payload_sha256'] = pipeline.payload_sha256(query[:1])
                elif mutation == 'output':
                    save(result_path, {'hypotheses': [{'fabricated': True}]})
                save(directory / 'comparator_status.json', status)
                actual = pipeline.eligible_comparator(directory, current, query[::-1] + query[:1])
                if mutation is None:
                    self.assertIsNotNone(actual)
                else:
                    self.assertIsNone(actual)

    def test_missing_comparator_returns_unavailable(self):
        with tempfile.TemporaryDirectory() as temp:
            self.assertIsNone(pipeline.eligible_comparator(Path(temp), {'vcf': 'synthetic'}, ['HP:0000001']))


class ExportTests(unittest.TestCase):
    def build(self, args):
        with mock.patch.object(package, 'source_bundle', side_effect=fake_source_bundle):
            return package.build_package(*args[:4], args[1], *args[4:])

    def test_pair_dedup_tie_disclosure_schema_and_bundle(self):
        with tempfile.TemporaryDirectory() as temp:
            args = export_fixture(Path(temp)); output = args[0]
            result = self.build(args)
            self.assertTrue(result['completed']); self.assertFalse(result['submitted'])
            self.assertFalse(result['clinical_accuracy_established'])
            with (output / 'predictions.csv').open(newline='') as handle:
                reader = csv.DictReader(handle); rows = list(reader)
                self.assertEqual(tuple(reader.fieldnames), package.FIELDS)
            self.assertEqual(len(rows), 2)
            self.assertEqual((rows[0]['pos_1'], rows[0]['pos_2']), ('10', '20'))
            self.assertEqual(rows[1]['pos_2'], '')
            self.assertEqual([float(r['epcr']) for r in rows], [.5, 1/3])
            self.assertTrue(all(r['finding_type'] == 'primary' for r in rows))
            self.assertTrue(package.validate_csv(output / 'predictions.csv')['valid'])
            comparison = (output / 'comparison_v3_v4.md').read_text()
            self.assertIn('Utility tie group size', comparison)
            self.assertIn('| 2 | unknown |', comparison)
            self.assertIn('do not establish confidence differences', comparison)
            selection = json.loads((output / 'selection.automated.json').read_text())
            self.assertEqual(selection['analysis_sha256'], pipeline.digest(output.parent / 'run_manifest.json'))
            self.assertFalse(selection['reviewed_locally'])
            self.assertFalse(selection['case_level_human_review'])
            self.assertTrue(all(not x['evidence_reviewed'] for x in selection['choices']))
            with zipfile.ZipFile(output / 'submission_bundle.zip') as archive:
                self.assertEqual(set(archive.namelist()), {'predictions.csv', 'methods_submission.md', 'repository.zip'})

    def test_admin_fields_remain_explicit_and_provider_settings_not_attested(self):
        for supplied in (True, False):
            with self.subTest(supplied=supplied), tempfile.TemporaryDirectory() as temp:
                args = list(export_fixture(Path(temp)))
                if not supplied:
                    args[4] = {}
                result = self.build(args)
                report = (args[0] / 'methods_submission.md').read_text()
                readiness = json.loads((args[0] / 'readiness.json').read_text())
                self.assertIn('have not been independently verified', report)
                self.assertEqual(readiness['local_validation_passed'], supplied)
                if supplied:
                    self.assertIn('Synthetic team', report)
                    self.assertEqual(result['participant_fields_missing'], [])
                    self.assertNotIn('[TO CONFIRM]', report)
                else:
                    self.assertIn('[TO CONFIRM]', report)
                    self.assertEqual(set(result['participant_fields_missing']), {'team', 'github_url', 'ai_disclosure'})

    def test_unverified_or_nonboolean_validation_flags_fail_closed(self):
        for flag in ('original_genotypes_checked', 'normalized_reference_checked'):
            for value in (False, None, 'false', 1):
                with self.subTest(flag=flag, value=value), tempfile.TemporaryDirectory() as temp:
                    args = list(export_fixture(Path(temp)))
                    args[-1][flag] = value
                    with self.assertRaises(ValueError):
                        self.build(args)
                    self.assertFalse((args[0] / 'predictions.csv').exists())

    def test_incomplete_or_changed_local_manifest_refuses_export(self):
        for mutation in ('incomplete', 'different_manifest', 'changed_ranked_evidence'):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as temp:
                args = export_fixture(Path(temp)); output = args[0]
                local = copy.deepcopy(args[3])
                if mutation == 'incomplete':
                    local['completed'] = False
                    args[3]['completed'] = False
                elif mutation == 'different_manifest':
                    local['elapsed_seconds'] = 123
                else:
                    save(output.parent / 'ranked_hypotheses.json', [])
                save(output.parent / 'run_manifest.json', local)
                with self.assertRaises(ValueError):
                    self.build(args)
                self.assertFalse((output / 'predictions.csv').exists())

    def test_invalid_csv_is_removed_and_source_bundle_not_created(self):
        with tempfile.TemporaryDirectory() as temp:
            args = export_fixture(Path(temp))
            with mock.patch.object(package, 'validate_csv', return_value={'valid': False}), \
                    mock.patch.object(package, 'source_bundle') as archive:
                with self.assertRaises(ValueError):
                    package.build_package(*args[:4], args[1], *args[4:])
            self.assertFalse((args[0] / 'predictions.csv').exists())
            archive.assert_not_called()


class WrapperTests(unittest.TestCase):
    def test_wrapper_routes_native_output_only_to_local_log(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            argv = ['run_local', '--base', str(root / 'synthetic_base'), '--skip-comparator']
            def fake_run(command, **kwargs):
                self.assertEqual(kwargs['cwd'], root)
                self.assertEqual(kwargs['stderr'], run_local.subprocess.STDOUT)
                self.assertEqual(kwargs['stdin'], run_local.subprocess.DEVNULL)
                self.assertIn('--skip-comparator', command)
                kwargs['stdout'].write(b'SYNTHETIC_PRIVATE_OUTPUT_MUST_NOT_BE_RELAYED')
                return mock.Mock(returncode=0)
            captured = io.StringIO()
            with mock.patch.object(run_local, 'ROOT', root), mock.patch.object(sys, 'argv', argv), \
                    mock.patch.object(run_local.subprocess, 'run', side_effect=fake_run), \
                    mock.patch.object(sys, 'stdout', captured):
                self.assertEqual(run_local.main(), 0)
            self.assertNotIn('SYNTHETIC_PRIVATE_OUTPUT', captured.getvalue())
            self.assertIn('No upload', captured.getvalue())


class OfficialScorerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = PUBLIC_ROOT / 'research/official/evaluation.py'
        if not path.is_file():
            raise unittest.SkipTest('Optional public official scorer snapshot absent from source archive')
        name = '_track1_v4_public_scorer_test'
        spec = importlib.util.spec_from_file_location(name, path)
        cls.scorer = importlib.util.module_from_spec(spec)
        sys.modules[name] = cls.scorer
        spec.loader.exec_module(cls.scorer)

    def test_complete_pair_gets_full_rank_credit_and_two_single_rows_do_not(self):
        m = self.scorer
        a, b = ('chr1', 10, 'A', 'C'), ('chr1', 20, 'A', 'C')
        truth = frozenset((a, b))
        pair = m.score_proband('SYNTHETIC', [m.SubmissionRow(truth, .5, 1)], truth)
        singles = m.score_proband('SYNTHETIC', [m.SubmissionRow(frozenset((a,)), .5, 1),
                                               m.SubmissionRow(frozenset((b,)), 1/3, 2)], truth)
        self.assertEqual(pair.rank_points, 100); self.assertEqual(pair.f_max, 1)
        self.assertEqual(singles.rank_points, 50); self.assertEqual(singles.f_max, 1)

    def test_late_complete_match_overrides_earlier_partial_rank_credit(self):
        m = self.scorer
        a, b = ('chr1', 10, 'A', 'C'), ('chr1', 20, 'A', 'C')
        truth = frozenset((a, b))
        partial = [m.SubmissionRow(frozenset((a,)), .5, 1)]
        full = partial + [m.SubmissionRow(frozenset((('chr1', 100+n, 'A', 'C'),)), 1/(n+1), n)
                          for n in range(2, 10)] + [m.SubmissionRow(truth, 1/11, 10)]
        first = m.score_proband('SYNTHETIC', partial, truth)
        second = m.score_proband('SYNTHETIC', full, truth)
        self.assertEqual(first.rank_points, 50)
        self.assertEqual(second.rank_points, 10)
        self.assertGreaterEqual(second.f_max, first.f_max)

    def test_f_max_uses_union_of_alleles_and_tied_threshold_groups(self):
        m = self.scorer
        a, b, c = (('chr1', n, 'A', 'C') for n in (10, 20, 30))
        true_pair = frozenset((a, b))
        overlap = [m.SubmissionRow(frozenset((a, c)), .5, 1),
                   m.SubmissionRow(frozenset((b, c)), 1/3, 2)]
        self.assertAlmostEqual(m.score_proband('SYNTHETIC', overlap, true_pair).f_max, .8)
        tied = [m.SubmissionRow(frozenset((a,)), .5, 1), m.SubmissionRow(frozenset((c,)), .5, 2)]
        ordered = [tied[0], m.SubmissionRow(frozenset((c,)), 1/3, 2)]
        self.assertAlmostEqual(m.score_proband('SYNTHETIC', tied, frozenset((a,))).f_max, 2/3)
        self.assertEqual(m.score_proband('SYNTHETIC', ordered, frozenset((a,))).f_max, 1)


if __name__ == '__main__':
    unittest.main()
