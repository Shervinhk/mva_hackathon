import csv
import argparse
import contextlib
import hashlib
import importlib.util
import json
import io
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('exomiser_wrapper', ROOT / 'scripts/run_track1_exomiser.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class FakeResponse:
    status_code = 206
    def __init__(self, data):
        self.content = data
    def raise_for_status(self):
        pass


class FakeClient:
    def __init__(self, data):
        self.data, self.calls = data, []
    def get(self, url, headers):
        begin, end = map(int, headers['Range'][6:].split('-'))
        self.calls.append((begin, end))
        return FakeResponse(self.data[begin:end + 1])


class ExomiserWrapperTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_build_analysis_uses_all_genes_and_valid_proband_terms(self):
        preset = self.root / 'preset.yml'
        preset.write_text('analysisMode: PASS_ONLY\nsteps:\n  - frequencyFilter: {maxFrequency: 2.0}\n')
        result = m.build_analysis(self.root / 'toy.vcf', ['HP:0000002', 'HP:0000001', 'HP:0000002'],
                                  self.root, preset, proband='TOY')
        self.assertEqual(result['analysis']['hpoIds'], ['HP:0000001', 'HP:0000002'])
        self.assertNotIn('genePanelFilter', str(result))
        self.assertEqual(result['analysis']['genomeAssembly'], 'hg38')
        self.assertEqual(result['outputOptions']['numGenes'], 0)
        with self.assertRaises(ValueError):
            m.build_analysis(self.root / 'toy.vcf', ['not-HPO'], self.root, preset)

    def test_provenance_tracks_phenotype_changes_and_stable_payload_serialization(self):
        files = [self.root / name for name in ('toy.vcf', 'toy.docx', 'preset.yml')]
        for path in files:
            path.write_text('synthetic input')
        before = m.input_hashes(*files)
        files[1].write_text('changed synthetic phenotype')
        after = m.input_hashes(*files)
        self.assertNotEqual(before['phenotype'], after['phenotype'])
        self.assertEqual(before['vcf'], after['vcf'])
        self.assertEqual(m.payload_sha256({'a': 1, 'b': 2}), m.payload_sha256({'b': 2, 'a': 1}))

    def test_failed_rerun_invalidates_previous_completion_before_reference_checks(self):
        out = self.root / 'results'
        out.mkdir()
        status_path = out / 'comparator_status.json'
        m.dump(status_path, {'completed': True, 'run_id': 'old'})
        def fail_verification(references):
            interim = json.loads(status_path.read_text())
            self.assertFalse(interim['completed'])
            self.assertNotEqual(interim['run_id'], 'old')
            raise RuntimeError('Synthetic missing reference')
        args = argparse.Namespace(output=out, references=self.root)
        with mock.patch.object(m, 'verify_data', side_effect=fail_verification), contextlib.redirect_stdout(io.StringIO()):
            result = m.run_comparator(args)
        self.assertEqual(result, 1)
        self.assertFalse(json.loads(status_path.read_text())['completed'])

    def test_only_local_inputs(self):
        with self.assertRaises(ValueError):
            m.local_path('https://example.test/private.vcf')

    def test_cloud_placeholder_is_rejected_before_content_read(self):
        with mock.patch.object(Path, 'stat', return_value=SimpleNamespace(st_flags=0x40000000)):
            with self.assertRaisesRegex(RuntimeError, 'cloud placeholder'):
                m.require_local_reference(self.root / 'public.db')
        with mock.patch.object(Path, 'stat', return_value=SimpleNamespace(st_flags=0)):
            m.require_local_reference(self.root / 'public.db')

    def test_archive_path_traversal(self):
        for name in ['../../outside', '/tmp/outside']:
            with self.assertRaises(ValueError):
                m.safe_member(self.root, name)

    def test_normalization_is_minimal_not_repeat_left_alignment(self):
        self.assertEqual(m.normalize_allele('1', 100, 'AC', 'AT'), 'chr1:101:C:T')
        self.assertEqual(m.normalize_allele('MT', 50, 'A', 'G'), 'chrM:50:A:G')
        with self.assertRaises(ValueError):
            m.normalize_allele('1', 100, 'N', '<DEL>')

    def test_native_and_python_logging_stays_local(self):
        log = self.root / 'private.log'
        with m.private_output(log):
            print('SYNTHETIC_PYTHON_LOG')
            os.write(2, b'SYNTHETIC_NATIVE_LOG\n')
        text = log.read_text()
        self.assertIn('SYNTHETIC_PYTHON_LOG', text)
        self.assertIn('SYNTHETIC_NATIVE_LOG', text)

    def test_public_download_resume_verifies_existing_chunks(self):
        content = b'synthetic-public-reference-contents'
        remote = type('FakeRemote', (), {'url': m.BASE_URL + 'toy.zip', 'client': FakeClient(content)})()
        target = self.root / 'toy.db'
        path, receipt = m.download_member(remote, 0, len(content), target, workers=2, chunk_size=8)
        self.assertEqual(path.read_bytes(), content)
        calls = len(remote.client.calls)
        m.download_member(remote, 0, len(content), target, workers=2, chunk_size=8)
        self.assertEqual(len(remote.client.calls), calls)
        with path.open('r+b') as handle:
            handle.write(b'BROKEN!!')
        m.download_member(remote, 0, len(content), target, workers=2, chunk_size=8)
        self.assertEqual(len(remote.client.calls), calls + 1)
        self.assertEqual(path.read_bytes(), content)

    def test_output_keeps_pair_together_without_claiming_phase(self):
        gene_fields = ['#RANK', 'GENE_SYMBOL', 'MOI', 'EXOMISER_GENE_COMBINED_SCORE', 'EXOMISER_GENE_PHENO_SCORE']
        variant_fields = ['GENE_SYMBOL', 'MOI', 'CONTRIBUTING_VARIANT', 'FILTER', 'CONTIG', 'START', 'REF', 'ALT']
        gene = self.root / 'genes.tsv'
        variants = self.root / 'variants.tsv'
        with gene.open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t'); writer.writerow(gene_fields)
            writer.writerow([1, 'TOY', 'AR', .8, .7])
        with variants.open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t'); writer.writerow(variant_fields)
            writer.writerows([['TOY', 'AR', '1', 'PASS', '1', '10', 'A', 'C'],
                              ['TOY', 'AR', '1', 'PASS', '1', '20', 'G', 'T'],
                              ['TOY', 'AR', '0', 'PASS', '1', '30', 'G', 'A']])
        output = self.root / 'hypotheses.json'
        self.assertTrue(m.export_hypotheses(gene, variants, output))
        result = json.loads(output.read_text())
        self.assertEqual(len(result['hypotheses']), 1)
        self.assertEqual(len(result['hypotheses'][0]['variant_ids']), 2)
        self.assertFalse(result['phase_confirmation'])
        self.assertFalse(result['hypotheses'][0]['score_is_calibrated_probability'])

    def test_reference_properties_disable_clinvar_whitelist(self):
        result = m.properties_text(self.root)
        self.assertIn('exomiser.hg38.use-clinvar-white-list=false', result)
        self.assertNotIn('exomiser.hg19.data-version=', result)

    def test_reference_verification_rejects_missing_size_and_hash_mismatches(self):
        member = 'toy/toy_variants.mv.db'
        path = self.root / member
        path.parent.mkdir()
        path.write_bytes(b'SYNTHETIC')
        manifest = {'url': m.BASE_URL + 'toy.zip', 'members': {
            member: {'bytes': path.stat().st_size, 'sha256': m.sha256(path)}}}
        m.dump(self.root / 'toy.zip.manifest.json', manifest)
        result = m.verify_archive_data(self.root, 'toy.zip', [member])
        self.assertIn(member, result['members'])
        with self.assertRaisesRegex(ValueError, 'entries are missing'):
            m.verify_archive_data(self.root, 'toy.zip', [member, 'toy/missing'])
        path.write_bytes(b'CORRUPTED')
        with self.assertRaisesRegex(ValueError, 'SHA-256'):
            m.verify_archive_data(self.root, 'toy.zip', [member])
        path.write_bytes(b'SHORT')
        with self.assertRaisesRegex(ValueError, 'size'):
            m.verify_archive_data(self.root, 'toy.zip', [member])


if __name__ == '__main__':
    unittest.main()
