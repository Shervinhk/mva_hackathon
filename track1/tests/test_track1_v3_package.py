"""Fabricated examples only; no challenge input or output is opened."""
import copy
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
import zipfile
import pysam

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import build_track1_v3_submission as package
import track1_v3 as analysis


def toy():
    row = {'rank': 1, 'gene': 'TOY1', 'disease_id': 'OMIM:000001', 'inheritance': 'AR',
           'members': ['chr1:10:A:C', 'chr1:20:G:T'], 'phase': 'unknown',
           'score': .71, 'flags': ['synthetic_example'], 'review_tier': 'review'}
    evidence = {
        'chr1:10:A:C': {'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'C'},
        'chr1:20:G:T': {'chrom': 'chr1', 'pos': 20, 'ref': 'G', 'alt': 'T'},
    }
    return [row], evidence


def reviewed(rows):
    selection = package.review_manifest(rows, 'toyhash')
    selection['reviewed_locally'] = True
    selection['choices'][0].update(selected=True, epcr=.4, evidence_reviewed=True,
                                    reason='Fabricated reviewed evidence for testing only.')
    return selection


class ReviewGateTests(unittest.TestCase):
    def test_automated_package_uses_identity_and_limits_bundle_contents(self):
        rows, evidence = toy()
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            analysis_dir, output = root / 'analysis', root / 'package'
            analysis_dir.mkdir(); output.mkdir(); (root / 'docs').mkdir()
            (root / 'docs/TRACK1_V3.md').write_text('Synthetic method description.')
            inputs = {}
            for name in ('vcf', 'phenotype', 'fasta'):
                path = root / name
                path.write_text('fabricated ' + name)
                inputs[name] = path
            records = [{**value, 'id': key} for key, value in evidence.items()]
            payloads = {'candidate_evidence': records, 'ranked_hypotheses': rows,
                        'review_hypotheses': [], 'disease_phenotype': {}, 'sensitivity': {},
                        'review_counts': {}, 'annotation_counts': {}}
            hashes = {}
            for name, value in payloads.items():
                path = analysis_dir / (name + '.json')
                path.write_text(json.dumps(value)); hashes[path.name] = package.digest(path)
            manifest = {'completed': True, 'provenance_schema_version': 1, 'output_hashes': hashes,
                        'input_hashes': {k: package.digest(inputs[k]) for k in ('vcf', 'fasta')},
                        'phenotype_sha256': package.digest(inputs['phenotype']), 'config': {}}
            (analysis_dir / 'run_manifest.json').write_text(json.dumps(manifest))
            identity = {'team': 'Synthetic team', 'github_url': 'https://github.com/example/synthetic-mva',
                        'ai_disclosure': 'Fabricated provider and tier; no real account setting is asserted.'}
            identity_path = root / 'identity.json'; identity_path.write_text(json.dumps(identity))
            def fake_archive(destination):
                with zipfile.ZipFile(destination / 'repository.zip', 'w') as archive:
                    archive.writestr('README.md', 'Fabricated public source only.')
                return {'files': 1, 'sha256': package.digest(destination / 'repository.zip')}
            argv = ['package', '--analysis', str(analysis_dir), '--output', str(output),
                    '--comparator-output', str(root / 'missing_comparator'), '--automated',
                    '--identity', str(identity_path)]
            for key, path in inputs.items():
                argv += ['--' + key, str(path)]
            with mock.patch.object(sys, 'argv', argv), mock.patch.object(package, 'ROOT', root), \
                    mock.patch.object(package, 'input_path', return_value=root / 'unused'), \
                    mock.patch.object(package, 'install_network_guard'), \
                    mock.patch.object(package, 'verify_calls', return_value={'synthetic_validation': True}), \
                    mock.patch.object(package, 'write_comparison', side_effect=lambda path, *a: path.write_text('Synthetic comparison.')), \
                    mock.patch.object(package, 'source_archive', side_effect=fake_archive):
                package.main()
            status = json.loads((output / 'package_status.json').read_text())
            self.assertTrue(status['automated_csv_exported_this_run'])
            self.assertFalse(status['reviewed_csv_exported_this_run'])
            self.assertTrue(status['identity_and_disclosure_complete'])
            self.assertFalse(status['provider_settings_independently_verified'])
            self.assertTrue(package.validate_csv(output / 'predictions.csv')['valid'])
            self.assertFalse(json.loads((output / 'selection.automated.json').read_text())['reviewed_locally'])
            report = (output / 'methods_submission.md').read_text()
            self.assertIn('Synthetic team', report)
            self.assertNotIn('[TO CONFIRM]', report)
            self.assertIn('not estimated or calibrated probabilities', report)
            with zipfile.ZipFile(output / 'submission_bundle.zip') as archive:
                self.assertEqual(set(archive.namelist()), {'predictions.csv', 'methods_submission.md', 'repository.zip'})

    def test_source_archive_excludes_private_inputs_and_legacy_document(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / 'source'
            output = Path(temp) / 'output'
            output.mkdir()
            for folder in ('scripts', 'tests', 'config', 'docs', 'research', 'data', 'work', 'examples'):
                (root / folder).mkdir(parents=True)
            for relative in ('scripts/example.py', 'tests/example.py', 'config/example.json',
                             'docs/TRACK1_V3.md', 'docs/AI_USE_AND_PRIVACY.md',
                             'requirements-track1.txt', '.gitignore', 'examples/synthetic_candidates.json'):
                (root / relative).write_text('fabricated public source')
            for relative in ('data/fixture.vcf', 'work/local_inputs.json', 'docs/TRACK1_V2.md'):
                (root / relative).write_text('fabricated excluded content')
            stale = output / 'repository/docs/TRACK1_V2.md'
            stale.parent.mkdir(parents=True)
            stale.write_text('fabricated superseded document')
            with mock.patch.object(package, 'ROOT', root):
                package.source_archive(output)
            self.assertFalse(stale.exists())
            with zipfile.ZipFile(output / 'repository.zip') as archive:
                names = archive.namelist()
            self.assertTrue(any(name.endswith('/scripts/example.py') for name in names))
            self.assertTrue(any(name.endswith('/examples/synthetic_candidates.json') for name in names))
            self.assertFalse(any('/data/' in name or '/work/' in name or 'TRACK1_V2.md' in name for name in names))

    def test_interrupted_rerun_invalidates_previous_analysis_manifest(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            manifest_path = directory / 'run_manifest.json'
            manifest_path.write_text(json.dumps({'completed': True, 'old_run': True}))
            argv = ['analysis', '--output', str(directory), '--annotated-vcf', str(directory / 'toy.bcf'),
                    '--gff', str(directory / 'toy.gff'), '--config', str(directory / 'missing_config.json')]
            with mock.patch.object(sys, 'argv', argv), \
                    mock.patch.object(analysis, 'install_network_guard'), \
                    mock.patch.object(analysis, 'input_path', return_value=directory / 'unused'):
                with self.assertRaises(FileNotFoundError):
                    analysis.main()
            manifest = json.loads(manifest_path.read_text())
            self.assertFalse(manifest['completed'])
            self.assertNotIn('old_run', manifest)

    def test_packaging_requires_complete_unchanged_evidence_and_reference(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            fasta = directory / 'toy.fa'
            fasta.write_text('fabricated reference bytes')
            expected = {'vcf': 'synthetic_vcf_hash', 'phenotype': 'synthetic_phenotype_hash'}
            names = ('candidate_evidence', 'ranked_hypotheses', 'review_hypotheses',
                     'disease_phenotype', 'sensitivity', 'review_counts', 'annotation_counts')
            hashes = {}
            for name in names:
                path = directory / (name + '.json')
                path.write_text('[]')
                hashes[path.name] = package.digest(path)
            manifest = {'provenance_schema_version': 1, 'completed': True,
                        'input_hashes': {'vcf': expected['vcf'], 'fasta': package.digest(fasta)},
                        'phenotype_sha256': expected['phenotype'], 'output_hashes': hashes}
            manifest_path = directory / 'run_manifest.json'
            manifest_path.write_text(json.dumps(manifest))
            self.assertEqual(package.validate_analysis_artifacts(directory, expected, fasta), manifest)
            for completed in (False, None):
                manifest_path.write_text(json.dumps({**manifest, 'completed': completed}))
                with self.assertRaisesRegex(ValueError, 'incomplete'):
                    package.validate_analysis_artifacts(directory, expected, fasta)
            manifest_path.write_text(json.dumps(manifest))
            for name in ('candidate_evidence.json', 'ranked_hypotheses.json'):
                path = directory / name
                path.write_text('["altered fabricated evidence"]')
                with self.assertRaisesRegex(ValueError, 'evidence changed'):
                    package.validate_analysis_artifacts(directory, expected, fasta)
                path.write_text('[]')
            fasta.write_text('different fabricated reference bytes')
            with self.assertRaisesRegex(ValueError, 'input provenance is stale'):
                package.validate_analysis_artifacts(directory, expected, fasta)

    def test_review_rebuild_archives_old_exports_even_when_analysis_is_incomplete(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            source, output = directory / 'analysis', directory / 'package'
            source.mkdir(); output.mkdir()
            (source / 'run_manifest.json').write_text(json.dumps({'completed': False}))
            for name in ('predictions.csv', 'selection.reviewed.json'):
                (output / name).write_text('previous fabricated export')
            for name in ('toy.vcf', 'toy.docx', 'toy.fa'):
                (directory / name).write_text('fabricated input')
            argv = ['package', '--analysis', str(source), '--output', str(output),
                    '--vcf', str(directory / 'toy.vcf'), '--phenotype', str(directory / 'toy.docx'),
                    '--fasta', str(directory / 'toy.fa')]
            with mock.patch.object(sys, 'argv', argv), \
                    mock.patch.object(package, 'install_network_guard'), \
                    mock.patch.object(package, 'input_path', return_value=directory / 'unused'):
                with self.assertRaisesRegex(ValueError, 'incomplete'):
                    package.main()
            for name in ('predictions.csv', 'selection.reviewed.json'):
                self.assertFalse((output / name).exists())
                preserved = list((output / 'superseded_exports').glob('*/' + name))
                self.assertEqual(len(preserved), 1)
                self.assertEqual(preserved[0].read_text(), 'previous fabricated export')
            self.assertFalse(json.loads((output / 'package_status.json').read_text())['completed'])

    def test_comparison_rejects_stale_inputs_and_changed_export(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'toy_comparator.json'
            path.write_text(json.dumps({'hypotheses': []}))
            expected = {'vcf': 'synthetic_vcf_hash', 'phenotype': 'synthetic_phenotype_hash'}
            status = {'completed': True, 'provenance_schema_version': 1,
                      'analysis_input_hashes': dict(expected), 'phenotype_input_kind': 'docx',
                      'inputs_unchanged_during_run': True, 'reference_size_and_hash_verified': True,
                      'output_hashes': {'hypotheses': package.digest(path)}}
            assessed, result = package.checked_comparator(status, path, expected)
            self.assertTrue(assessed['eligible_for_current_comparison'])
            self.assertEqual(result, {'hypotheses': []})
            for key in expected:
                stale = copy.deepcopy(status)
                stale['analysis_input_hashes'][key] = 'old'
                assessed, result = package.checked_comparator(stale, path, expected)
                self.assertFalse(assessed['eligible_for_current_comparison'])
                self.assertIsNone(result)
            path.write_text(json.dumps({'hypotheses': ['changed']}))
            assessed, result = package.checked_comparator(status, path, expected)
            self.assertIn('missing_or_changed_comparator_export', assessed['comparison_exclusion_reasons'])
            self.assertIsNone(result)

    def test_normalized_allele_and_original_source_call_are_checked_separately(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            fasta = directory / 'toy.fa'
            fasta.write_text('>chr1\n' + 'A' * 100 + '\n')
            pysam.faidx(str(fasta))
            vcf = directory / 'toy.vcf'
            vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=1,length=100>\n'
                          '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                          '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTOY\n'
                          '1\t20\t.\tAA\tA\t50\tPASS\t.\tGT\t0/1\n')
            compressed = directory / 'toy.vcf.gz'
            pysam.tabix_compress(str(vcf), str(compressed))
            pysam.tabix_index(str(compressed), preset='vcf')
            member = 'chr1:1:AA:A'
            row = {'members': [member]}
            evidence = {member: {'chrom': 'chr1', 'pos': 1, 'ref': 'AA', 'alt': 'A',
                                'source_contig': '1', 'source_pos': 20, 'source_ref': 'AA',
                                'source_alt': 'A', 'quality_pass': True}}
            self.assertTrue(package.verify_calls([row], evidence, compressed, fasta)['original_genotypes_checked'])
            evidence[member]['source_alt'] = 'AAA'
            with self.assertRaises(ValueError):
                package.verify_calls([row], evidence, compressed, fasta)

    def test_comparison_recognizes_equivalent_repeat_indels(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            fasta = directory / 'toy.fa'
            fasta.write_text('>chr1\n' + 'A' * 100 + '\n')
            pysam.faidx(str(fasta))
            rows, _ = toy()
            rows[0]['members'] = ['chr1:1:AA:A']
            previous = copy.deepcopy(rows)
            previous[0]['members'] = ['chr1:20:AA:A']
            output = directory / 'comparison.md'
            comparator = {'hypotheses': [{'gene': 'TOY1', 'rank': 2, 'variant_ids': ['chr1:10:AA:A']}]}
            package.write_comparison(output, rows, previous, comparator, fasta)
            result = output.read_text()
            self.assertIn('| 1 | TOY1: chr1:20:AA:A | 1 |', result)
            self.assertIn('| 1 | TOY1 | 2 | 2 |', result)

    @unittest.skipUnless(shutil.which('node'), 'Node is needed for generated review-page behavior checks')
    def test_generated_javascript_exports_reviewed_valid_csv(self):
        rows, evidence = toy()
        with tempfile.TemporaryDirectory() as temp:
            page = Path(temp) / 'review.html'
            package.create_review_html(page, rows, evidence, package.review_manifest(rows, 'toyhash'))
            result = subprocess.run(['node', str(ROOT / 'tests/review_page_harness.cjs'), str(page)],
                                    check=False, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            generated = json.loads(result.stdout)
            self.assertEqual(len(package.validate_choices(generated['selection'], rows, 'toyhash')), 1)
            csv_path = Path(temp) / 'predictions.csv'
            csv_path.write_text(generated['csv'])
            self.assertTrue(package.validate_csv(csv_path)['valid'])

    def test_template_does_not_invent_selection_or_probability(self):
        rows, _ = toy()
        selection = package.review_manifest(rows, 'toyhash')
        self.assertFalse(selection['reviewed_locally'])
        self.assertIsNone(selection['choices'][0]['epcr'])
        self.assertFalse(selection['choices'][0]['selected'])
        with self.assertRaises(ValueError):
            package.validate_choices(selection, rows, 'toyhash')

    def test_explicit_reviewed_choice_can_be_exported(self):
        rows, evidence = toy()
        selection = reviewed(rows)
        result = package.validate_choices(selection, rows, 'toyhash')
        self.assertEqual(len(result), 1)
        exported = package.csv_row(*result[0], evidence)
        self.assertEqual(exported['epcr'], .4)
        self.assertNotEqual(exported['epcr'], rows[0]['score'])
        self.assertEqual(exported['pos_2'], 20)

    def test_stale_manifest_is_rejected(self):
        rows, _ = toy()
        with self.assertRaises(ValueError):
            package.validate_choices(reviewed(rows), rows, 'changed_analysis')

    def test_incomplete_genotype_cannot_be_exported_as_diagnosis(self):
        rows, _ = toy()
        rows[0]['eligible_for_manual_selection'] = False
        with self.assertRaises(ValueError):
            package.validate_choices(reviewed(rows), rows, 'toyhash')

    def test_unknown_and_duplicate_choices_rejected(self):
        rows, _ = toy()
        selection = reviewed(rows)
        selection['choices'].append(copy.deepcopy(selection['choices'][0]))
        with self.assertRaises(ValueError):
            package.validate_choices(selection, rows, 'toyhash')
        selection = reviewed(rows)
        selection['choices'][0]['hypothesis_id'] = 'invented'
        with self.assertRaises(ValueError):
            package.validate_choices(selection, rows, 'toyhash')

    def test_bad_confidence_or_incomplete_review_rejected(self):
        rows, _ = toy()
        for value in [None, True, '0.4', 0, -1, 1.1, float('nan'), float('inf')]:
            selection = reviewed(rows)
            selection['choices'][0]['epcr'] = value
            with self.subTest(value=value), self.assertRaises(ValueError):
                package.validate_choices(selection, rows, 'toyhash')
        for key, value in [('evidence_reviewed', False), ('reason', ''), ('finding_type', 'incidental')]:
            selection = reviewed(rows)
            selection['choices'][0][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                package.validate_choices(selection, rows, 'toyhash')

    def test_pair_identity_order_does_not_change_manifest_key(self):
        rows, _ = toy()
        row = copy.deepcopy(rows[0])
        row['members'].reverse()
        self.assertEqual(package.hypothesis_id(row), package.hypothesis_id(rows[0]))

    def test_different_disease_profiles_do_not_duplicate_submission(self):
        rows, _ = toy()
        alternative = copy.deepcopy(rows[0])
        alternative['disease_id'] = 'OMIM:000002'
        rows.append(alternative)
        selection = reviewed(rows)
        selection['choices'][1].update(selected=True, epcr=.3, evidence_reviewed=True,
            reason='A different synthetic disease profile for the same allele pair.')
        with self.assertRaises(ValueError):
            package.validate_choices(selection, rows, 'toyhash')
        choices = package.distinct_review_rows(rows)
        self.assertEqual(len(choices), 1)
        self.assertEqual(choices[0]['alternative_disease_interpretations'][0]['disease_id'], 'OMIM:000002')

    def test_html_does_not_interpolate_evidence_as_markup(self):
        rows, evidence = toy()
        rows[0]['flags'] = ['</script><script>alert("synthetic")</script>']
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'review.html'
            package.create_review_html(path, rows, evidence, package.review_manifest(rows, 'toyhash'))
            text = path.read_text()
            self.assertNotIn('</script><script>alert', text)
            self.assertIn("connect-src 'none'", text)
            self.assertIn('textContent', text)
            self.assertNotIn('fetch(', text)
            self.assertNotIn('XMLHttpRequest', text)

    def test_csv_escaping_and_schema(self):
        rows, evidence = toy()
        selection = reviewed(rows)
        selection['choices'][0]['reason'] = 'Reviewed synthetic evidence, with a comma and "quotation".'
        row, choice = package.validate_choices(selection, rows, 'toyhash')[0]
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'toy.csv'
            with path.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=package.FIELDS)
                writer.writeheader()
                writer.writerow(package.csv_row(row, choice, evidence))
            self.assertTrue(package.validate_csv(path)['valid'])


if __name__ == '__main__':
    unittest.main()
