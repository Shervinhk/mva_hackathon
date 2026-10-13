"""Fabricated entry-preparation fixtures only; never open private analyses."""
import copy
import csv
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import prepare_track1_v4_entry as prep
import build_track1_v3_submission as v3package
from track1_baseline import digest
from track1_automated_export import automated_choices
from validate_track1_submission import FIELDS


def supplied_identity(**details):
    return {'team': 'Synthetic team', 'github_url': 'https://github.com/example/synthetic-v4',
            'ai_disclosure': '', 'ai_details': details}


def toy_prediction():
    members = ['chr1:10:A:C', 'chr1:20:G:T']
    row = {'rank': 1, 'score': .71, 'gene': 'SYNTHETIC', 'disease_id': 'OMIM:000001',
           'inheritance': 'AR', 'phase': 'unknown', 'members': members,
           'eligible_for_manual_selection': True, 'automatic_shortlist_eligible': True}
    evidence = {members[0]: {'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'C'},
                members[1]: {'chrom': 'chr1', 'pos': 20, 'ref': 'G', 'alt': 'T'}}
    selected, _ = automated_choices([row], 'fabricated-manifest', v3package.hypothesis_id)
    data = io.StringIO(newline='')
    writer = csv.DictWriter(data, fieldnames=FIELDS)
    writer.writeheader()
    for hypothesis, choice in selected:
        choice = dict(choice, reason='Automated V4 order; EPCR is an ordinal weight, not a causal probability; no human review claimed.')
        writer.writerow(v3package.csv_row(hypothesis, choice, evidence))
    return data.getvalue().encode(), selected, evidence


class ProvenanceTests(unittest.TestCase):
    def test_current_inventory_passes_and_stale_source_fails(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / 'source.py').write_text('# Fabricated source')
            hashes = {'source.py': digest(root / 'source.py')}
            prep.check_hashes(root, hashes)
            (root / 'source.py').write_text('# Different fabricated source')
            with self.assertRaises(ValueError):
                prep.check_hashes(root, hashes)

    def test_empty_inventory_missing_files_and_traversal_fail_closed(self):
        with tempfile.TemporaryDirectory() as temp:
            top = Path(temp); root = top / 'root'; root.mkdir()
            outside = top / 'outside.py'; outside.write_text('# Synthetic outside fixture')
            for hashes in ({}, None, {'missing.py': 'synthetic'},
                           {'../outside.py': digest(outside)}, {str(outside): digest(outside)}):
                with self.subTest(hashes=hashes), self.assertRaises(ValueError):
                    prep.check_hashes(root, hashes)

    def test_symlink_to_external_source_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            top = Path(temp); root = top / 'root'; root.mkdir()
            outside = top / 'outside.py'; outside.write_text('# Synthetic outside fixture')
            (root / 'source.py').symlink_to(outside)
            with self.assertRaises(ValueError):
                prep.check_hashes(root, {'source.py': digest(outside)})

    def test_symlink_to_internal_source_is_also_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            real = root / 'real.py'; real.write_text('# Synthetic source')
            (root / 'link.py').symlink_to(real)
            with self.assertRaises(ValueError):
                prep.check_hashes(root, {'link.py': digest(real)})


class IdentityTests(unittest.TestCase):
    def test_actual_plan_and_explicit_other_tools_are_required(self):
        identity, missing = prep.entry_identity(supplied_identity())
        self.assertIn('ai_details.plan_tier', missing)
        self.assertIn('ai_details.other_ai_tools', missing)
        self.assertIn('ai_disclosure', missing)
        self.assertEqual(identity['ai_details']['plan_tier'], '')
        self.assertIn('not independently inspected', identity['ai_details']['data_controls'])

    def test_freeform_old_disclosure_does_not_replace_missing_structured_facts(self):
        supplied = supplied_identity()
        supplied['ai_disclosure'] = 'Fabricated old statement with no plan or other tools.'
        identity, missing = prep.entry_identity(supplied)
        self.assertEqual(identity['ai_disclosure'], supplied['ai_disclosure'])
        self.assertIn('ai_details.plan_tier', missing)
        self.assertIn('ai_details.other_ai_tools', missing)

    def test_complete_participant_facts_generate_disclosure_without_verification_claim(self):
        identity, missing = prep.entry_identity(supplied_identity(plan_tier='Synthetic test tier', other_ai_tools='None'))
        self.assertFalse(missing)
        self.assertIn('Plan/tier: Synthetic test tier.', identity['ai_disclosure'])
        self.assertIn('Other AI tools: None.', identity['ai_disclosure'])
        self.assertIn('account controls were not independently inspected', identity['ai_disclosure'])
        self.assertEqual(identity['ai_details']['provider'], 'OpenAI')

    def test_unknown_or_placeholder_facts_are_not_complete(self):
        for key in ('plan_tier', 'other_ai_tools'):
            for value in ('', 'TODO', '[TO CONFIRM]', 'Unknown', 'not sure'):
                with self.subTest(key=key, value=value):
                    details = {'plan_tier': 'Synthetic test tier', 'other_ai_tools': 'None'}
                    details[key] = value
                    _, missing = prep.entry_identity(supplied_identity(**details))
                    self.assertIn('ai_details.' + key, missing)

    def test_nontext_ai_details_rejected(self):
        for supplied in (dict(supplied_identity(), ai_details=[]), supplied_identity(plan_tier=True)):
            with self.subTest(supplied=supplied), self.assertRaises(ValueError):
                prep.entry_identity(supplied)


class PredictionBindingTests(unittest.TestCase):
    def test_exact_existing_pair_csv_is_preserved_byte_for_byte(self):
        contents, selected, evidence = toy_prediction()
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'predictions.csv'; path.write_bytes(contents)
            self.assertEqual(prep.checked_prediction_bytes(path, selected, evidence, digest(path)), contents)

    def test_hash_or_model_order_change_is_rejected_even_with_valid_schema(self):
        contents, selected, evidence = toy_prediction()
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'predictions.csv'; path.write_bytes(contents)
            with self.assertRaises(ValueError):
                prep.checked_prediction_bytes(path, selected, evidence, 'stale-recorded-hash')
            changed = copy.deepcopy(selected)
            changed[0][1]['epcr'] = .4
            with self.assertRaises(ValueError):
                prep.checked_prediction_bytes(path, changed, evidence, digest(path))
            changed_evidence = copy.deepcopy(evidence)
            changed_evidence['chr1:10:A:C']['pos'] = 11
            with self.assertRaises(ValueError):
                prep.checked_prediction_bytes(path, selected, changed_evidence, digest(path))

    def test_schema_failure_rejected_even_if_hash_matches(self):
        _, selected, evidence = toy_prediction()
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'predictions.csv'; path.write_text('wrong,header\n1,2\n')
            with self.assertRaises(ValueError):
                prep.checked_prediction_bytes(path, selected, evidence, digest(path))


class MethodsTests(unittest.TestCase):
    def test_replace_identity_and_insert_rationale_without_changing_scientific_body(self):
        original = ('# Original\n\n- Team: Old synthetic identity\n\n## Abstract\n'
                    'Preserved fabricated abstract.\n\n## Track 1 method\nPreserved method.\n\n'
                    '## Validation and limitations\nPreserved limitations.\n')
        identity, _ = prep.entry_identity(supplied_identity(plan_tier='Synthetic tier', other_ai_tools='None'))
        report = prep.prepared_methods(original, identity, '## Candidate rationale\nFabricated rationale.')
        self.assertNotIn('Old synthetic identity', report)
        self.assertIn('Synthetic team', report)
        self.assertEqual(report.count('## Abstract'), 1)
        self.assertEqual(report.count('## Candidate rationale'), 1)
        self.assertIn('Preserved fabricated abstract.', report)
        self.assertIn('Preserved method.', report)
        self.assertIn('Preserved limitations.', report)
        self.assertLess(report.index('## Candidate rationale'), report.index('## Validation and limitations'))

    def test_unexpected_methods_structure_fails_closed(self):
        for original in ('Missing abstract entirely', '## Abstract\nMissing validation section'):
            with self.subTest(original=original), self.assertRaises(ValueError):
                prep.prepared_methods(original, supplied_identity(), 'Synthetic rationale')


class SourceBundleTests(unittest.TestCase):
    def fixture(self, root):
        for name in ('scripts', 'tests', 'config'):
            (root / name).mkdir()
        for name in ('README.md', 'run.py', 'config.json'):
            path = root / 'track1_v4' / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('Fabricated public ' + name)
        path = root / 'docs/TRACK1_V4_ENTRY_PREPARATION.md'
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('Fabricated preparation instructions.')
        # This is an excluded synthetic canary, never actual local metadata.
        (root / 'track1_v4/entry_identity.json').write_text('SYNTHETIC_EXCLUDED_CANARY')

    def fake_archive(self, stage, extra=False, traversal=False):
        repo = stage / 'repository'
        (repo / 'docs').mkdir(parents=True)
        (repo / 'README.md').write_text('Fabricated V3 readme.')
        (repo / 'docs/public.md').write_text('Fabricated V3 public documentation.')
        if extra:
            (repo / 'unlisted.json').write_text('SYNTHETIC_STALE_CANARY')
        with zipfile.ZipFile(stage / 'repository.zip', 'w') as archive:
            for relative in ('README.md', 'docs/public.md'):
                archive.write(repo / relative, 'mva-track1-v3/' + relative)
            if traversal:
                archive.writestr('mva-track1-v3/../../outside.py', 'Untrusted public test entry')

    def test_only_explicit_fresh_source_is_packaged_with_hash_inventory(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / 'source'; root.mkdir(); self.fixture(root)
            stage = Path(temp) / 'stage'; stage.mkdir()
            files = ('README.md', 'run.py', 'config.json')
            with mock.patch.object(prep, 'ROOT', root), mock.patch.object(prep, 'V4_FILES', files), \
                    mock.patch.object(v3package, 'source_archive', side_effect=self.fake_archive):
                inventory = prep.source_bundle(stage)
            expected = {'README.md', 'docs/public.md', 'docs/TRACK1_V4_ENTRY_PREPARATION.md',
                        'track1_v4/README.md', 'track1_v4/run.py', 'track1_v4/config.json'}
            self.assertEqual(set(inventory), expected)
            with zipfile.ZipFile(stage / 'repository.zip') as archive:
                self.assertEqual(set(archive.namelist()), {'mva-track1-v4/' + name for name in expected})
                self.assertFalse(any(b'SYNTHETIC_EXCLUDED_CANARY' in archive.read(name) for name in archive.namelist()))
            self.assertTrue(all(digest(stage / 'repository' / name) == value for name, value in inventory.items()))

    def test_preexisting_source_directory_is_rejected_before_copy(self):
        with tempfile.TemporaryDirectory() as temp:
            stage = Path(temp); (stage / 'repository').mkdir()
            with mock.patch.object(v3package, 'source_archive') as archive, self.assertRaises(ValueError):
                prep.source_bundle(stage)
            archive.assert_not_called()

    def test_unlisted_staging_files_and_symlinked_v4_sources_are_rejected(self):
        for mutation in ('unlisted', 'symlink'):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as temp:
                root = Path(temp) / 'source'; root.mkdir(); self.fixture(root)
                stage = Path(temp) / 'stage'; stage.mkdir()
                if mutation == 'symlink':
                    path = root / 'track1_v4/run.py'; path.unlink()
                    path.symlink_to(root / 'track1_v4/README.md')
                with mock.patch.object(prep, 'ROOT', root), mock.patch.object(prep, 'V4_FILES', ('run.py',)), \
                        mock.patch.object(v3package, 'source_archive', side_effect=lambda p: self.fake_archive(p, extra=mutation == 'unlisted')):
                    with self.assertRaises(ValueError):
                        prep.source_bundle(stage)

    def test_archive_allowlist_path_traversal_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / 'source'; root.mkdir(); self.fixture(root)
            stage = Path(temp) / 'stage'; stage.mkdir()
            (Path(temp) / 'outside.py').write_text('SYNTHETIC_OUTSIDE_CANARY')
            with mock.patch.object(prep, 'ROOT', root), mock.patch.object(prep, 'V4_FILES', ('run.py',)), \
                    mock.patch.object(v3package, 'source_archive', side_effect=lambda p: self.fake_archive(p, traversal=True)):
                with self.assertRaises(ValueError):
                    prep.source_bundle(stage)
            # A rejected path must not first be written into a shareable ZIP.
            with zipfile.ZipFile(stage / 'repository.zip') as archive:
                self.assertFalse(any(b'SYNTHETIC_OUTSIDE_CANARY' in archive.read(name) for name in archive.namelist()))


if __name__ == '__main__':
    unittest.main()
