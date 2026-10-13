"""Fabricated document/ontology regressions; no challenge data is read."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from xml.sax.saxutils import escape
from zipfile import ZipFile

from track1_v4.phenotypes import DiseaseMatcherV4, extract_observations, load_observations, positive_terms


def docx(path, rows=(), before='', after='', raw=''):
    def paragraph(value):
        return '<w:p><w:r><w:t>' + escape(value) + '</w:t></w:r></w:p>'
    table = '<w:tbl>' + ''.join('<w:tr>' + ''.join('<w:tc>' + paragraph(cell) + '</w:tc>' for cell in row) + '</w:tr>' for row in rows) + '</w:tbl>'
    xml = '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:body>' + paragraph(before) + table + raw + paragraph(after) + '</w:body></w:document>'
    with ZipFile(path, 'w') as archive:
        archive.writestr('word/document.xml', xml)
    return path


class ObservationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / 'fabricated.docx'

    def test_present_absent_paternal_and_sibling(self):
        docx(self.path, [('HP:0000001', 'Seizure'), ('HP:0000002', 'Absent microcephaly'),
                         ('HP:0000003', 'Paternal finding'), ('HP:0000004', 'Sibling finding')])
        result = extract_observations(self.path)
        self.assertEqual(positive_terms(result), ['HP:0000001'])
        self.assertIs(result['observations'][1]['present'], False)
        self.assertEqual([row['subject'] for row in result['observations'][2:]], ['family', 'family'])
        self.assertIs(result['reviewed_locally'], False)
        self.assertEqual(result['source_sha256'], hashlib.sha256(self.path.read_bytes()).hexdigest())

    def test_formatting_split_hpo_id(self):
        raw = '<w:tbl><w:tr><w:tc><w:p><w:r><w:t>HP:</w:t></w:r><w:r><w:t>1234567</w:t></w:r></w:p></w:tc><w:tc><w:p><w:r><w:t>Seizure</w:t></w:r></w:p></w:tc></w:tr></w:tbl>'
        docx(self.path, raw=raw)
        self.assertEqual(positive_terms(extract_observations(self.path)), ['HP:1234567'])

    def test_ambiguous_scope_and_presence(self):
        docx(self.path, [('HP:0000001', 'Absent in father; present in child'),
                         ('HP:0000002', 'Seizures cannot be excluded'),
                         ('HP:0000003', 'Unknown status'), ('HP:0000004', 'Possible finding')])
        result = extract_observations(self.path)
        self.assertEqual(positive_terms(result), [])
        self.assertEqual(result['observations'][0]['subject'], 'unknown')
        self.assertTrue(all(row['present'] is None for row in result['observations']))

    def test_negated_present_and_absence_seizures(self):
        docx(self.path, [('HP:0000001', 'Seizures not present'), ('HP:0000002', 'Absence seizures')])
        result = extract_observations(self.path)
        self.assertIs(result['observations'][0]['present'], False)
        self.assertEqual(positive_terms(result), ['HP:0000002'])

    def test_no_status_and_not_documented(self):
        docx(self.path, [('HP:0000001', 'Seizures', 'No'), ('HP:0000002', 'Finding not documented')])
        result = extract_observations(self.path)
        self.assertIs(result['observations'][0]['present'], False)
        self.assertIsNone(result['observations'][1]['present'])
        self.assertEqual(positive_terms(result), [])

    def test_alternate_id_contradictions_do_not_become_positive(self):
        value = [{'hpo_id': 'HP:9999991', 'present': True, 'subject': 'proband'},
                 {'hpo_id': 'HP:1000001', 'present': False, 'subject': 'proband'}]
        self.assertEqual(positive_terms(value, {'HP:9999991': 'HP:1000001'}), [])

    def test_family_heading_inherits_scope(self):
        docx(self.path, [('HP:0000001', 'Seizures')], before='Family history', after='Clinical findings')
        self.assertEqual(positive_terms(extract_observations(self.path)), [])

    def test_mixed_subject_headers_are_unresolved(self):
        docx(self.path, [('HPO ID', 'Patient', 'Father'), ('HP:0000001', 'Yes', 'No')])
        result = extract_observations(self.path)
        self.assertEqual(positive_terms(result), [])
        self.assertEqual(result['observations'][0]['subject'], 'unknown')

    def test_single_family_column_is_not_proband(self):
        docx(self.path, [('HPO ID', 'Symptom', 'Father'), ('HP:0000001', 'Seizure', 'Yes')])
        result = extract_observations(self.path)
        self.assertEqual(positive_terms(result), [])
        self.assertEqual(result['observations'][0]['subject'], 'family')

    def test_unscoped_prose_and_example_are_not_observations(self):
        docx(self.path, [('HP:0000001', 'Example HPO identifier')], after='See also HP:0000002')
        result = extract_observations(self.path)
        self.assertEqual(positive_terms(result), [])

    def test_explicit_proband_prose_is_retained(self):
        docx(self.path, after='Patient has seizures HP:0000001')
        self.assertEqual(positive_terms(extract_observations(self.path)), ['HP:0000001'])

    def test_contradictory_rows_and_multiple_qualified_terms(self):
        docx(self.path, [('HP:0000001', 'Seizure present'), ('HP:0000001', 'Seizure absent'),
                         ('HP:0000002 and HP:0000003', 'Absent finding')])
        result = extract_observations(self.path)
        self.assertEqual(positive_terms(result), [])
        self.assertEqual(result['contradictory_terms'], ['HP:0000001'])

    def test_onset_and_source_row(self):
        docx(self.path, [('HP:0000001', 'Seizure; onset: infancy')])
        observation = extract_observations(self.path)['observations'][0]
        self.assertEqual(observation['onset'], 'infancy')
        self.assertEqual(observation['source_row'], 1)

    def override(self, **updates):
        docx(self.path, [('HP:0000001', 'Seizure')])
        value = {'schema_version': 1, 'source_sha256': hashlib.sha256(self.path.read_bytes()).hexdigest(),
                 'reviewed_locally': False, 'observations': [{'hpo_id': 'HP:0000001', 'subject': 'proband',
                 'present': False, 'source_row': 1}]}
        value.update(updates)
        override = self.path.with_suffix('.json')
        override.write_text(json.dumps(value))
        return override

    def test_explicit_override_has_no_assumed_review(self):
        result = load_observations(self.path, self.override())
        self.assertFalse(result['reviewed_locally'])
        self.assertEqual(positive_terms(result), [])
        self.assertIn('override_sha256', result)

    def test_invalid_override_hash_and_review_status(self):
        for updates in ({'source_sha256': '0' * 64}, {'reviewed_locally': 'yes'}, {'schema_version': 2}):
            with self.subTest(updates=updates), self.assertRaises(ValueError):
                load_observations(self.path, self.override(**updates))

    def test_invalid_observation_shapes(self):
        original = {'hpo_id': 'HP:0000001', 'subject': 'proband', 'present': True, 'source_row': 1}
        for updates in ({'present': 1}, {'subject': 'child'}, {'hpo_id': 'HP:abc'}, {'source_row': True}, {'onset': []}):
            with self.subTest(updates=updates), self.assertRaises(ValueError):
                load_observations(self.path, self.override(observations=[{**original, **updates}]))

    def test_empty_document_fails_explicitly(self):
        docx(self.path, [('No identifiers', 'No observations')])
        with self.assertRaises(ValueError):
            extract_observations(self.path)


class SimilarityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        root = Path(self.temp.name)
        obo = root / 'toy.obo'
        obo.write_text('''[Term]
id: HP:1000000
name: Root

[Term]
id: HP:1000001
name: Parent
is_a: HP:1000000
alt_id: HP:9999991

[Term]
id: HP:1000002
name: Unannotated child
is_a: HP:1000001

[Term]
id: HP:1000003
name: Other symptom
is_a: HP:1000000

[Term]
id: HP:0001419
name: X-linked recessive inheritance

[Term]
id: HP:0001423
name: X-linked dominant inheritance

[Term]
id: HP:0001417
name: X-linked inheritance
''')
        associations = root / 'associations.tsv'
        associations.write_text('gene_symbol\thpo_id\tdisease_id\nA\tHP:9999991\tOMIM:1\nA\tHP:0001419\tOMIM:1\nA\tHP:0001417\tOMIM:1\nB\tHP:1000003\tOMIM:2\nB\tHP:0001423\tOMIM:2\nC\tHP:1000000\tOMIM:3\nC\tHP:0001417\tOMIM:3\n')
        self.matcher = DiseaseMatcherV4(obo, associations, {'ALIAS_A': 'A'})

    def test_zero_ic_does_not_add_numerator(self):
        result = self.matcher.score_terms({'HP:1000001', 'HP:1000003'}, {'HP:1000002', 'HP:1000003'})
        self.assertEqual(result['score'], 1.0)
        self.assertEqual(result['unscored_terms'], ['HP:1000002'])
        unsupported = next(row for row in result['terms'] if row['term'] == 'HP:1000002')
        self.assertEqual(unsupported['matched_information'], 0)
        self.assertFalse(unsupported['scored'])

    def test_unscored_child_does_not_remove_supported_parent(self):
        result = self.matcher.score_terms({'HP:1000001'}, {'HP:1000001', 'HP:1000002'})
        self.assertEqual(result['score'], 1.0)
        self.assertEqual(result['scored_terms'], ['HP:1000001'])

    def test_aliases_canonical_for_query_disease_and_gene(self):
        direct = self.matcher.score_terms({'HP:9999991'}, {'HP:1000001'})
        alternate = self.matcher.score_terms({'HP:1000001'}, {'HP:9999991'})
        self.assertEqual(direct['score'], alternate['score'])
        self.assertEqual(alternate['score'], 1.0)
        self.assertEqual(self.matcher.gene('ALIAS_A', ['HP:1000001'], [])[0]['score'], 1.0)

    def test_inheritance_subtypes_are_preserved(self):
        self.assertEqual(self.matcher.gene('A', [], [])[0]['modes'], ['XLR'])
        self.assertEqual(self.matcher.gene('B', [], [])[0]['modes'], ['XLD'])
        self.assertEqual(self.matcher.gene('C', [], [])[0]['modes'], ['XL'])

    def test_unknown_and_zero_information_are_explicit(self):
        result = self.matcher.score_terms({'HP:1000001'}, {'HP:8888888', 'HP:1000000'})
        self.assertEqual(result['score'], 0.0)
        self.assertEqual(result['unscored_terms'], ['HP:1000000', 'HP:8888888'])

    def test_inheritance_term_does_not_become_symptom(self):
        result = self.matcher.score_terms({'HP:0001419'}, {'HP:0001419'})
        self.assertEqual(result['score'], 0)
        self.assertEqual(result['unscored_terms'], ['HP:0001419'])


if __name__ == '__main__':
    unittest.main()
