import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('submission_validator', ROOT / 'scripts/validate_track1_submission.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def toy(**changes):
    row = dict(zip(m.FIELDS, ['PROBAND01', 'chr1', '123', 'A', 'C',
                             'chr1', '456', 'G', 'T', '.8', 'primary', 'Synthetic test only']))
    row.update(changes)
    return row


class SubmissionValidationTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name) / 'synthetic.csv'

    def write(self, rows):
        with self.path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=m.FIELDS)
            writer.writeheader()
            writer.writerows(rows)

    def test_valid_pair_and_single_with_comma_in_notes(self):
        self.write([toy(), toy(chrom_1='chrX', pos_1='789', chrom_2='', pos_2='', ref_2='', alt_2='',
                               epcr='.2', finding_type='secondary', notes='Quoted, note')])
        result = m.validate_csv(self.path)
        self.assertTrue(result['valid'], result)
        self.assertEqual(result['rows'], 2)
        self.assertTrue(any('secondary' in w for w in result['warnings']))

    def test_incomplete_pair_is_rejected(self):
        self.write([toy(alt_2='')])
        self.assertFalse(m.validate_csv(self.path)['valid'])

    def test_nonfinite_and_out_of_range_confidence(self):
        for value in ['nan', 'inf', '-inf', '0', '-.1', '1.1', '']:
            with self.subTest(value=value):
                self.write([toy(epcr=value)])
                self.assertFalse(m.validate_csv(self.path)['valid'])

    def test_duplicate_reversed_pair_and_identical_members(self):
        self.write([toy(), toy(pos_1='456', ref_1='G', alt_1='T', pos_2='123', ref_2='A', alt_2='C')])
        self.assertTrue(any('Duplicate hypothesis' in e for e in m.validate_csv(self.path)['errors']))
        self.write([toy(pos_2='123', ref_2='A', alt_2='C')])
        self.assertFalse(m.validate_csv(self.path)['valid'])

    def test_wrong_build_style_and_bad_alleles(self):
        for changes in [{'chrom_1': '15'}, {'chrom_1': 'chr23'}, {'pos_1': '0'},
                        {'pos_1': '12.5'}, {'alt_1': '<DEL>'}, {'alt_1': 'C,G'},
                        {'alt_1': 'A'}, {'proband_id': 'OTHER'}, {'finding_type': 'incidental'}]:
            with self.subTest(changes=changes):
                self.write([toy(**changes)])
                self.assertFalse(m.validate_csv(self.path)['valid'])

    def test_row_limits_and_malformed_records(self):
        for rows in [[], [toy(pos_1=str(i + 1)) for i in range(11)]]:
            self.write(rows)
            self.assertFalse(m.validate_csv(self.path)['valid'])
        self.path.write_text(','.join(m.FIELDS) + '\nPROBAND01,chr1,123\n')
        self.assertFalse(m.validate_csv(self.path)['valid'])

    def test_duplicate_headers_rejected(self):
        self.path.write_text(','.join(m.FIELDS) + ',notes\n')
        self.assertFalse(m.validate_csv(self.path)['valid'])

    def test_order_ties_and_bom(self):
        self.write([toy(epcr='.5'), toy(pos_1='124', epcr='.8'), toy(pos_1='125', epcr='.8')])
        self.path.write_bytes(b'\xef\xbb\xbf' + self.path.read_bytes())
        result = m.validate_csv(self.path)
        self.assertTrue(result['valid'])
        self.assertTrue(any('descending' in w for w in result['warnings']))
        self.assertTrue(any('ties' in w for w in result['warnings']))

    def test_package_never_claims_upload_or_verified_repo(self):
        self.write([toy()])
        report = Path(self.directory.name) / 'methods.md'
        report.write_text('# Synthetic methods\nLocal test report.\n')
        result = m.validate_package(self.path, report, 'https://github.com/example/analysis', True)
        self.assertTrue(result['local_validation_passed'])
        self.assertTrue(result['package_files_present'])
        self.assertEqual(result['remote_submission_status'], 'not_submitted_by_this_tool')
        self.assertTrue(any('GitHub' in item for item in result['not_verified']))

    def test_missing_package_placeholder_and_fake_pdf(self):
        self.write([toy()])
        self.assertFalse(m.validate_package(self.path, require_package=True)['local_validation_passed'])
        report = Path(self.directory.name) / 'methods.md'
        report.write_text('AI plan: [TO CONFIRM]\n')
        self.assertFalse(m.validate_package(self.path, report)['local_validation_passed'])
        report = report.with_suffix('.pdf')
        report.write_text('not a PDF')
        self.assertFalse(m.validate_package(self.path, report)['local_validation_passed'])
        self.assertFalse(m.validate_package(self.path, github='https://github.com/your-org/your-repo')['local_validation_passed'])


if __name__ == '__main__':
    unittest.main()
