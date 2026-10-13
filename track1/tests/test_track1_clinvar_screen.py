import importlib.util
from pathlib import Path
import tempfile
import unittest

import pysam

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('screen', ROOT / 'scripts/track1_clinvar_screen.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class ClinVarIntersectionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)

    def records(self, lines, sample=True):
        path = Path(self.temporary.name) / ('patient.vcf' if sample else 'reference.vcf')
        header = '##fileformat=VCFv4.2\n##reference=GRCh38\n##contig=<ID=1,length=10000>\n'
        if sample:
            header += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' + ('\tFORMAT\tPATIENT' if sample else '') + '\n'
        path.write_text(header + '\n'.join(lines) + '\n')
        with pysam.VariantFile(path) as f:
            return list(f)

    def intersection(self, patient, reference):
        return list(m.intersect_alleles(m.normalized_alleles(patient, 'PATIENT'),
                                       m.normalized_alleles(reference)))

    def test_exact_allele_and_not_same_position_different_alt(self):
        patient = self.records(['1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/1'])
        wrong = self.records(['1\t10\twrong\tA\tT\t.\t.\t.'], False)
        self.assertEqual(self.intersection(patient, wrong), [])
        right = self.records(['1\t10\tright\tA\tG\t.\t.\t.'], False)
        self.assertEqual(self.intersection(patient, right)[0][0], (10, 'A', 'G'))

    def test_only_called_multiallelic_alt_is_matched(self):
        patient = self.records(['1\t10\t.\tA\tG,T\t.\tPASS\t.\tGT\t0/2'])
        reference = self.records(['1\t10\tr1\tA\tG\t.\t.\t.', '1\t10\tr2\tA\tT\t.\t.\t.'], False)
        hits = self.intersection(patient, reference)
        self.assertEqual([x[0] for x in hits], [(10, 'A', 'T')])
        self.assertEqual(hits[0][1][0][1], 2)

    def test_prefix_trimming_reordering_preserves_all_matches(self):
        patient = self.records(['1\t10\t.\tAAAA\tAAAG\t.\tPASS\t.\tGT\t0/1',
                                '1\t11\t.\tA\tT\t.\tPASS\t.\tGT\t0/1'])
        reference = self.records(['1\t11\tr1\tA\tT\t.\t.\t.', '1\t13\tr2\tA\tG\t.\t.\t.'], False)
        self.assertEqual([x[0] for x in self.intersection(patient, reference)], [(11, 'A', 'T'), (13, 'A', 'G')])

    def test_ref_mismatch_does_not_match_and_missing_gt_is_excluded(self):
        patient = self.records(['1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t./1'])
        reference = self.records(['1\t10\tr1\tA\tG\t.\t.\t.'], False)
        self.assertEqual(self.intersection(patient, reference), [])
        patient = self.records(['1\t10\t.\tC\tG\t.\tPASS\t.\tGT\t0/1'])
        self.assertEqual(self.intersection(patient, reference), [])

    def test_duplicate_reference_rows_are_preserved(self):
        patient = self.records(['1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t1/1'])
        reference = self.records(['1\t10\tr1\tA\tG\t.\t.\t.', '1\t10\tr2\tA\tG\t.\t.\t.'], False)
        self.assertEqual(len(self.intersection(patient, reference)[0][2]), 2)

    def test_contig_aliases(self):
        self.assertEqual(m.canonical_contig('chr1'), '1')
        self.assertEqual(m.canonical_contig('chrM'), 'MT')

    def test_conflict_label_is_not_itself_a_pathogenic_assertion(self):
        self.assertFalse(m.has_pathogenic_conflict({
            'CLNSIG': ['Conflicting_classifications_of_pathogenicity'],
            'CLNSIGCONF': ['Likely_benign(1)', 'Uncertain_significance(1)']}))
        self.assertTrue(m.has_pathogenic_conflict({'CLNSIGCONF': ['Pathogenic(1)', 'Benign(1)']}))
        self.assertTrue(m.has_pathogenic_conflict({'CLNSIGCONF': ['Likely_pathogenic(1)', 'Uncertain_significance(1)']}))

    def test_exact_review_status_tokens(self):
        statuses = [(['no_assertion_criteria_provided'], 0),
                    (['criteria_provided', '_single_submitter'], 1),
                    (['criteria_provided', '_multiple_submitters', '_no_conflicts'], 2),
                    (['reviewed_by_expert_panel'], 3),
                    (['practice_guideline'], 4)]
        for status, expected in statuses:
            self.assertEqual(m.clinical_support({'CLNSIG': ['Pathogenic'], 'CLNREVSTAT': status})['review_stars'], expected)

    def test_full_screen_preserves_quality_tiers_and_exhausts_input_counts(self):
        directory = Path(self.temporary.name)
        patient = directory / 'synthetic_patient.vcf'
        reference = directory / 'synthetic_clinvar.vcf'
        patient.write_text('''##fileformat=VCFv4.2
##reference=GRCh38
##contig=<ID=chr1,length=10000>
##contig=<ID=chr2,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPATIENT
chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT:DP:GQ:AD\t0/1:20:50:10,10
chr1\t11\t.\tC\tT\t.\tPASS\t.\tGT:DP:GQ:AD\t0/1:5:10:2,3
chr1\t12\t.\tA\tT\t.\tPASS\t.\tGT:DP:GQ:AD\t1/1:20:50:0,20
chr2\t99\t.\tC\tG\t.\tPASS\t.\tGT:DP:GQ:AD\t0/1:20:50:10,10
''')
        reference.write_text('''##fileformat=VCFv4.2
##contig=<ID=1,length=10000>
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="Review">
##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Gene">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t10\t100\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter;GENEINFO=TOY:999
1\t11\t101\tC\tT\t.\t.\tCLNSIG=Likely_pathogenic;CLNREVSTAT=no_assertion_criteria_provided;GENEINFO=TOY:999
1\t12\t102\tA\tT\t.\t.\tCLNSIG=Uncertain_significance;CLNREVSTAT=criteria_provided,_single_submitter;GENEINFO=TOY:999
''')
        for path in [patient, reference]:
            pysam.tabix_compress(str(path), str(path) + '.gz', force=True)
            pysam.tabix_index(str(path) + '.gz', preset='vcf', force=True)
        result = m.screen(str(patient) + '.gz', str(reference) + '.gz')
        self.assertEqual(result['patient_counts']['records'], 4)
        self.assertEqual(result['clinvar_counts_in_intersecting_contigs']['records'], 3)
        self.assertEqual(result['match_counts']['distinct_minimal_alleles_matched'], 3)
        self.assertEqual(len(result['prioritized_hits']), 2)
        self.assertTrue(result['prioritized_hits'][0]['quality_pass'])
        self.assertEqual(result['prioritized_hits'][0]['clinical_priority_class'], 'unconflicted_PLP_at_least_one_star')
        self.assertFalse(result['prioritized_hits'][1]['quality_pass'])
        self.assertEqual(result['prioritized_hits'][1]['clinical_priority_class'], 'unconflicted_PLP_zero_star')
        self.assertEqual(result['patient_contigs_absent_from_clinvar'], [{'contig': 'chr2', 'records': 1}])


if __name__ == '__main__':
    unittest.main()
