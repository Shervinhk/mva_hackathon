"""Fabricated fixtures only: never reads a patient callset or derived output."""
import copy
import importlib.util
import json
from pathlib import Path
import sqlite3
import sys
import tempfile
import unittest
from unittest.mock import patch

import pysam

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
spec = importlib.util.spec_from_file_location('v3', ROOT / 'scripts/track1_v3.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def config():
    return json.loads((ROOT / 'config/track1_v3.json').read_text())


def annotation(gene='TOY', transcript='ENSTTOY.1', consequence='missense', score=.35):
    return {'gene': gene, 'transcript': transcript, 'consequences': [consequence],
            'effect_score': score, 'amino_acid_change': '2E>2K', 'dna_change': '4G>A', 'strand': '+'}


def variant(pos=4, genotype='het', chrom='chr1', anns=None):
    v = {'id': f'{chrom}:{pos}:G:A', 'chrom': chrom, 'source_contig': chrom, 'pos': pos, 'ref': 'G', 'alt': 'A',
         'alt_index': 1, 'record_id': f'{chrom}:{pos}:G:A', 'gt': [1, 1] if genotype == 'hom_alt' else [0, 1],
         'genotype': genotype, 'phase_set': None, 'haplotype': None, 'quality_pass': True,
         'quality_flags': [], 'dp': 30, 'gq': 99, 'ad': [15, 15], 'allele_balance': .5,
         'annotations': anns or [annotation()], 'clinvar_records': [], 'public_annotation': {}}
    v['population'] = m.population_evidence(v)
    return v


def disease(mode='AR', score=.8, disease_id='OMIM:999999'):
    return {'disease_id': disease_id, 'modes': [mode], 'score': score, 'terms': []}


def assertion(classification='Pathogenic', disease_id='OMIM:999999', kind='SCV', **extra):
    return {'accession': 'SCV_SYNTHETIC', 'disease_ids': [disease_id], 'assertion_type': kind,
            'classification': classification, 'review_stars': 1, 'review_status': 'criteria provided, single submitter', **extra}


def synthetic_vcf(path, records):
    path.write_text('''##fileformat=VCFv4.2
##reference=GRCh38
##contig=<ID=chr1,length=100>
##INFO=<ID=BCSQ,Number=.,Type=String,Description="Synthetic local consequences">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC
''' + '\n'.join(records) + '\n')


class ReferenceJoinTests(unittest.TestCase):
    def test_repeated_absent_and_family_terms_do_not_inflate_proband_similarity(self):
        present = {'hpo_id': 'HP:0000001', 'subject': 'proband', 'present': True}
        terms = [present, dict(present),
                 {'hpo_id': 'HP:0000002', 'subject': 'proband', 'present': False},
                 {'hpo_id': 'HP:0000003', 'subject': 'family', 'present': True}]
        self.assertEqual(m.positive_proband_terms(terms), ['HP:0000001'])

    def test_left_alignment_and_source_alleles_retained(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'synthetic.fa'
            path.write_text('>chr1\nCAAAAAGT\n')
            pysam.faidx(str(path))
            with pysam.FastaFile(str(path)) as fasta:
                self.assertEqual(m.normalize_allele('1', 4, 'AA', 'A', fasta), ('1', 1, 'CA', 'C'))
                self.assertEqual(m.normalize_allele('chr1', 2, 'AA', 'A', fasta), ('1', 1, 'CA', 'C'))
                v = variant(4)
                v.update(ref='AA', alt='A')
                normalized, _ = m.normalized_variants([v], fasta)
                got = next(iter(normalized.values()))
                self.assertEqual((got['pos'], got['ref'], got['alt']), (1, 'CA', 'C'))
                self.assertEqual((got['source_pos'], got['source_ref'], got['source_alt']), (4, 'AA', 'A'))
                self.assertEqual(got['gt'], [0, 1])
                with self.assertRaises(ValueError):
                    m.normalize_allele('1', 1, 'T', 'A', fasta)

    def test_read_only_sqlite_view_and_exact_allele_key(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'public.sqlite'
            db = sqlite3.connect(path)
            db.execute('CREATE TABLE raw(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,payload TEXT)')
            db.execute('INSERT INTO raw VALUES(?,?,?,?,?)', ('1', 4, 'G', 'A', json.dumps([assertion()])))
            db.execute('CREATE VIEW clinvar_rcv AS SELECT * FROM raw')
            db.execute('CREATE TABLE metadata(key TEXT,value TEXT)')
            db.execute('INSERT INTO metadata VALUES(?,?)', ('release', json.dumps('synthetic')))
            db.commit(); db.close()
            index = m.PublicIndex(path)
            self.assertEqual(index.lookup(('1', 4, 'G', 'A'))['rcv_assertions'][0]['assertion_type'], 'SCV')
            self.assertEqual(index.lookup(('1', 4, 'G', 'T')), {})
            self.assertEqual(index.metadata['release'], 'synthetic')
            with self.assertRaises(sqlite3.OperationalError):
                index.connection.execute('DELETE FROM raw')
            index.close()

    def test_merge_empty_mvstore_channel_does_not_erase_transcript_predictor(self):
        legacy = {'alphamissense_annotations': [{'transcript_id': 'ENSTTOY.1'}]}
        mvstore = {'alphamissense_annotations': [], 'other_predictors': {'ALPHA_MISSENSE': {'score': .8}}}
        result = m.merge_annotations(legacy, mvstore)
        self.assertEqual(len(result['alphamissense_annotations']), 1)
        self.assertIn('ALPHA_MISSENSE', result['other_predictors'])

    def test_original_equivalent_records_are_preserved_before_genotype_conflict_check(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'synthetic.vcf'
            synthetic_vcf(path, [
                'chr1\t4\t.\tG\tA\t.\tPASS\tBCSQ=missense|TOY|ENSTTOY.1|protein_coding|+|2E>2K|4G>A\tGT:DP:GQ:AD\t0/1:30:99:15,15',
                'chr1\t3\t.\tAG\tAA\t.\tPASS\tBCSQ=missense|TOY|ENSTTOY.1|protein_coding|+|2E>2K|3AG>AA\tGT:DP:GQ:AD\t1/1:30:99:0,30'])
            variants, _ = m.extract_variants(path, config())
            self.assertEqual(len(variants), 2)
            normalized, counts = m.normalized_variants(variants)
            self.assertEqual(len(normalized), 1)
            allele = next(iter(normalized.values()))
            self.assertEqual([r['source_pos'] for r in allele['source_alleles']], [4, 3])
            self.assertEqual([r['gt'] for r in allele['source_alleles']], [[0, 1], [1, 1]])
            self.assertFalse(allele['quality_pass'])
            self.assertIn('equivalent_records_genotype_disagreement', allele['quality_flags'])
            self.assertEqual(counts['equivalent_source_alleles_merged'], 1)

    def test_equivalent_alt_indexes_use_dosage_and_quality_flags_are_combined(self):
        a, b = variant(), variant()
        b.update(gt=[0, 2], alt_index=2, record_id='chr1:4:G:T,A')
        normalized, _ = m.normalized_variants([a, b])
        self.assertTrue(next(iter(normalized.values()))['quality_pass'])
        b.update(quality_pass=False, quality_flags=['low_or_missing_GQ'])
        normalized, _ = m.normalized_variants([a, b])
        got = next(iter(normalized.values()))
        self.assertFalse(got['quality_pass'])
        self.assertIn('low_or_missing_GQ', got['quality_flags'])

    def test_conflicting_duplicate_phase_is_not_silently_promoted_to_trans(self):
        a, duplicate, b = variant(4), variant(4), variant(8)
        a.update(phase_set='1:PS:1', haplotype=0)
        duplicate.update(phase_set='1:PS:1', haplotype=1)
        b.update(phase_set='1:PS:1', haplotype=0)
        normalized, _ = m.normalized_variants([a, duplicate, b])
        phase, reason = m.phase_call(*normalized.values())
        self.assertEqual(phase, 'unknown')
        self.assertIn('Conflicting', reason)

    def test_incomplete_public_build_is_rejected_but_population_placeholder_is_allowed(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'public.sqlite'
            db = sqlite3.connect(path)
            db.execute('CREATE TABLE metadata(key TEXT,value TEXT)')
            db.execute('INSERT INTO metadata VALUES(?,?)', ('population', '{"complete":false}'))
            db.commit()
            m.PublicIndex(path).close()
            db.execute('INSERT INTO metadata VALUES(?,?)', ('clinvar_scv', '{"complete":false}'))
            db.commit(); db.close()
            with self.assertRaisesRegex(ValueError, 'incomplete'):
                m.PublicIndex(path)


class PopulationTests(unittest.TestCase):
    def test_unknown_is_not_zero_frequency_and_zero_requires_denominator(self):
        v = variant()
        self.assertFalse(v['population']['known'])
        self.assertIsNone(v['population']['af'])
        v['public_annotation'] = {'population': [{'af': 0, 'filters': []}]}
        self.assertFalse(m.population_evidence(v)['known'])
        v['public_annotation']['population'][0]['an'] = 1000
        measured = m.population_evidence(v)
        self.assertTrue(measured['known'])
        self.assertEqual(measured['af'], 0)
        self.assertEqual(measured['rarity_credit'], 0)

    def test_filtered_reference_rejected_and_source_quality_not_invented(self):
        v = variant()
        v['public_annotation'] = {'population': [{'af': .3, 'an': 100, 'filters': ['FAIL'], 'source': 'bad'},
            {'af': .002, 'an': None, 'filters': None, 'source': 'Exomiser:PUBLIC', 'quality_scope': 'source preprocessing'}]}
        result = m.population_evidence(v)
        self.assertTrue(result['known'])
        self.assertEqual(result['af'], .002)
        self.assertFalse(result['all_measurement_quality_verified'])
        self.assertEqual(result['unreliable_callsets'], ['bad'])

    def test_group_frequency_can_exclude_observed_common_allele(self):
        v = variant(genotype='hom_alt')
        v['public_annotation'] = {'population': [{'af': .0001, 'an': 100000, 'filters': [], 'groupmax_af': .1}]}
        v['population'] = m.population_evidence(v)
        rows, counts = m.model([v], {'TOY': [disease()]}, config())
        self.assertEqual(rows, [])
        self.assertEqual(counts['observed_common_alleles_excluded'], 1)

    def test_invalid_values_do_not_become_observed_absence_or_crash(self):
        v = variant()
        for row in [{'af': 'missing'}, {'af': 0, 'an': float('nan')}, {'AF': 0, 'AN': '0'}, {'ac': 0, 'an': -1}]:
            v['public_annotation'] = {'population': [row]}
            self.assertFalse(m.population_evidence(v)['known'])
        v['public_annotation'] = {'population': [{'AC': 2, 'AN': 1000, 'groupmax_af': 'missing',
            'populations': [{'ac': 200, 'an': 100}, {'ac': 'unknown', 'an': 100}]}]}
        self.assertEqual(m.population_evidence(v)['af'], .002)

    def test_population_missingness_changes_flags_not_utility_or_tie_order(self):
        missing = variant(4, genotype='hom_alt')
        measured = variant(8, genotype='hom_alt')
        measured['public_annotation'] = {'population': [{'af': 0, 'an': 10000, 'filters': []}]}
        measured['population'] = m.population_evidence(measured)
        rows, _ = m.model([measured, missing], {'TOY': [disease()]}, config())
        self.assertEqual(rows[0]['score'], rows[1]['score'])
        self.assertEqual(rows[0]['members'], [missing['id']])


class ConditionAssertionTests(unittest.TestCase):
    def test_aggregate_label_not_condition_specific_evidence(self):
        v = variant()
        v['clinvar_records'] = [{'CLNSIG': 'Pathogenic', 'CLNDISDB': 'OMIM:999999', 'CLNREVSTAT': 'reviewed_by_expert_panel'}]
        e = m.evidence(v, v['annotations'][0], disease(), config())
        self.assertEqual(e['clinical_support'], 0)
        self.assertIn('aggregate_ClinVar_retained_for_review_not_condition_specific_support', e['flags'])

    def test_scv_condition_match_and_wrong_condition_do_not_share_support(self):
        v = variant()
        v['public_annotation'] = {'rcv_assertions': [assertion()]}
        yes = m.rcv_evidence(v, disease())
        no = m.rcv_evidence(v, disease(disease_id='OMIM:888888'))
        self.assertEqual(yes['support'], .75)
        self.assertEqual(yes['source'], 'condition-specific SCV')
        self.assertEqual(no['support'], 0)

    def test_opposing_condition_specific_assertions_and_noncontributors(self):
        v = variant()
        v['public_annotation'] = {'rcv_assertions': [assertion(), assertion('Uncertain significance')]}
        self.assertTrue(m.rcv_evidence(v, disease())['conflict'])
        self.assertEqual(m.rcv_evidence(v, disease())['support'], 0)
        v['public_annotation']['rcv_assertions'][1]['contributes_to_aggregate'] = False
        self.assertEqual(m.rcv_evidence(v, disease())['support'], .75)
        v['public_annotation']['rcv_assertions'].append(assertion('Benign', disease_id='OMIM:888888'))
        self.assertEqual(m.rcv_evidence(v, disease())['support'], .75)

    def test_verified_condition_aliases_and_identifier_spellings(self):
        public = {'rcv_assertions': [assertion(disease_id='MEDGEN:C123456')]}
        public = m.attach_condition_aliases(public, {'MEDGEN:C123456': ['OMIM:999999']})
        v = variant(); v['public_annotation'] = public
        self.assertEqual(m.rcv_evidence(v, disease())['support'], .75)
        self.assertEqual(m.condition_ids('MONDO:MONDO:0000001,Orphanet:123'), {'MONDO:0000001', 'ORPHA:123'})

    def test_classification_is_exact_and_submission_context_is_retained(self):
        v = variant()
        v['public_annotation'] = {'rcv_assertions': [assertion('not pathogenic')]}
        self.assertEqual(m.rcv_evidence(v, disease())['support'], 0)
        v['public_annotation'] = {'rcv_assertions': [assertion(classification_context='somatic clinical impact')]}
        self.assertEqual(m.rcv_evidence(v, disease())['assertions'], [])
        submitted = assertion('Pathogenic/Likely pathogenic', submitter='Fabricated lab',
            origin_counts='germline:1', collection_method='clinical testing', date_last_evaluated='2020-01-01')
        v['public_annotation'] = {'rcv_assertions': [submitted]}
        result = m.rcv_evidence(v, disease())
        self.assertEqual(result['support'], .75)
        for key in ['submitter', 'origin_counts', 'collection_method', 'date_last_evaluated']:
            self.assertEqual(result['assertions'][0][key], submitted[key])


class PredictorAndMechanismTests(unittest.TestCase):
    def test_transcript_evidence_is_evaluated_before_selection(self):
        v = variant(anns=[annotation(transcript='ENSTA.1'), annotation(transcript='ENSTB.1')])
        v['public_annotation'] = {'alphamissense_annotations': [{'gene': 'TOY', 'transcript_id': 'ENSTB.1',
            'protein_variant': 'E2K', 'am_pathogenicity': .9}]}
        ann, e = m.best_transcript_evidence(v, 'TOY', disease(), 'AR', config(), predictor_weight=1)
        self.assertEqual(ann['transcript'], 'ENSTB.1')
        self.assertTrue(e['alphamissense_used'])
        self.assertAlmostEqual(e['score'], .65)

    def test_missing_predictor_has_no_default_score_penalty_or_credit(self):
        a = variant(); b = copy.deepcopy(a)
        b['public_annotation'] = {'alphamissense_annotations': [{'gene': 'TOY', 'transcript_id': 'ENSTTOY.1',
            'protein_variant': 'E2K', 'am_pathogenicity': .02}]}
        ea = m.evidence(a, a['annotations'][0], disease(), config())
        eb = m.evidence(b, b['annotations'][0], disease(), config())
        self.assertEqual(ea['score'], eb['score'])
        self.assertFalse(ea['predictor_evidence']['available'])
        self.assertTrue(eb['predictor_evidence']['available'])
        b['public_annotation']['alphamissense_annotations'][0]['protein_variant'] = 'E2R'
        self.assertFalse(m.predictor_evidence(b, b['annotations'][0])['available'])

    def test_genomic_score_is_not_fabricated_transcript_annotation_and_splice_rescue_is_scoped(self):
        v = variant(anns=[annotation(consequence='splice_region', score=.2)])
        v['public_annotation'] = {'other_predictors': {'ALPHA_MISSENSE': {'score': .9, 'transcript_unresolved': True},
                                                       'SPLICE_AI': {'score': .9, 'transcript_unresolved': True}}}
        ev = m.evidence(v, v['annotations'][0], disease(), config())
        self.assertFalse(ev['alphamissense_used'])
        self.assertTrue(ev['splice_prediction_rescue'])
        self.assertEqual(ev['score'], .4)
        self.assertIn('splice_prediction_rescue_requires_transcript_and_RNA_review', ev['flags'])

    def test_curated_mechanism_applies_only_to_exact_disease_and_does_not_promote_CDC20_truncation(self):
        cfg = config()
        a = annotation(gene='BUB1B')
        evidence = m.mechanism_evidence(a, disease(disease_id='OMIM:257300'), 'AR', cfg)
        self.assertEqual(evidence['status'], 'plausible')
        self.assertTrue(evidence['source'].startswith('https://'))
        self.assertEqual(m.mechanism_evidence(a, disease(disease_id='OMIM:114500'), 'AD', cfg)['status'], 'unreviewed')
        a = annotation(gene='CDC20', consequence='stop_gained', score=.8)
        evidence = m.mechanism_evidence(a, disease('AD', disease_id='LITERATURE:PMID33094908'), 'AD', cfg)
        self.assertEqual(evidence['status'], 'uncertain')
        self.assertIn('heterozygous_knockout_not_equivalent_to_reported_missense_mechanism', evidence['flags'])


class HypothesisTests(unittest.TestCase):
    def test_unknown_phase_retained_without_fabricated_probability(self):
        rows, _ = m.model([variant(4), variant(8)], {'TOY': [disease()]}, config())
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row['phase'], 'unknown')
        self.assertEqual(row['components']['phase_factor'], 1)
        self.assertIn('phase_unresolved_not_probability', row['flags'])
        self.assertIsNone(row['confidence'])
        self.assertFalse(row['submission_eligible'])
        self.assertTrue(row['eligible_for_manual_selection'])

    def test_known_cis_excluded_and_distinct_multiallelic_genotype_is_trans(self):
        a, b = variant(4), variant(8)
        a.update(phase_set='1:PS:100', haplotype=0)
        b.update(phase_set='1:PS:100', haplotype=0)
        rows, counts = m.model([a, b], {'TOY': [disease()]}, config())
        self.assertEqual(rows, [])
        self.assertEqual(counts['cis_pairs_excluded'], 1)
        a, b = variant(4), variant(4)
        a.update(gt=[1, 2], record_id='1:4:G:A,T')
        b.update(gt=[1, 2], alt_index=2, record_id='1:4:G:A,T', id='chr1:4:G:T', alt='T')
        rows, _ = m.model([a, b], {'TOY': [disease()]}, config())
        self.assertEqual(rows[0]['phase'], 'trans')

    def test_novel_dominant_and_heterozygous_X_linked_are_reviewed_without_sex_or_de_novo_inference(self):
        for mode, chrom in [('AD', 'chr1'), ('XL', 'chrX')]:
            rows, _ = m.model([variant(chrom=chrom)], {'TOY': [disease(mode)]}, config())
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]['review_tier'], 'variant_mechanism_unresolved_review')
            self.assertIn('sex and penetrance are not inferred', rows[0]['phase_evidence'])
            self.assertIn('novel_single_allele_hypothesis_requires_variant_specific_validation', rows[0]['flags'])

    def test_incomplete_recessive_and_unknown_inheritance_go_to_review_only(self):
        rows, _ = m.model([variant()], {'TOY': [disease(), disease('unknown', disease_id='OMIM:888888')]}, config())
        ranked, review = m.split_review_hypotheses(rows)
        self.assertEqual(ranked, [])
        self.assertEqual(len(review), 2)
        self.assertEqual({r['review_tier'] for r in review}, {'carrier_or_unresolved_second_allele', 'inheritance_unresolved_review'})

    def test_large_pair_set_preserves_all_alleles_for_review(self):
        cfg = config(); cfg['maximum_pair_alleles_per_gene_disease'] = 2
        rows, counts = m.model([variant(4), variant(8), variant(12)], {'TOY': [disease()]}, cfg)
        ranked, review = m.split_review_hypotheses(rows)
        self.assertEqual(ranked, [])
        self.assertEqual(len(review), 3)
        self.assertEqual(counts['pair_enumeration_deferred_excess_alleles'], 1)

    def test_evidence_tier_does_not_override_phenotype_utility(self):
        a = variant(4, genotype='hom_alt'); b = variant(8, genotype='hom_alt', anns=[annotation(gene='OTHER')])
        a['public_annotation'] = {'rcv_assertions': [assertion()]}
        profiles = {'TOY': [disease(score=0)], 'OTHER': [disease(score=1)]}
        rows, _ = m.model([a, b], profiles, config())
        self.assertEqual(rows[0]['gene'], 'OTHER')
        self.assertEqual(rows[1]['review_tier'], 'condition_supported_review')

    def test_overlapping_gene_hypotheses_are_not_deduplicated_across_genes(self):
        v = variant(genotype='hom_alt', anns=[annotation(), annotation(gene='OTHER')])
        rows, _ = m.model([v], {'TOY': [disease()], 'OTHER': [disease()]}, config())
        self.assertEqual({r['gene'] for r in rows}, {'TOY', 'OTHER'})


if __name__ == '__main__':
    unittest.main()
