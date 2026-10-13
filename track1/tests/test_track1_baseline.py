import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import zipfile

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('track1_baseline',ROOT/'scripts/track1_baseline.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def transcript(strand='+'):
    dna='ATGGAAGGTTAA'
    if strand=='-':dna=m.reverse_complement(dna)
    meta={'gene':'TOY','transcript':'TOY.1','strand':strand,'cds':[[0,12]],
          'exons':[[0,12]],'reference_start':0}
    return m.Transcript(meta,dna)


def variant(name,pos=4,gt=None,consequence='missense_variant',clinvar=None):
    return {'id':name,'gene':'TOY','record_id':f'1:{pos}:G:A','alt_index':1,
            'gt':gt or [0,1],'genotype':'hom_alt' if gt==[1,1] else 'het',
            'phase_set':None,'haplotype':None,'quality_pass':True,'quality_flags':[],
            'consequence':consequence,'clinvar':clinvar or {}}


class TranscriptTests(unittest.TestCase):
    def test_positive_strand_stop_synonymous_and_start(self):
        t=transcript()
        self.assertEqual(t.annotate(4,'G','T')['protein_change'],'p.E2*')
        self.assertEqual(t.annotate(6,'A','G')['consequence'],'synonymous_variant')
        self.assertEqual(t.annotate(1,'A','G')['consequence'],'start_lost')
        self.assertEqual(t.annotate(10,'T','C')['consequence'],'stop_lost')

    def test_negative_strand_complements_and_reverses_coordinates(self):
        t=transcript('-')
        self.assertEqual(t.cds,'ATGGAAGGTTAA')
        self.assertEqual(t.annotate(9,'C','A')['protein_change'],'p.E2*')
        self.assertEqual(t.annotate(9,'C','A')['coding_change'],'c.4G>T')

    def test_reference_mismatch_fails_annotation(self):
        self.assertEqual(transcript().annotate(4,'A','T')['consequence'],'reference_mismatch')

    def test_frameshift_vs_inframe(self):
        t=transcript()
        self.assertEqual(t.annotate(4,'G','GA')['consequence'],'frameshift_variant')
        self.assertEqual(t.annotate(4,'G','GAAA')['consequence'],'inframe_indel')

    def test_splice_boundary_is_not_any_intron(self):
        meta={'gene':'TOY','transcript':'TOY.1','strand':'+','cds':[[0,6],[16,22]],
              'exons':[[0,6],[16,22]],'reference_start':0}
        t=m.Transcript(meta,'ATGGAA'+'GTAAAAAAAG'+'GGTTAA')
        self.assertEqual(t.annotate(7,'G','A')['consequence'],'splice_donor_variant')
        self.assertEqual(t.annotate(16,'G','A')['consequence'],'splice_acceptor_variant')
        self.assertNotEqual(t.annotate(11,'A','T')['consequence'],'splice_donor_variant')

    def test_invalid_reference_protein_rejected(self):
        with self.assertRaises(ValueError):
            m.Transcript({'gene':'TOY','transcript':'T','strand':'+','cds':[[0,9]],
                          'exons':[[0,9]],'reference_start':0},'ATGTAATAA')


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.config=json.loads((ROOT/'config/track1_panel.json').read_text())
        self.pheno={'TOY':{'score':.8}}

    def test_unknown_frequency_not_rare_and_patient_AF_ignored(self):
        r=m.rarity_evidence({'AF':0.5})
        self.assertFalse(r['known']);self.assertEqual(r['score'],.5)
        self.assertGreater(m.rarity_evidence({'AF_EXAC':[0.000001]})['score'],r['score'])

    def test_conflicting_clinvar_not_pathogenic_support(self):
        c=m.clinical_support({'CLNSIG':['Conflicting_classifications_of_pathogenicity'],
                              'CLNSIGCONF':['Pathogenic(1)','Benign(1)']})
        self.assertTrue(c['conflicting']);self.assertEqual(c['support'],0)

    def test_clinvar_review_stars_require_exact_criteria_token(self):
        statuses=[(['no_assertion_criteria_provided'],0),
                  (['no_assertion_provided'],0),
                  (['criteria_provided','_single_submitter'],1),
                  (['criteria_provided','_multiple_submitters','_no_conflicts'],2),
                  (['reviewed_by_expert_panel'],3),
                  (['practice_guideline'],4)]
        for review,stars in statuses:
            with self.subTest(review=review):
                result=m.clinical_support({'CLNSIG':['Pathogenic'],'CLNREVSTAT':review})
                self.assertEqual(result['review_stars'],stars)
                self.assertEqual(result['support'],[.65,.75,.85,.95,.95][stars])

    def test_recessive_pair_required_and_unknown_phase_flagged(self):
        a=variant('a');b=variant('b',8)
        self.assertEqual(m.generate_hypotheses([a],self.pheno,self.config)[0],[])
        r=m.generate_hypotheses([a,b],self.pheno,self.config)[0][0]
        self.assertEqual(r['phase'],'unknown')

    def test_phase_blocks_and_multiallelic_phase(self):
        a=variant('a');b=variant('b',8)
        a.update(phase_set='PS1',haplotype=0);b.update(phase_set='PS1',haplotype=0)
        self.assertEqual(m.phase_call(a,b)[0],'cis')
        self.assertEqual(m.generate_hypotheses([a,b],self.pheno,self.config)[0],[])
        b['haplotype']=1;self.assertEqual(m.phase_call(a,b)[0],'trans')
        b['phase_set']='PS2';self.assertEqual(m.phase_call(a,b)[0],'unknown')
        a=variant('a',gt=[1,2]);b=variant('b',gt=[1,2]);b['alt_index']=2
        self.assertEqual(m.phase_call(a,b)[0],'trans')

    def test_qc_and_common_benign_variants_are_reviewed(self):
        a=variant('a');b=variant('b',8)
        b['quality_pass']=False;b['quality_flags']=['GQ_missing_or_low']
        self.assertFalse(m.generate_hypotheses([a,b],self.pheno,self.config)[0])
        b=variant('b',8,clinvar={'CLNSIG':['Benign'],'AF_EXAC':[.1]})
        rows,review=m.generate_hypotheses([a,b],self.pheno,self.config)
        self.assertFalse(rows);self.assertIn('population_af_above_threshold',review[0]['reasons'])

    def test_weaker_allele_limits_pair_and_inputs_are_deterministic(self):
        a=variant('a',consequence='stop_gained');b=variant('b',8)
        r=m.generate_hypotheses([a,b],self.pheno,self.config)[0][0]
        self.assertAlmostEqual(r['components']['weakest_allele_score'],b['features']['score'])


class PhenotypeTests(unittest.TestCase):
    def test_family_hpo_is_separated(self):
        xml='''<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:body><w:tbl>
        <w:tr><w:tc><w:p><w:r><w:t>Short stature</w:t></w:r></w:p></w:tc><w:tc><w:p><w:r><w:t>HP:0004322</w:t></w:r></w:p></w:tc></w:tr>
        <w:tr><w:tc><w:p><w:r><w:t>Parental miscarriages</w:t></w:r></w:p></w:tc><w:tc><w:p><w:r><w:t>HP:0200067</w:t></w:r></w:p></w:tc></w:tr>
        </w:tbl></w:body></w:document>'''
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'toy.docx'
            with zipfile.ZipFile(path,'w') as z:z.writestr('word/document.xml',xml)
            p=m.extract_phenotypes(path)
            self.assertEqual([i['subject'] for i in p],['proband','family'])

    def test_ontology_ancestors_and_redundant_terms(self):
        with tempfile.TemporaryDirectory() as d:
            obo=Path(d)/'hp.obo';assoc=Path(d)/'genes.txt'
            obo.write_text('[Term]\nid: HP:0000001\nname: Root\n\n[Term]\nid: HP:0000002\nname: Specific\nis_a: HP:0000001 ! Root\n\n[Term]\nid: HP:0000003\nname: Other\nis_a: HP:0000001 ! Root\n')
            assoc.write_text('gene_symbol\thpo_id\nA\tHP:0000002\nB\tHP:0000003\n')
            matcher=m.PhenotypeMatcher(obo,assoc,{})
            good=matcher.score('A',['HP:0000001','HP:0000002'])
            bad=matcher.score('B',['HP:0000002'])
            self.assertEqual(good['removed_redundant_terms'],['HP:0000001'])
            self.assertAlmostEqual(good['score'],1)
            self.assertGreater(good['score'],bad['score'])


if __name__=='__main__':unittest.main()
