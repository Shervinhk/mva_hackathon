"""Synthetic-only regressions for the genome-wide v2 pipeline."""
import json,sys,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from track1_v2 import read_csq,evidence,available_population,model,extract_variants
from prepare_track1_v2_references import convert_mane
import pysam,pysam.bcftools,gzip
CONFIG=json.loads((Path(__file__).resolve().parents[1]/'config/track1_v2.json').read_text())


def variant():
 return {'id':'chr1:4:G:A','chrom':'chr1','pos':4,'ref':'G','alt':'A','alt_index':1,
  'record_id':'one','gt':[0,1],'genotype':'het','quality_pass':True,'phase_set':'one',
  'haplotype':0,'clinvar_records':[],'public_annotation':{},'population':{'known':False,'af':None},
  'annotations':[{'gene':'SYN','effect_score':.35,'consequences':['missense'],'transcript':'ENST00001.1','amino_acid_change':'2A>2T'}]}

class V2(unittest.TestCase):
 def test_source_contig_is_preserved_and_export_uses_chr_prefix(self):
  with tempfile.TemporaryDirectory() as td:
   path=Path(td)/'synthetic.vcf'
   path.write_text('##fileformat=VCFv4.2\n##reference=GRCh38\n##contig=<ID=1,length=100>\n'
    '##INFO=<ID=BCSQ,Number=.,Type=String,Description="Consequences">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYN\n'
    '1\t4\t.\tG\tT\t100\tPASS\tBCSQ=stop_gained|SYN|TX|protein_coding|+|2E>2*|4G>T\tGT:DP:GQ\t0/1:30:99\n')
   rows,_=extract_variants(path,CONFIG);v=next(iter(rows.values()))
   self.assertEqual(v['chrom'],'chr1');self.assertEqual(v['source_contig'],'1')
   self.assertEqual(v['id'],'chr1:4:G:T')
 def test_multiallelic_effect_separation(self):
  r=SimpleNamespace(pos=4,ref='G',alts=('A','T'),info={'BCSQ':('missense|SYN|ENST00001.1|protein_coding|+|2A>2T|4G>A','stop_gained|SYN|ENST00001.1|protein_coding|+|2A>2*|4G>T')})
  self.assertEqual(read_csq(r,'A',CONFIG)[0]['consequences'],['missense'])
  self.assertEqual(read_csq(r,'T',CONFIG)[0]['consequences'],['stop_gained'])
 def test_multiallelic_missing_dna_not_guessed(self):
  r=SimpleNamespace(pos=4,ref='G',alts=('A','T'),info={'BCSQ':('splice_donor|SYN|TX|protein_coding',)})
  self.assertEqual(read_csq(r,'A',CONFIG),[])
 def test_exact_condition_identifier(self):
  v=variant();v['clinvar_records']=[{'CLNSIG':'Pathogenic','CLNREVSTAT':'criteria_provided,_single_submitter','CLNDISDB':'OMIM:123456'}]
  self.assertEqual(evidence(v,v['annotations'][0],{'disease_id':'OMIM:12345'})['clinical_support'],0)
  self.assertGreater(evidence(v,v['annotations'][0],{'disease_id':'OMIM:123456'})['clinical_support'],0)
 def test_missing_frequency_has_no_invented_af(self):
  p=available_population(variant());self.assertFalse(p['known']);self.assertIsNone(p['af'])
 def test_filtered_population_not_used(self):
  v=variant();v['public_annotation']={'gnomad':{'exome':{'ac':5,'an':100,'af':.05,'filters':['AC0']}}}
  self.assertFalse(available_population(v)['known'])
 def test_missense_wrong_isoform_not_used(self):
  v=variant();v['public_annotation']={'alphamissense':{'gene':'SYN','transcript_id':'ENST00002.1','protein_variant':'A2T','am_pathogenicity':.999}}
  self.assertFalse(evidence(v,v['annotations'][0],{'disease_id':'none'})['alphamissense_used'])
 def test_missense_transcript_version_difference_flagged(self):
  v=variant();v['public_annotation']={'alphamissense':{'gene':'SYN','transcript_id':'ENST00001.2','protein_variant':'A2T','am_pathogenicity':.999}}
  r=evidence(v,v['annotations'][0],{'disease_id':'none'})
  self.assertTrue(r['alphamissense_used']);self.assertTrue(any('version_differs' in f for f in r['flags']))
 def test_cis_is_rejected_trans_is_ranked(self):
  a=variant();b=variant();b.update(id='chr1:7:G:T',pos=7,record_id='two')
  profiles={'SYN':[{'disease_id':'SYN:1','modes':['AR'],'score':.5}]}
  self.assertEqual(model([a,b],profiles,CONFIG)[0],[])
  b['haplotype']=1
  self.assertEqual(model([a,b],profiles,CONFIG)[0][0]['phase'],'trans')
 def test_csq_and_gff_conversion_with_synthetic_sequence(self):
  with tempfile.TemporaryDirectory() as td:
   d=Path(td);fa=d/'ref.fa';fa.write_text('>chr1\nATGGAATTTTAA'+'A'*88+'\n');pysam.faidx(str(fa))
   gff=d/'input.gff.gz'
   with gzip.open(gff,'wt') as f:
    f.write('##gff-version 3\nchr1\tsynthetic\tgene\t1\t12\t.\t+\t.\tID=ENSG1;gene_type=protein_coding;gene_name=SYN\n')
    f.write('chr1\tsynthetic\ttranscript\t1\t12\t.\t+\t.\tID=ENST1;Parent=ENSG1;gene_type=protein_coding;transcript_type=protein_coding;gene_name=SYN\n')
    f.write('chr1\tsynthetic\tCDS\t1\t12\t.\t+\t0\tParent=ENST1;gene_type=protein_coding;transcript_type=protein_coding\n')
   outg=d/'out.gff';convert_mane(gff,outg)
   vc=d/'input.vcf';vc.write_text('##fileformat=VCFv4.2\n##reference=GRCh38\n##contig=<ID=chr1,length=100>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYN\nchr1\t4\t.\tG\tT\t100\tPASS\t.\tGT\t0/1\n')
   ov=d/'out.vcf';pysam.bcftools.csq('-f',str(fa),'-g',str(outg),'-l','-p','a','-Ov','-o',str(ov),str(vc),catch_stdout=False)
   with pysam.VariantFile(str(ov)) as f:
    rec=next(f);annotations=read_csq(rec,'T',CONFIG)
    self.assertEqual(annotations[0]['gene'],'SYN');self.assertEqual(annotations[0]['consequences'],['stop_gained'])

if __name__=='__main__':unittest.main()
