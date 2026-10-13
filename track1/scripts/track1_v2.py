"""Offline, genome-wide disease-aware variant prioritization.

Public reference acquisition is separate. No model is trained and no networking
or submission is performed. Real findings must stay in local output files.
"""
from __future__ import annotations

if __package__:
    from .local_inputs import input_path
else:
    from local_inputs import input_path
import argparse, csv, gzip, hashlib, itertools, json, math, re, socket, sys, time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
import pysam
import pysam.bcftools
from track1_baseline import (PhenotypeMatcher, extract_phenotypes, clinical_support,
    minimal_allele, phase_call, text_field, serial, digest, write_json)

ROOT=Path(__file__).resolve().parents[1]
AR='HP:0000007'; AD='HP:0000006'; XL={'HP:0001417','HP:0001419','HP:0001423'}
CV_FIELDS=['CLNSIG','CLNSIGCONF','CLNREVSTAT','CLNDN','CLNDISDB','MC','AF_EXAC','AF_ESP','AF_TGP','GENEINFO']


def deny_network(*args,**kwargs):
    raise RuntimeError('Network access is disabled in patient analysis')


def install_network_guard():
    socket.create_connection=deny_network
    socket.socket=deny_network


def local_path(value):
    if '://' in str(value):raise ValueError('Only local files are accepted')
    return Path(value)


def safe_key(chrom,pos,ref,alt):
    p,r,a=minimal_allele(pos,ref.upper(),alt.upper())
    return (chrom.removeprefix('chr').replace('MT','M'),p,r,a)


class DiseaseMatcher:
    """Same ontology/IC as baseline, but separate disease profiles and inheritance."""
    def __init__(self,obo,associations,aliases):
        self.base=PhenotypeMatcher(obo,associations,aliases)
        self.profiles=defaultdict(lambda:defaultdict(set))
        with Path(associations).open() as f:
            for r in csv.DictReader(f,delimiter='\t'):
                gene=aliases.get(r['gene_symbol'],r['gene_symbol'])
                if not gene or gene=='-':continue
                self.profiles[gene][r['disease_id']].add(r['hpo_id'])
    def score_terms(self,terms,query):
        q=set(self.base.alts.get(t,t) for t in query)
        q={t for t in q if not any(t!=r and t in self.base.ancestors(r) for r in q)}
        closure=set().union(*(self.base.ancestors(t) for t in terms)) if terms else set()
        numerator=denominator=0.0;details=[]
        for t in sorted(q):
            weight=self.base.ic.get(t,0);common=self.base.ancestors(t)&closure
            match=max(common,key=lambda x:(self.base.ic.get(x,0),x)) if common else None
            support=self.base.ic.get(match,0);numerator+=support;denominator+=weight
            details.append({'term':t,'matched_ancestor':match,'information':weight,'matched_information':support})
        return {'score':numerator/denominator if denominator else 0,'terms':details}
    def gene(self,gene,query,core):
        rows=[]
        for disease,terms in self.profiles.get(gene,{}).items():
            modes=[]
            if AR in terms:modes.append('AR')
            if AD in terms:modes.append('AD')
            if terms&XL:modes.append('XL')
            if not modes:modes=['unknown']
            rows.append({'disease_id':disease,'modes':modes,**self.score_terms(terms,query)})
        if gene in core and not any('AR' in r['modes'] for r in rows):
            terms=self.base.direct.get(gene,set())
            rows.append({'disease_id':'literature_mva_recessive','modes':['AR'],
                'inheritance_source':'literature panel; not inferred from genotype',**self.score_terms(terms,query)})
        if not rows:rows=[{'disease_id':'unannotated','modes':['unknown'],'score':0,'terms':[]}]
        return sorted(rows,key=lambda r:(-r['score'],r['disease_id']))


def read_csq(rec,alt,config):
    """Local csq's DNA change disambiguates consequences at multiallelic sites."""
    rows=[]; target=f'{rec.pos}{rec.ref}>{alt}'
    for text in rec.info.get('BCSQ',[]):
        parts=text.split('|');parts+=['']*(7-len(parts))
        effect,gene,tx,biotype,strand,aa,dna=parts[:7]
        if effect.startswith('@'):continue
        if dna and target not in dna.split('+'):continue
        if not dna and len(rec.alts)>1:continue
        if not gene:continue
        types=effect.lstrip('*').split('&')
        score=max([config['effect_scores'].get(t,0) for t in types],default=0)
        if score==0:continue
        rows.append({'gene':config['aliases'].get(gene,gene),'transcript':tx,'consequences':types,
                     'effect_score':score,'amino_acid_change':aa,'dna_change':dna,'strand':strand,
                     'annotation':'bcftools csq local per record; no haplotype inference'})
    return rows


def extract_variants(path,config):
    variants={};counts=Counter()
    with pysam.VariantFile(str(path)) as f:
        if len(f.header.samples)!=1:raise ValueError('Single explicitly selected sample required')
        sn=list(f.header.samples)[0]
        if 'GRCh38' not in str(f.header):raise ValueError('GRCh38 header required')
        for rec in f:
            counts['records']+=1;s=rec.samples[sn];gt=s.get('GT')
            if not gt or None in gt:continue
            for ai in sorted(set(gt)-{0}):
                alt=rec.alts[ai-1]
                if not set(rec.ref+alt)<=set('ACGT'):continue
                anns=read_csq(rec,alt,config)
                if not anns:continue
                k=safe_key(rec.contig,rec.pos,rec.ref,alt)
                dp,gq,ad=s.get('DP'),s.get('GQ'),s.get('AD')
                ab=ad[ai]/sum(ad) if ad and len(ad)>ai and all(x is not None for x in ad) and sum(ad)>0 else None
                genotype='hom_alt' if len(gt)==2 and gt.count(ai)==2 else 'het' if len(gt)==2 else 'hemizygous'
                flags=[]
                if list(rec.filter)!=['PASS']:flags.append('FILTER_not_PASS')
                if dp is None or dp<config['minimum_dp']:flags.append('low_or_missing_DP')
                if gq is None or gq<config['minimum_gq']:flags.append('low_or_missing_GQ')
                if genotype=='het' and ab is not None and not config['minimum_het_balance']<=ab<=config['maximum_het_balance']:
                    flags.append('heterozygous_allele_balance_outside_0.2_to_0.8')
                phase_set=hap=None
                if genotype=='het':
                    if s.phased and s.get('PS') is not None:
                        phase_set=f'{rec.contig}:PS:{s["PS"]}';hap=gt.index(ai)
                    elif s.get('PID') and re.fullmatch(r'\d+\|\d+',s.get('PGT') or ''):
                        pgt=tuple(map(int,s['PGT'].split('|')))
                        if sorted(pgt)==sorted(gt):phase_set=f'{rec.contig}:PID:{s["PID"]}';hap=pgt.index(ai)
                chrom='chr'+rec.contig.removeprefix('chr').replace('MT','M')
                id=f'{chrom}:{rec.pos}:{rec.ref}:{alt}'
                if k not in variants:
                    variants[k]={'id':id,'chrom':chrom,'source_contig':rec.contig,'pos':rec.pos,'ref':rec.ref,'alt':alt,'alt_index':ai,
                        'record_id':f'{rec.contig}:{rec.pos}:{rec.ref}:{",".join(rec.alts)}','gt':list(gt),'genotype':genotype,
                        'dp':dp,'gq':gq,'ad':serial(ad),'allele_balance':ab,'phase_set':phase_set,'haplotype':hap,
                        'quality_pass':not flags,'quality_flags':flags,'annotations':[],'clinvar_records':[]}
                variants[k]['annotations'].extend(anns)
    counts['functional_alleles']=len(variants)
    return variants,dict(counts)


def join_clinvar(variants,path):
    with pysam.VariantFile(str(path)) as f:
        for r in f:
            for a in r.alts or []:
                k=safe_key(r.contig,r.pos,r.ref,a)
                if k in variants:
                    variants[k]['clinvar_records'].append({'id':r.id,**{x:serial(r.info[x]) for x in CV_FIELDS if x in r.info}})


def load_public_functional(path):
    """Optional schema standardized by the public-only downloader."""
    p=Path(path)/'functional_catalogue.json'
    if not p.exists():raise FileNotFoundError('Prepare the required public functional catalogue before analysis')
    return json.loads(p.read_text())


def available_population(v):
    public=v.get('public_annotation',{});gn=public.get('gnomad',{})
    usable=[];unreliable=[]
    for typ in ['exome','genome']:
        x=gn.get(typ) or {}
        if not x:continue
        if x.get('filters'):unreliable.append(typ);continue
        if x.get('an',0)>0 and x.get('af') is not None:usable.append(float(x['af']))
    if usable:return {'known':True,'af':max(usable),'source':'gnomAD public gene catalogue PASS','unreliable_callsets':unreliable}
    # Legacy ClinVar population fields are intentionally NOT used as positive rarity evidence.
    legacy=[]
    for c in v['clinvar_records']:
        for field in ['AF_EXAC','AF_ESP','AF_TGP']:
            values=c.get(field,[]);values=values if isinstance(values,(list,tuple)) else [values]
            legacy.extend(float(z) for z in values if z is not None and 0<=float(z)<=1)
    return {'known':False,'af':None,'source':'unknown or filtered; no rarity credit',
            'unreliable_callsets':unreliable,'legacy_af_for_common_variant_exclusion_only':max(legacy) if legacy else None}


def evidence(v,ann,disease):
    cs=[clinical_support(c) for c in v['clinvar_records']]
    any_benign=any(x['benign'] for x in cs);any_path=any(x['pathogenic_assertion'] for x in cs)
    conflict=any(x['conflicting'] for x in cs) or (any_benign and any_path)
    disease_id=disease['disease_id'];matched=[]
    for cv,c in zip(v['clinvar_records'],cs):
        condition_ids=set(re.findall(r'(?:OMIM|ORPHA|MONDO):[A-Za-z0-9_]+',text_field(cv.get('CLNDISDB'))))
        if disease_id in condition_ids and c['pathogenic_assertion'] and not conflict:
            matched.append(c['support'])
    clinical=max(matched,default=0)
    score=ann['effect_score'];public=v.get('public_annotation',{})
    possible=public.get('alphamissense_annotations') or ([public['alphamissense']] if public.get('alphamissense') else [])
    am=next((a for a in possible if a.get('gene',ann['gene'])==ann['gene']
             and a.get('transcript_id','').split('.')[0]==ann['transcript'].split('.')[0]),None)
    am_used=False
    if 'missense' in ann['consequences'] and am:
        # Match effect to catalogue protein change. Scores from another isoform must not be transferred.
        aa=ann.get('amino_acid_change','');m=re.fullmatch(r'(\d+)([A-Z*])>\d+([A-Z*])',aa)
        expected=f'{m[2]}{m[1]}{m[3]}' if m else None
        if expected and am.get('protein_variant')==expected:
            value=float(am['am_pathogenicity'])
            if not math.isfinite(value) or not 0<=value<=1:raise ValueError('Invalid public missense score')
            score=.20+.50*value;am_used=True
    score=max(score,clinical)
    flags=[]
    if not v['population']['known']:flags.append('population_frequency_unknown_or_filtered')
    if conflict:flags.append('conflicting_ClinVar_evidence')
    if not clinical:flags.append('no_condition_matched_pathogenic_ClinVar_support')
    if 'missense' in ann['consequences'] and not am_used:flags.append('no_transcript_matched_missense_predictor')
    if not am_used and 'missense' in ann['consequences']:flags.append('missense_function_unconfirmed')
    if am_used:flags.append('AlphaMissense_computational_only_not_functional_validation')
    if am_used and am.get('transcript_id')!=ann['transcript']:
        flags.append('AlphaMissense_transcript_version_differs_genomic_and_protein_change_match')
    return {'score':score,'effect_score':ann['effect_score'],'clinical_support':clinical,
        'has_pathogenic_assertion_any_condition':any_path and not conflict,'benign':any_benign and not any_path and not conflict,
        'clinical_conflict':conflict,'alphamissense_used':am_used,'flags':flags}


def model(variants,profiles,config,*,pheno_weight=None,mechanism_weight=None,no_clinvar=False,no_am=False):
    grouped=defaultdict(list);review=Counter();rows=[]
    for v in variants:
        if not v['quality_pass']:review['quality_excluded']+=1;continue
        best={}
        for a in v['annotations']:
            if a['gene'] not in best or a['effect_score']>best[a['gene']]['effect_score']:best[a['gene']]=a
        for gene,ann in best.items():grouped[gene].append((v,ann))
    pw=config['phenotype_weight'] if pheno_weight is None else pheno_weight
    mw=config['mechanism_weight'] if mechanism_weight is None else mechanism_weight
    for gene,items in grouped.items():
        for disease in profiles[gene]:
            eligible=[]
            for v,ann in items:
                evv=v
                if no_clinvar:evv={**evv,'clinvar_records':[]}
                if no_am:evv={**evv,'public_annotation':{k:v for k,v in evv.get('public_annotation',{}).items()
                                                      if k not in {'alphamissense','alphamissense_annotations'}}}
                ev=evidence(evv,ann,disease)
                if ev['benign']:review['benign_excluded']+=1;continue
                if ev['score']<.3:continue
                eligible.append((v,ann,ev))
            if len(eligible)>config['maximum_eligible_alleles_per_gene']:
                review['gene_disease_pair_enumeration_skipped_excess_alleles']+=1
                # Retain anchored variants only. No arbitrary truncation to top N hides specific alleles.
                eligible=[x for x in eligible if x[2]['clinical_support']>0 or x[1]['effect_score']>=.65]
            for mode in disease['modes']:
                if mode=='unknown':continue
                threshold=config['maximum_population_af_ar'] if mode=='AR' else config['maximum_population_af_ad']
                vs=[]
                for x in eligible:
                    pop=x[0]['population'];af=pop['af'] if pop['known'] else pop.get('legacy_af_for_common_variant_exclusion_only')
                    if af is not None and af>threshold:review['common_variant_excluded']+=1;continue
                    vs.append(x)
                combos=[]
                if mode=='AR':
                    combos.extend(([x],'homozygous','Homozygous called genotype') for x in vs if x[0]['genotype']=='hom_alt')
                    for a,b in itertools.combinations([x for x in vs if x[0]['genotype']=='het'],2):
                        ph,why=phase_call(a[0],b[0])
                        if ph=='cis':review['cis_pair_excluded']+=1;continue
                        combos.append(([a,b],ph,why))
                else:
                    for x in vs:
                        # Dominant/X-linked generic missense predictions alone are insufficient for a primary hypothesis.
                        if x[2]['clinical_support']<=0:continue
                        if mode=='XL' and x[0]['genotype'] not in {'hemizygous','hom_alt'}:continue
                        combos.append(([x],'not_applicable','Single allele; clinical assertion must match disease'))
                for members,phase,why in combos:
                    weakest=min(x[2]['score'] for x in members)
                    strongest=max(x[2]['score'] for x in members)
                    phenotype=disease['score'];mechanism=1.0 if gene in config['mva_recessive_genes'] and mode=='AR' else 0
                    # The strongest allele contributes only a little: two-copy evidence remains the bottleneck.
                    genetic=.85*weakest+.15*strongest
                    score=(1-pw-mw)*genetic+pw*phenotype+mw*mechanism
                    phase_factor=config['unknown_phase_factor'] if phase=='unknown' else 1
                    score*=phase_factor
                    rows.append({'gene':gene,'disease_id':disease['disease_id'],'inheritance':mode,
                        'members':[x[0]['id'] for x in members],'phase':phase,'phase_evidence':why,'score':score,
                        'components':{'weakest_allele':weakest,'strongest_allele':strongest,'phenotype':phenotype,'literature_mechanism':mechanism,'phase_factor':phase_factor},
                        'allele_evidence':[x[2] for x in members],
                        'flags':sorted({f for x in members for f in x[2]['flags']})})
    rows.sort(key=lambda r:(-r['score'],r['gene'],r['members'],r['disease_id']))
    unique=[];seen=set()
    for r in rows:
        k=tuple(sorted(r['members']))
        if k in seen:continue
        seen.add(k);r['rank']=len(unique)+1;unique.append(r)
    return unique,dict(review)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--vcf',type=local_path,default=input_path('vcf'))
    p.add_argument('--annotated-vcf',type=local_path,default=ROOT/'outputs/track1_v2/annotated.bcf')
    p.add_argument('--fasta',type=local_path,default=ROOT/'references/track1_v2/hg38.fa')
    p.add_argument('--gff',type=local_path,required=True)
    p.add_argument('--phenotype',type=local_path,default=input_path('phenotype'))
    p.add_argument('--references',type=local_path,default=ROOT/'references/track1')
    p.add_argument('--functional',type=local_path,default=ROOT/'references/track1_v2_functional')
    p.add_argument('--config',type=local_path,default=ROOT/'config/track1_v2.json')
    p.add_argument('--output',type=local_path,default=ROOT/'outputs/track1_v2')
    p.add_argument('--reuse-annotation',action='store_true')
    args=p.parse_args();install_network_guard();started=time.perf_counter();args.output.mkdir(parents=True,exist_ok=True)
    config=json.loads(args.config.read_text())
    annotation_inputs={'vcf':digest(args.vcf),'fasta':digest(args.fasta),'gff':digest(args.gff),'pysam':pysam.__version__}
    receipt=args.annotated_vcf.with_suffix('.inputs.json')
    if args.reuse_annotation and args.annotated_vcf.exists():
        if not receipt.exists() or json.loads(receipt.read_text())!=annotation_inputs:raise ValueError('Annotation cache provenance mismatch')
    else:
        # -l is intentionally local: no fabricated phase is used to join distant calls.
        # The supplied callset can include decoy/viral contigs absent from UCSC.
        # Analyze the primary human chromosomes; disclose excluded contigs locally.
        primary={str(i) for i in range(1,23)}|{'X','Y','M','MT'}
        with pysam.VariantFile(str(args.vcf)) as source:
            targets=[c for c in source.header.contigs if c.removeprefix('chr') in primary]
            excluded=[c for c in source.header.contigs if c not in targets]
        if not targets:raise ValueError('No primary human contigs in input')
        write_json(args.output/'contig_scope.json',{'included':targets,'excluded':excluded})
        pysam.bcftools.csq('-f',str(args.fasta),'-g',str(args.gff),'-l','-p','a','-B','30',
            '-t',','.join(targets),'-Ob','-o',str(args.annotated_vcf),str(args.vcf),catch_stdout=False)
        write_json(receipt,annotation_inputs)
    variants,counts=extract_variants(args.annotated_vcf,config)
    join_clinvar(variants,args.references/'clinvar.vcf.gz')
    catalogue=load_public_functional(args.functional)
    for k,v in variants.items():
        public=catalogue.get(':'.join(map(str,k)),{})
        v['public_annotation']=public;v['population']=available_population(v)
    phenotypes=extract_phenotypes(args.phenotype);query=[x['hpo_id'] for x in phenotypes if x['subject']=='proband']
    matcher=DiseaseMatcher(args.references/'hp.obo',args.references/'genes_to_phenotype.txt',config['aliases'])
    genes={ann['gene'] for v in variants.values() for ann in v['annotations']}
    profiles={g:matcher.gene(g,query,config['mva_recessive_genes']) for g in sorted(genes)}
    ranked,review=model(list(variants.values()),profiles,config)
    scenarios={}
    for name,kwargs in [('without_literature_prior',{'mechanism_weight':0}),('without_phenotype',{'pheno_weight':0}),
        ('phenotype_0.2',{'pheno_weight':.2}),('phenotype_0.5',{'pheno_weight':.5}),
        ('without_ClinVar',{'no_clinvar':True}),('without_AlphaMissense',{'no_am':True})]:
        rr,_=model(list(variants.values()),profiles,config,**kwargs);scenarios[name]=rr[:20]
    for filename,obj in [('candidate_evidence.json',list(variants.values())),('ranked_hypotheses.json',ranked),
        ('disease_phenotype.json',profiles),('sensitivity.json',scenarios),('review_counts.json',review),('annotation_counts.json',counts)]:
        write_json(args.output/filename,obj)
    manifest={'model':config['model'],'run_utc':datetime.now(timezone.utc).isoformat(),
        'config':config,'input_hashes':annotation_inputs,'phenotype_sha256':digest(args.phenotype),
        'config_sha256':digest(args.config),'code_sha256':digest(Path(__file__)),
        'shared_baseline_sha256':digest(ROOT/'scripts/track1_baseline.py'),
        'reference_hashes':{name:digest(args.references/name) for name in ['clinvar.vcf.gz','hp.obo','genes_to_phenotype.txt']},
        'samtools_bcftools_version':pysam.__samtools_version__,
        'functional_catalogue_sha256':digest(args.functional/'functional_catalogue.json') if catalogue else None,
        'elapsed_seconds':time.perf_counter()-started,'training':'none','network_during_analysis':False,
        'genome_wide_csq':True,'annotation_model':'MANE transcript set; local per-record consequences',
        'score_interpretation':'uncalibrated utility; separate provisional manual EPCR',
        'limitations':['unknown phase','only provided sequence calls; structural and repeat causes not excluded',
         'single-record csq can miss multi-variant codon or frame rescue','ClinVar match minimally trimmed, not repeat-left-aligned',
         'AlphaMissense and gnomAD catalogue limited to downloaded public genes','no clinical classification or diagnostic validation',
         'disease-specific inheritance annotations incomplete; XL heterozygous and mitochondrial hypotheses not ranked',
         'MANE omits some alternative-transcript effects','no inference of mosaic chromosome counts from allele balance',
         'prior public participant candidate exposure; not blinded']}
    write_json(args.output/'run_manifest.json',manifest)
    print('Offline v2 analysis completed; private results saved locally.')

if __name__=='__main__':main()
