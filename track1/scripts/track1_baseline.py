"""Offline Track 1 panel baseline with local sequence annotation and HPO matching.

Research ranking only. No training, network requests, answer key, or submissions.
"""
from __future__ import annotations

if __package__:
    from .local_inputs import input_path
else:
    from local_inputs import input_path
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import statistics
import time
from zipfile import ZipFile
import xml.etree.ElementTree as ET
import pysam

ROOT = Path(__file__).resolve().parents[1]
COMP = str.maketrans('ACGTNacgtn','TGCANtgcan')
CODONS = dict(zip((''.join(p) for p in itertools.product('TCAG', repeat=3)),
                 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))


def digest(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()


def write_json(path, obj):
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')


def reverse_complement(seq):
    return seq.translate(COMP)[::-1].upper()


def minimal_allele(pos, ref, alt):
    """Trim redundant bases; does NOT left-align repeats. pos is 1-based."""
    while len(ref)>1 and len(alt)>1 and ref[-1]==alt[-1]:
        ref,alt=ref[:-1],alt[:-1]
    while len(ref)>1 and len(alt)>1 and ref[0]==alt[0]:
        ref,alt,pos=ref[1:],alt[1:],pos+1
    return pos,ref,alt


def changed_interval(pos, ref, alt):
    """Zero-based half-open affected REF interval; insertion is empty."""
    start=pos-1
    while ref and alt and ref[0]==alt[0]:
        start+=1;ref=ref[1:];alt=alt[1:]
    while ref and alt and ref[-1]==alt[-1]:
        ref=ref[:-1];alt=alt[:-1]
    return start,start+len(ref),ref,alt


class Transcript:
    def __init__(self, meta, dna):
        self.meta=meta;self.dna=dna.upper()
        intervals=sorted(meta['cds'],reverse=meta['strand']=='-')
        coords=[]
        for s,e in intervals:
            coords.extend(range(s,e) if meta['strand']=='+' else range(e-1,s-1,-1))
        self.coding_index={g:i for i,g in enumerate(coords)}
        self.cds=''.join(self.dna[g-meta['reference_start']] for g in coords)
        if meta['strand']=='-':
            self.cds=self.cds.translate(COMP)
        self.protein=''.join(CODONS.get(self.cds[i:i+3],'X') for i in range(0,len(self.cds)-2,3))
        if len(self.cds)%3 or not self.cds.startswith('ATG') or not self.protein.endswith('*') or '*' in self.protein[:-1]:
            raise ValueError(f"Reference CDS validation failed: {meta['gene']} {meta['transcript']}")

    def annotate(self,pos,ref,alt):
        m=self.meta;start=pos-1
        observed=self.dna[start-m['reference_start']:start-m['reference_start']+len(ref)]
        if observed!=ref:
            return {'consequence':'reference_mismatch','expected_ref':observed,'transcript':m['transcript']}
        s,e,changed_ref,changed_alt=changed_interval(pos,ref,alt)
        if s==e:
            impacted={s-1,s}
        else:
            impacted=set(range(s,e))
        exons=sorted(m['exons'])
        splice=None
        for (_,left_end),(right_start,_) in zip(exons,exons[1:]):
            donor=set(range(left_end,left_end+2))
            acceptor=set(range(right_start-2,right_start))
            if m['strand']=='-':
                donor,acceptor=acceptor,donor
            if impacted & donor:
                splice='splice_donor_variant'
            if impacted & acceptor:
                splice='splice_acceptor_variant'
        coding=[self.coding_index[g] for g in range(s,e) if g in self.coding_index]
        result={'transcript':m['transcript'],'consequence':'intron_variant',
                'protein_change':None,'coding_change':None}
        if splice:
            result['consequence']=splice
            return result
        if len(ref)==len(alt)==1 and start in self.coding_index:
            idx=self.coding_index[start];offset=idx%3;codon_start=idx-offset
            codon=self.cds[codon_start:codon_start+3]
            oriented_alt=alt if m['strand']=='+' else reverse_complement(alt)
            oriented_ref=ref if m['strand']=='+' else reverse_complement(ref)
            mutant=codon[:offset]+oriented_alt+codon[offset+1:]
            aa,bb=CODONS[codon],CODONS[mutant]
            consequence='synonymous_variant' if aa==bb else 'missense_variant'
            if aa!='*' and bb=='*':consequence='stop_gained'
            elif aa=='*' and bb!='*':consequence='stop_lost'
            elif idx<3 and aa=='M' and bb!='M':consequence='start_lost'
            result.update(consequence=consequence,protein_change=f'p.{aa}{idx//3+1}{bb}',
                          coding_change=f'c.{idx+1}{oriented_ref}>{oriented_alt}')
            return result
        insertion_coding=s==e and all(g in self.coding_index for g in [s-1,s]) and abs(self.coding_index[s]-self.coding_index[s-1])==1
        full_coding=e>s and len(coding)==e-s and max(coding)-min(coding)+1==len(coding)
        if insertion_coding or full_coding:
            delta=len(changed_alt)-len(changed_ref)
            if delta:
                result['consequence']='frameshift_variant' if delta%3 else 'inframe_indel'
            else:
                result['consequence']='coding_complex'
            result['annotation_caveat']='Simple CDS-overlap consequence; no HGVS indel normalization, NMD or protein reconstruction'
            return result
        if coding:
            result['consequence']='coding_complex'
            result['annotation_caveat']='Partial CDS overlap; requires established annotator review'
            return result
        for (_,left_end),(right_start,_) in zip(exons,exons[1:]):
            region=set(range(left_end-3,left_end+8)) | set(range(right_start-8,right_start+3))
            if impacted & region:
                result['consequence']='splice_region_variant'
                return result
        if any(impacted & set(range(a,b)) for a,b in exons):
            result['consequence']='UTR_variant'
        return result


def extract_phenotypes(docx):
    with ZipFile(docx) as z:
        root=ET.fromstring(z.read('word/document.xml'))
    ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
    observations=[]
    for row in root.findall('.//w:tr',ns):
        cells=[' '.join(t.text or '' for t in c.findall('.//w:t',ns)) for c in row.findall('./w:tc',ns)]
        text=' '.join(cells)
        ids=re.findall(r'HP:\d{7}',text)
        for term in ids:
            family=bool(re.search(r'parent|maternal|family history',text,re.I))
            observations.append({'hpo_id':term,'label':cells[1] if len(cells)>1 else '',
                                 'subject':'family' if family else 'proband','present':True,
                                 'source':'Supplied clinical document table'})
    if not observations:
        raise ValueError('No HPO terms found in the clinical table; curate manually')
    return observations


class PhenotypeMatcher:
    def __init__(self, obo, associations, aliases):
        parents=defaultdict(set);labels={};alts={};term=None
        for line in Path(obo).read_text().splitlines():
            if line=='[Term]':term=None
            elif line.startswith('id: HP:'):term=line[4:]
            elif term and line.startswith('name: '):labels[term]=line[6:]
            elif term and line.startswith('is_a: '):parents[term].add(line[6:].split()[0])
            elif term and line.startswith('alt_id: '):alts[line[8:]]=term
        self.labels=labels;self.parents=parents;self.alts=alts
        direct=defaultdict(set)
        with Path(associations).open() as f:
            header=None
            for line in f:
                if not line.strip():continue
                parts=line.lstrip('#').rstrip('\n').split('\t')
                if 'gene_symbol' in parts and 'hpo_id' in parts:
                    header=parts;continue
                if line.startswith('#'):continue
                if header is None:raise ValueError('HPO association header not recognized')
                row=dict(zip(header,parts));g=aliases.get(row['gene_symbol'],row['gene_symbol'])
                direct[g].add(alts.get(row['hpo_id'],row['hpo_id']))
        self.direct=direct;self.closure={}
        counts=Counter()
        for g,terms in direct.items():
            inherited=set().union(*(self.ancestors(t) for t in terms))
            self.closure[g]=inherited;counts.update(inherited)
        total=len(direct)
        if not total:raise ValueError('Empty HPO associations')
        self.ic={t:-math.log((n+1)/(total+1)) for t,n in counts.items()}
        self.total_genes=total

    @lru_cache(maxsize=None)
    def ancestors(self,t):
        return frozenset({t}).union(*(self.ancestors(p) for p in self.parents[t]))

    def score(self,gene,query):
        query=set(self.alts.get(q,q) for q in query)
        # If a listed term is an ancestor of another listed term, count the more specific term only.
        redundant={q for q in query if any(q!=r and q in self.ancestors(r) for r in query)}
        details=[];denom=0;num=0
        for q in sorted(query-redundant):
            weight=self.ic.get(q,0)
            common=self.ancestors(q)&self.closure.get(gene,set())
            match=max(common,key=lambda t:(self.ic.get(t,0),t)) if common else None
            support=self.ic.get(match,0) if match else 0
            denom+=weight;num+=support
            details.append({'hpo_id':q,'label':self.labels.get(q),'information_content':weight,
                            'best_common_ancestor':match,'matched_information':support})
        return {'score':num/denom if denom else 0,'terms':details,'removed_redundant_terms':sorted(redundant),
                'annotation_present':gene in self.direct,'unscored_terms':[r['hpo_id'] for r in details if r['information_content']==0]}


def serial(value):
    if isinstance(value,tuple):return list(value)
    return value


def clinvar_index(path, transcripts):
    wanted=['CLNSIG','CLNREVSTAT','CLNDN','CLNDISDB','MC','AF_EXAC','AF_ESP','AF_TGP','CLNSIGCONF']
    result={}
    with pysam.VariantFile(str(path)) as vcf:
        for gene,t in transcripts.items():
            chrom=t['chrom'].removeprefix('chr')
            for rec in vcf.fetch(chrom,t['start'],t['end']):
                for alt in rec.alts or []:
                    pos,ref,alt=minimal_allele(rec.pos,rec.ref,alt)
                    key=(chrom,pos,ref,alt)
                    result[key]={'id':rec.id,**{k:serial(rec.info[k]) for k in wanted if k in rec.info}}
    return result


def text_field(value):
    return '|'.join(map(str,value)) if isinstance(value,(list,tuple)) else str(value or '')


def clinical_support(cv):
    sig=text_field(cv.get('CLNSIG')).lower();review=text_field(cv.get('CLNREVSTAT')).lower()
    conflict='conflict' in sig or bool(cv.get('CLNSIGCONF'))
    benign=('benign' in sig and 'pathogenic' not in sig and not conflict)
    # Exact tokens matter: no_assertion_criteria_provided is a zero-star status.
    review_tokens={part.strip('_ ') for part in re.split(r'[|,]',review)}
    stars=(4 if 'practice_guideline' in review_tokens
           else 3 if review_tokens & {'reviewed_by_expert_panel','expert_panel'}
           else 2 if {'criteria_provided','multiple_submitters','no_conflicts'}<=review_tokens
           else 1 if 'criteria_provided' in review_tokens else 0)
    pathogenic='pathogenic' in sig and not conflict
    value=([0.65,0.75,0.85,0.95,0.95][stars] if pathogenic else 0)
    return {'support':value,'benign':benign,'conflicting':conflict,'review_stars':stars,
            'pathogenic_assertion':pathogenic,'classification':text_field(cv.get('CLNSIG')) or 'not_in_ClinVar'}


def rarity_evidence(cv):
    values=[]
    for key in ['AF_EXAC','AF_ESP','AF_TGP']:
        data=cv.get(key,[])
        if not isinstance(data,(list,tuple)):data=[data]
        for val in data:
            if val is not None:
                try:
                    val=float(val)
                    if math.isfinite(val) and 0<=val<=1:values.append((key,val))
                except (TypeError,ValueError):pass
    if not values:
        return {'score':0.5,'known':False,'af':None,'sources':[],
                'caveat':'Population data missing; neutral placeholder, NOT evidence of rarity'}
    af=max(v for _,v in values)
    score=max(0,min(1,(-math.log10(max(af,1e-6))-2)/4))
    return {'score':score,'known':True,'af':af,'sources':values,
            'caveat':'Legacy population fields from ClinVar; not comprehensive gnomAD or ancestry filtering'}


def variant_features(v,config):
    clinical=clinical_support(v['clinvar']);rarity=rarity_evidence(v['clinvar'])
    effect=config['effect_scores'].get(v['consequence'],0)
    # Consequence and clinical evidence are overlapping signals: use the maximum, not their sum.
    functional=max(effect,clinical['support'])
    w=config['rarity_weight']
    value=(1-w)*functional+w*rarity['score']
    if clinical['benign']:value*=0.1
    flags=[]
    if not rarity['known']:flags.append('population_frequency_unknown')
    if clinical['conflicting']:flags.append('ClinVar_conflict_requires_review')
    if clinical['benign']:flags.append('ClinVar_benign_evidence')
    return {'score':value,'effect_support':effect,'clinical':clinical,'rarity':rarity,'flags':flags}


def phase_call(a,b):
    if a['record_id']==b['record_id'] and a['alt_index']!=b['alt_index']:
        if sorted(a['gt'])==sorted([a['alt_index'],b['alt_index']]):
            return 'trans','Distinct alternate alleles in one diploid 1/2-like genotype'
    if a['phase_set'] and a['phase_set']==b['phase_set'] and a['haplotype'] is not None and b['haplotype'] is not None:
        return ('cis' if a['haplotype']==b['haplotype'] else 'trans'),'Shared VCF phase block; read-based phase requires validation'
    return 'unknown','No shared informative VCF phase block; parental/orthogonal phasing required'


def generate_hypotheses(variants, phenotype, config):
    groups=defaultdict(list);review=[]
    for v in variants:
        features=variant_features(v,config);v['features']=features
        reason=[]
        if not v['quality_pass']:reason.extend(v['quality_flags'])
        if features['clinical']['benign']:reason.append('benign_ClinVar')
        if features['rarity']['known'] and features['rarity']['af']>config['maximum_population_af']:reason.append('population_af_above_threshold')
        if features['effect_support']<0.30 and not features['clinical']['pathogenic_assertion']:reason.append('no_functional_priority_in_this_model')
        if reason:
            review.append({'id':v['id'],'reasons':reason});continue
        groups[v['gene']].append(v)
    rows=[]
    for gene,vs in groups.items():
        combos=[([v],'not_applicable','Homozygous alternate genotype') for v in vs if v['genotype']=='hom_alt']
        for a,b in itertools.combinations([v for v in vs if v['genotype']=='het'],2):
            phase,evidence=phase_call(a,b)
            if phase=='cis':
                review.append({'members':[a['id'],b['id']],'reasons':['shared_block_cis']});continue
            combos.append(([a,b],phase,evidence))
        for members,phase,evidence in combos:
            vscore=min(v['features']['score'] for v in members)
            pscore=phenotype[gene]['score']
            w=config['phenotype_weight'];factor=config['unknown_phase_factor'] if phase=='unknown' else 1
            score=factor*((1-w)*vscore+w*pscore)
            rows.append({'gene':gene,'members':[v['id'] for v in members],
                'kind':'pair' if len(members)==2 else 'homozygous',
                'phase':phase,'phase_evidence':evidence,'score':score,
                'components':{'weakest_allele_score':vscore,'phenotype_score':pscore,'phase_factor':factor},
                'caveats':sorted(set(itertools.chain.from_iterable(v['features']['flags'] for v in members)))})
    rows.sort(key=lambda r:(-r['score'],r['gene'],r['members']))
    for i,row in enumerate(rows):row['rank']=i+1
    return rows,review


def inspect_variants(vcf_path,transcripts,cv,config):
    out=[];counts={}
    with pysam.VariantFile(str(vcf_path)) as vcf:
        header=str(vcf.header)
        if 'GRCh38' not in header:raise ValueError('VCF header does not establish GRCh38')
        if len(vcf.header.samples)!=1:raise ValueError('This baseline requires an explicitly selected single sample')
        sample_name=list(vcf.header.samples)[0]
        for gene,t in transcripts.items():
            dna=json.loads((t['_reference_dir']/t['reference_file']).read_text())['dna']
            transcript=Transcript(t,dna)
            counter=Counter();chrom=t['chrom'].removeprefix('chr')
            contig=t['chrom'] if t['chrom'] in vcf.header.contigs else chrom
            for rec in vcf.fetch(contig,t['start'],t['end']):
                counter['records_in_transcript_span']+=1
                sample=rec.samples[sample_name];gt=sample.get('GT')
                if not gt or len(gt)!=2 or None in gt:
                    counter['missing_or_non_diploid_gt']+=1;continue
                for ai in sorted(set(gt)-{0}):
                    alt=rec.alts[ai-1];ref=rec.ref
                    if not set(ref+alt)<=set('ACGT'):
                        counter['symbolic_or_nonsequence_alleles']+=1;continue
                    annotation=transcript.annotate(rec.pos,ref,alt)
                    flags=[]
                    if list(rec.filter)!=['PASS']:flags.append('FILTER_not_PASS')
                    dp=sample.get('DP');gq=sample.get('GQ');ad=sample.get('AD')
                    if dp is None or dp<config['minimum_dp']:flags.append('DP_missing_or_low')
                    if gq is None or gq<config['minimum_gq']:flags.append('GQ_missing_or_low')
                    if annotation['consequence']=='reference_mismatch':flags.append('reference_mismatch')
                    balance=ad[ai]/sum(ad) if ad and len(ad)>ai and all(i is not None for i in ad) and sum(ad)>0 else None
                    genotype='hom_alt' if gt==(ai,ai) else 'het'
                    phase_set=None;haplotype=None
                    if genotype=='het':
                        ps=sample.get('PS')
                        if sample.phased and ps is not None:
                            phase_set=f'{chrom}:PS:{ps}';haplotype=gt.index(ai)
                        elif sample.get('PID') and sample.get('PGT'):
                            pgt=sample['PGT']
                            if re.fullmatch(r'\d+\|\d+',pgt):
                                alleles=tuple(map(int,pgt.split('|')))
                                if sorted(alleles)==sorted(gt):
                                    phase_set=f"{chrom}:PID:{sample['PID']}";haplotype=alleles.index(ai)
                    pos,refkey,altkey=minimal_allele(rec.pos,ref,alt)
                    key=(chrom,pos,refkey,altkey)
                    out.append({'id':f'chr{chrom}:{rec.pos}:{ref}:{alt}', 'gene':gene,'chrom':f'chr{chrom}',
                        'pos':rec.pos,'ref':ref,'alt':alt,'alt_index':ai,
                        'record_id':f'{chrom}:{rec.pos}:{rec.ref}:{",".join(rec.alts)}',
                        'gt':list(gt),'genotype':genotype,'dp':dp,'gq':gq,'ad':serial(ad),
                        'allele_balance':balance,'phase_set':phase_set,'haplotype':haplotype,
                        'quality_pass':not flags,'quality_flags':flags,'filter':list(rec.filter),
                        'clinvar':cv.get(key,{}),**annotation})
                    counter['nonref_sequence_alleles']+=1
                    counter[annotation['consequence']]+=1
            counts[gene]=dict(counter)
    return out,counts,sample_name


def audit_vcf(path, transcripts=None):
    counts=Counter();depth=[];quality=[];panel_counts=Counter()
    by_contig=defaultdict(list)
    for gene,t in (transcripts or {}).items():
        by_contig[t['chrom'].removeprefix('chr')].append((gene,t['start'],t['end']))
    with pysam.VariantFile(str(path)) as vcf:
        for i,r in enumerate(vcf):
            counts['records']+=1
            if list(r.filter)==['PASS']:counts['PASS']+=1
            for gene,start,end in by_contig.get(r.contig.removeprefix('chr'),[]):
                if r.start<end and r.stop>start:panel_counts[gene]+=1
            if i%100==0:
                sample=next(iter(r.samples.values()))
                if sample.get('DP') is not None:depth.append(sample['DP'])
                if sample.get('GQ') is not None:quality.append(sample['GQ'])
    return {**counts,'sampled_every_n_records':100,'sampled_median_DP':statistics.median(depth) if depth else None,
            'sampled_median_GQ':statistics.median(quality) if quality else None,
            'panel_records_sequential':dict(panel_counts),
            'note':'Full VCF parse completed; QC statistics are not diagnostic accuracy'}


def sensitivity(variants,phenotype,config):
    scenarios={}
    for label,change in [('no_phenotype',{'phenotype_weight':0}),('phenotype_weight_0.15',{'phenotype_weight':.15}),
                         ('phenotype_weight_0.45',{'phenotype_weight':.45}),('no_phase_penalty',{'unknown_phase_factor':1}),
                         ('no_population_weight',{'rarity_weight':0})]:
        rows,_=generate_hypotheses(variants,phenotype,{**config,**change})
        scenarios[label]=[{'gene':r['gene'],'members':r['members'],'score':r['score']} for r in rows[:10]]
    # Removing ClinVar removes both classification and its legacy frequency fields: label explicitly.
    stripped=[{**v,'clinvar':{}} for v in variants]
    rows,_=generate_hypotheses(stripped,phenotype,config)
    scenarios['without_ClinVar_and_its_population_fields']=[{'gene':r['gene'],'members':r['members'],'score':r['score']} for r in rows[:10]]
    return scenarios


def report(out,ranked,variants,phenotypes,gene_pheno,audit,provenance,counts):
    by_id={v['id']:v for v in variants}
    text=['# Track 1 actual-data baseline report','',
          'Model: mva_panel_recessive_v1. Research candidate ranking, not a clinical diagnosis. Scores are uncalibrated utilities, not EPCR or disease probabilities.',
          '',f"Run: {provenance['run_utc']}",
          f"Scanned {audit['records']:,} VCF records ({audit['PASS']:,} PASS). Sampled median DP {audit['sampled_median_DP']}; median GQ {audit['sampled_median_GQ']}.",
          f'Annotated {len(variants):,} nonreference sequence alleles in six representative transcript spans; ranked {len(ranked)} eligible recessive hypotheses.',
          '', '## Scope and limitations','',
          'This is a literature-selected six-gene baseline. It cannot exclude other genes, structural variants, repeats, deep-intronic effects, alternative-transcript effects, or low-quality/mosaic calls. A correct result here would not establish genome-wide or patient-level generalization.',
          'The transcript is the longest complete curated RefSeq CDS per gene, not necessarily MANE. Reference CDS translation and candidate REF bases are checked locally. Indel annotations are limited; repeat left-alignment and formal HGVS are not performed.',
          'Population evidence uses only available legacy ClinVar AF fields. Missing frequency receives a neutral placeholder and is flagged; absent ClinVar is not evidence of rarity or benignity. A gnomAD-based follow-up remains necessary.',
          'Public participant methods revealing candidate alleles were consulted before this run. No participant candidate coordinate, answer key or leaderboard score is used by this model; it is not a blinded discovery.',
          '', '## Phenotype handling','',
          'HPO identifiers are read from the supplied clinical table. Family-history terms remain documented separately and are not treated as phenotypes of the child. Unmentioned findings are not coded absent.',
          '']
    for p in phenotypes:text.append(f"- {p['hpo_id']} {p['label']} — {p['subject']}")
    text.extend(['','## Gene phenotype scores','', '| Gene | HPO score | Annotation available |','|---|---:|---|'])
    for g,p in sorted(gene_pheno.items(),key=lambda x:-x[1]['score']):text.append(f"| {g} | {p['score']:.4f} | {p['annotation_present']} |")
    text.extend(['','## Leading hypotheses','', '| Rank | Gene | Alleles | Score | Phase |','|---:|---|---|---:|---|'])
    for row in ranked[:10]:text.append(f"| {row['rank']} | {row['gene']} | {' + '.join(row['members'])} | {row['score']:.4f} | {row['phase']} |")
    if ranked:
        top=ranked[0];text.extend(['','## Evidence for the leading hypothesis','',top['phase_evidence'],''])
        for key in top['members']:
            v=by_id[key];f=v['features'];cv=f['clinical']
            text.extend([f"### {v['id']}",'',
              f"Gene/transcript: {v['gene']} / {v['transcript']}. Predicted consequence: {v['consequence']}; {v.get('coding_change') or 'coding change not assigned'}; {v.get('protein_change') or 'protein change not assigned'}.",
              f"Genotype: {v['gt']}; DP {v['dp']}; GQ {v['gq']}; AD {v['ad']}; alternate-read fraction {v['allele_balance']}.",
              f"ClinVar: {cv['classification']}; review stars {cv['review_stars']}; accession {v['clinvar'].get('id','unmatched')}. Conditions: {text_field(v['clinvar'].get('CLNDN')) or 'unavailable'}.",
              f"Population AF: {f['rarity']['af']}; {f['rarity']['caveat']}.",
              f"Allele utility {f['score']:.4f}; consequence support {f['effect_support']:.2f}; clinical support {cv['support']:.2f}.",''])
        text.extend(['The pair (if any) remains a hypothesis. Functional effects, phase, transcript relevance, condition-specific ClinVar assertions and technical validity require review. No automated ACMG classification has been assigned.',''])
    text.extend(['## Reproduction and validation','',
                'See `run_manifest.json`, `gene_phenotype.json`, `candidate_evidence.json`, `sensitivity.json`, `review_queue.json` and `panel_counts.json` in this directory. The config, code, inputs and reference hashes are recorded. Software tests validate implementation behavior, not biological accuracy.',
                'The model was not trained on this child. No challenge score was requested and no submission was made.',
                '', '## Next evidence to obtain','',
                '1. Comprehensive population-frequency and variant-effect annotation, plus MANE/transcript review.',
                '2. Parental segregation or validated shared-block phasing of the leading pair.',
                '3. Independent genome-wide prioritization (for example Exomiser) and alternative-mechanism review.',
                '4. Read-level/orthogonal validation where needed and functional evidence for uncertain alleles.',
                '', 'See `docs/TRACK1_BASELINE.md` for biology, formulas, selection rationale and cited sources.'])
    (out/'model_report.md').write_text('\n'.join(text)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--vcf',type=Path,default=input_path('vcf'))
    p.add_argument('--phenotype',type=Path,default=input_path('phenotype'))
    p.add_argument('--references',type=Path,default=ROOT/'references/track1')
    p.add_argument('--config',type=Path,default=ROOT/'config/track1_panel.json')
    p.add_argument('--output',type=Path,default=ROOT/'outputs/track1_baseline')
    a=p.parse_args();start=time.perf_counter();a.output.mkdir(parents=True,exist_ok=True)
    config=json.loads(a.config.read_text())
    manifest=json.loads((a.references/'manifest.json').read_text())
    for name,expected in manifest.items():
        if digest(a.references/name)!=expected:raise ValueError(f'Reference hash mismatch: {name}')
    transcripts=json.loads((a.references/'transcripts.json').read_text())
    if set(transcripts)!=set(config['genes']):raise ValueError('Reference and configured gene panels differ')
    for t in transcripts.values():t['_reference_dir']=a.references
    phenotypes=extract_phenotypes(a.phenotype)
    matcher=PhenotypeMatcher(a.references/'hp.obo',a.references/'genes_to_phenotype.txt',config['aliases'])
    query=[x['hpo_id'] for x in phenotypes if x['subject']=='proband']
    gene_pheno={g:matcher.score(g,query) for g in config['genes']}
    cv=clinvar_index(a.references/'clinvar.vcf.gz',transcripts)
    variants,counts,sample=inspect_variants(a.vcf,transcripts,cv,config)
    print('Local panel annotation complete; auditing full VCF',flush=True)
    audit=audit_vcf(a.vcf,transcripts)
    for gene,c in counts.items():
        if c['records_in_transcript_span']!=audit['panel_records_sequential'].get(gene,0):
            raise ValueError(f'Indexed/sequential panel record counts disagree: {gene}; rebuild the index')
    audit['indexed_panel_counts_verified_against_sequential_scan']=True
    scenarios=sensitivity(variants,gene_pheno,config)
    ranked,review=generate_hypotheses(variants,gene_pheno,config)
    provenance={'run_utc':datetime.now(timezone.utc).isoformat(),'model':config['model'],
       'config':config,'config_sha256':digest(a.config),'code_sha256':digest(__file__),
       'vcf_sha256':digest(a.vcf),'index_sha256':digest(str(a.vcf)+'.tbi'),
       'phenotype_sha256':digest(a.phenotype),'reference_hashes':manifest,'pysam_version':pysam.__version__,
       'elapsed_seconds':time.perf_counter()-start,'sample':sample,
       'training':'none','network_during_analysis':False,'score_interpretation':'uncalibrated heuristic utility'}
    for name,obj in [('ranked_hypotheses',ranked),('candidate_evidence',variants),('gene_phenotype',gene_pheno),
                     ('phenotypes',phenotypes),('vcf_audit',audit),('panel_counts',counts),('review_queue',review),
                     ('sensitivity',scenarios),('run_manifest',provenance)]:write_json(a.output/f'{name}.json',obj)
    report(a.output,ranked,variants,phenotypes,gene_pheno,audit,provenance,counts)
    print(f'Finished: {len(ranked)} research hypotheses; report at {a.output / "model_report.md"}',flush=True)


if __name__=='__main__':
    main()
