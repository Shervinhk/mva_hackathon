"""Prepare cached public genome, MANE GFF and gene-level functional catalogues.

This script does not open any patient file. Network downloads are handled by
prepare_track1_references.py and prepare_functional_reference.py separately.
"""
from __future__ import annotations
import argparse, csv, gzip, json, shutil
from collections import defaultdict
from pathlib import Path
import pysam
from track1_baseline import digest, minimal_allele, write_json

ROOT=Path(__file__).resolve().parents[1]


def key(chrom,pos,ref,alt):
    pos,ref,alt=minimal_allele(int(pos),ref,alt)
    return ':'.join(map(str,[chrom.removeprefix('chr'),pos,ref,alt]))


def convert_mane(source,target):
    """Convert NCBI-distributed MANE attributes to bcftools Ensembl GFF dialect.

    Coordinates, strand, CDS phases and CDS inclusion of stop codons are kept.
    Noncoding MANE transcripts are outside this coding analysis.
    """
    counts=defaultdict(int)
    with gzip.open(source,'rt') as src,target.open('w') as dst:
        dst.write('##gff-version 3\n')
        for line in src:
            if line.startswith('#'):continue
            fields=line.rstrip('\n').split('\t')
            attrs=dict(x.split('=',1) for x in fields[8].split(';') if '=' in x)
            kind=fields[2]
            if attrs.get('gene_type')!='protein_coding':continue
            if kind=='gene':
                fields[8]=f'ID=gene:{attrs["ID"]};Name={attrs["gene_name"]};biotype=protein_coding'
            elif kind=='transcript':
                if attrs.get('transcript_type')!='protein_coding':continue
                fields[2]='mRNA'
                fields[8]=f'ID=transcript:{attrs["ID"]};Parent=gene:{attrs["Parent"]};biotype=protein_coding'
            elif kind in {'exon','CDS','five_prime_UTR','three_prime_UTR'}:
                if attrs.get('transcript_type')!='protein_coding':continue
                fields[8]=f'Parent=transcript:{attrs["Parent"]}'
            else:continue
            dst.write('\t'.join(fields)+'\n');counts[kind]+=1
    write_json(target.with_suffix('.provenance.json'),{'source_sha256':digest(source),
        'output_sha256':digest(target),'conversion':'Ensembl attribute dialect only; no coordinate or CDS phase edits',
        'public_feature_counts':dict(counts)})


def build_catalogue(directory):
    catalogue={}
    for path in sorted(directory.glob('*.gnomad_r4.json')):
        raw=json.loads(path.read_text())['data']['gene']
        for v in raw['variants']:
            k=key(raw['chrom'],v['pos'],v['ref'],v['alt'])
            item=catalogue.setdefault(k,{})
            item['gnomad']={t:v.get(t) for t in ['exome','genome','joint']}
            item.setdefault('genes',[]).append(raw['symbol'])
    for path in sorted(directory.glob('*.alphamissense_hg38.csv')):
        gene=path.name.split('.')[0]
        with path.open() as f:
            for row in csv.DictReader(f):
                if row['genome']!='hg38':raise ValueError('Public annotation build mismatch')
                k=key(row['CHROM'],row['POS'],row['REF'],row['ALT'])
                am={**row,'gene':gene,'am_pathogenicity':float(row['am_pathogenicity'])}
                item=catalogue.setdefault(k,{})
                item.setdefault('alphamissense_annotations',[]).append(am)
                item.setdefault('alphamissense',am)
    write_json(directory/'functional_catalogue.json',catalogue)
    write_json(directory/'functional_catalogue.provenance.json',{
        'source_hashes':{p.name:digest(p) for p in sorted(directory.iterdir())
            if p.name.endswith(('.gnomad_r4.json','.alphamissense_hg38.csv'))},
        'output_sha256':digest(directory/'functional_catalogue.json'),
        'public_reference_alleles':len(catalogue),'patient_inputs':False})


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--genome-dir',type=Path,default=ROOT/'references/track1_v2')
    p.add_argument('--functional-dir',type=Path,default=ROOT/'references/track1_v2_functional')
    p.add_argument('--download-genome',action='store_true',help='Download public UCSC hg38 only; never opens patient files')
    a=p.parse_args();fa=a.genome_dir/'hg38.fa'
    a.genome_dir.mkdir(parents=True,exist_ok=True)
    if a.download_genome:
        from prepare_track1_references import download
        download('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz',fa.with_suffix('.fa.gz'))
    # Validate cached public download receipts before using reference content.
    for directory in [a.genome_dir,a.functional_dir]:
        for receipt in directory.glob('*.source.json'):
            meta=json.loads(receipt.read_text());src=receipt.with_name(receipt.name.removesuffix('.source.json'))
            if not src.exists() or digest(src)!=meta['sha256']:raise ValueError('Public reference hash mismatch')
    if not fa.exists():
        temp=fa.with_suffix('.fa.partial')
        with gzip.open(fa.with_suffix('.fa.gz'),'rb') as src,temp.open('wb') as dst:shutil.copyfileobj(src,dst,4*1024*1024)
        temp.replace(fa)
    pysam.faidx(str(fa))
    gff=a.functional_dir/'MANE.GRCh38.v1.5.ensembl_genomic.gff.gz'
    convert_mane(gff,a.genome_dir/'MANE.v1.5.bcfcsq.gff3')
    build_catalogue(a.functional_dir)
    print('Public reference preparation completed')

if __name__=='__main__':main()
