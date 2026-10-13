"""Download public references only. This script never opens patient files."""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
from pathlib import Path
import time
import httpx

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def download(url, path):
    receipt = path.with_name(path.name + '.source.json')
    if path.exists() and receipt.exists():
        if sha256(path) != json.loads(receipt.read_text())['sha256']:
            raise ValueError(f'Cached reference changed: {path}')
        return
    temporary = path.with_name(path.name + '.partial')
    for attempt in range(3):
        try:
            with httpx.stream('GET', url, timeout=60, follow_redirects=True) as r:
                r.raise_for_status()
                with temporary.open('wb') as f:
                    for chunk in r.iter_bytes():
                        f.write(chunk)
                final_url = str(r.url)
            temporary.replace(path)
            receipt.write_text(json.dumps({'url':url,'resolved_url':final_url,
                'retrieved_utc':datetime.now(timezone.utc).isoformat(),
                'sha256':sha256(path),'bytes':path.stat().st_size},indent=2)+'\n')
            print('Downloaded', path.name, path.stat().st_size, flush=True)
            return
        except httpx.HTTPError:
            if attempt == 2:
                raise
            time.sleep(2)


def parse_transcripts(path, genes, aliases):
    by_gene = {g:[] for g in genes}
    with gzip.open(path,'rt') as f:
        for line in f:
            x = line.rstrip().split('\t')
            gene = aliases.get(x[12],x[12])
            if gene not in by_gene or not x[1].startswith('NM_') or '_' in x[2]:
                continue
            if x[13:15] != ['cmpl','cmpl']:
                continue
            starts = [int(i) for i in x[9].strip(',').split(',')]
            ends = [int(i) for i in x[10].strip(',').split(',')]
            cs,ce = int(x[6]),int(x[7])
            cds = [(max(s,cs),min(e,ce)) for s,e in zip(starts,ends) if max(s,cs)<min(e,ce)]
            by_gene[gene].append({'gene':gene,'transcript':x[1],'chrom':x[2], 'strand':x[3],
                'start':int(x[4]),'end':int(x[5]),'exons':list(zip(starts,ends)),
                'cds':cds,'coding_bases':sum(e-s for s,e in cds)})
    selected = {}
    for gene, rows in by_gene.items():
        if not rows:
            raise ValueError(f'No complete curated coding RefSeq transcript: {gene}')
        ordered = sorted(rows,key=lambda r:(-r['coding_bases'],r['transcript']))
        selected[gene] = ordered[0]
        selected[gene]['selection'] = 'Longest complete curated RefSeq CDS; transcript ID breaks ties; not necessarily MANE'
        selected[gene]['other_transcripts'] = [r['transcript'] for r in ordered[1:]]
    return selected


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'references/track1')
    p.add_argument('--config',type=Path,default=ROOT/'config/track1_panel.json')
    p.add_argument('--clinvar-date',default='20260923')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    config=json.loads(a.config.read_text())
    sources = {
      'ncbiRefSeqCurated.txt.gz':'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqCurated.txt.gz',
      'genes_to_phenotype.txt':'https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt',
      'hp.obo':'https://purl.obolibrary.org/obo/hp.obo',
      'clinvar.vcf.gz':f'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_{a.clinvar_date}.vcf.gz',
      'clinvar.vcf.gz.tbi':f'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_{a.clinvar_date}.vcf.gz.tbi'
    }
    for name,url in sources.items():
        download(url,a.output/name)
    transcripts=parse_transcripts(a.output/'ncbiRefSeqCurated.txt.gz',config['genes'],config['aliases'])
    for gene,t in transcripts.items():
        # Whole public transcript region plus flanks, irrespective of any patient calls.
        start=max(0,t['start']-20);end=t['end']+20
        url=f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={t['chrom']};start={start};end={end}"
        path=a.output/f'{gene}.sequence.json'
        download(url,path)
        seq=json.loads(path.read_text())
        if seq.get('chrom') != t['chrom'] or len(seq['dna']) != end-start:
            raise ValueError(f'Reference region mismatch: {gene}')
        t['reference_start']=start;t['reference_end']=end
        t['reference_file']=path.name
        print('Reference ready:',gene,t['transcript'],flush=True)
    (a.output/'transcripts.json').write_text(json.dumps(transcripts,indent=2)+'\n')
    manifest={p.name:sha256(p) for p in sorted(a.output.iterdir()) if p.is_file() and not p.name.endswith(('.partial','.source.json')) and p.name!='manifest.json'}
    (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Public reference bundle complete',flush=True)


if __name__=='__main__':
    main()
