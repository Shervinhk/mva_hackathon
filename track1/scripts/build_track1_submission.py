"""Build a private local Track 1 submission from offline evidence.

No network, account login, publication, or leaderboard action is performed.
Outputs contain findings and must not be returned to a hosted assistant.
"""
from __future__ import annotations

if __package__:
    from .local_inputs import input_path
else:
    from local_inputs import input_path
import argparse,csv,json,math,shutil,sys,zipfile
from datetime import datetime,timezone
from pathlib import Path
import pysam
from track1_baseline import digest,write_json,Transcript
from validate_track1_submission import FIELDS,validate_csv,validate_package
from track1_v2 import install_network_guard
ROOT=Path(__file__).resolve().parents[1]

ACK='This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.'

def cell(x):return str(x).replace('|',' / ').replace('\n',' ')

def check_findings(rows,by_id,vcf_path,fasta_path):
    checks={'ranked_candidates_present':bool(rows),'finite_ordered_scores':True,
            'proposed_alleles_in_original_calls':True,'reference_bases_match':True,'sequence_effect_crosscheck':[]}
    if not rows:raise ValueError('No ranked hypotheses; do not invent a submission')
    scores=[r['score'] for r in rows]
    if any(not math.isfinite(s) or not 0<=s<=1 for s in scores) or scores!=sorted(scores,reverse=True):raise ValueError('Ranking invalid')
    with pysam.VariantFile(str(vcf_path)) as original,pysam.FastaFile(str(fasta_path)) as fasta:
        for row in rows[:10]:
            for id in row['members']:
                v=by_id[id]
                if not v['quality_pass']:raise ValueError('Proposed allele failed call quality')
                actual=fasta.fetch(v['chrom'],v['pos']-1,v['pos']-1+len(v['ref'])).upper()
                if actual!=v['ref']:raise ValueError('Proposed allele reference mismatch')
                matched=False
                source_contig=v.get('source_contig',v['chrom'])
                if source_contig not in original.header.contigs:
                    source_contig=next(c for c in original.header.contigs if c.removeprefix('chr')==v['chrom'].removeprefix('chr'))
                for rec in original.fetch(source_contig,v['pos']-1,v['pos']):
                    if rec.pos!=v['pos'] or rec.ref!=v['ref'] or v['alt'] not in (rec.alts or []):continue
                    sample=next(iter(rec.samples.values()));idx=rec.alts.index(v['alt'])+1
                    if idx in sample.get('GT',()):matched=True
                if not matched:raise ValueError('Proposed allele not observed in original sample calls')
    # Independent original local sequence annotator: only compare the same available transcript coding effect.
    meta_path=ROOT/'references/track1/transcripts.json'
    if meta_path.exists():
        metas=json.loads(meta_path.read_text())
        for row in rows[:10]:
            if row['gene'] not in metas:continue
            m=metas[row['gene']];dna=json.loads((meta_path.parent/m['reference_file']).read_text())['dna']
            tx=Transcript(m,dna)
            for id in row['members']:
                v=by_id[id]
                if len(v['ref'])!=1 or len(v['alt'])!=1:continue
                a=tx.annotate(v['pos'],v['ref'],v['alt'])
                checks['sequence_effect_crosscheck'].append({'id':id,'refseq':m['transcript'],
                    'independent_annotation':a,'bcftools_annotations':v['annotations'],
                    'interpretation':'Two annotation engines; transcript versions must be compared before calling agreement'})
    return checks


def report(out,rows,by_id,manifest,scenarios,checks,screen):
    top=rows[0];variants=[by_id[k] for k in top['members']]
    lines=['# Track 1 methods and results','',
        '**Local submission draft.** Confirm the participant fields and confidence before uploading. No official evaluation has been requested.','',
        '## Identity and disclosure','',
        '- Team / participant: [TO CONFIRM AT FINAL REVIEW].',
        '- GitHub repository URL: [TO CONFIRM AT FINAL REVIEW].',
        f'- Approach: {manifest["model"]}.',
        '- AI assistance: OpenAI Codex for source-code development, public-source research and drafting; exact model history and product tier: [TO CONFIRM].',
        '- Training control: participant reports that “Improve the model for everyone” had been disabled long before this work. This is participant attestation; the assistant cannot inspect the account setting.',
        '- Feedback and any credit-program terms: [TO CONFIRM].',
        '- Patient-specific material appeared in the earlier conversation before the local-only restriction. Subsequent private analyses and outputs remained in local programs/files, with only generic execution status returned to the assistant. No retrospective guarantee is made.',
        '- No patient data were used to train or fine-tune a model by this pipeline. The use of pretrained AlphaMissense does not train it on this child.',
        '- Manual review: [TO CONFIRM BY PARTICIPANT]. The v2 results below were assembled by deterministic local code; they were not displayed to the hosted assistant for case-level review.','',
        '## Abstract','',
        'We developed a reproducible local variant-prioritization workflow for the MVA Track 1 challenge. '
        'BCFtools csq annotates primary human chromosome calls using the MANE 1.5 transcript set and GRCh38. '
        'Observed alleles are linked locally to ClinVar and to public gnomAD and AlphaMissense catalogues for a prespecified chromosome-instability gene set. '
        'Clinical features are compared with individual gene-disease HPO profiles, with inheritance used to construct single-allele and recessive hypotheses. '
        'Quality, population filters, condition-specific clinical support and available phase information determine candidate eligibility. '
        'A transparent score combines genetic evidence, phenotype similarity and literature mechanism support. '
        'A separate genome-wide ClinVar intersection checks for established annotated alleles outside the main candidate route. '
        'Ablations assess dependence on component choices. The resulting ranking is a hypothesis list, with unresolved phase and missing evidence explicitly reported. '
        'No clinical answer key, leaderboard tuning, new functional experiments or patient-trained model was used. '
        'Important limitations include incomplete population coverage outside the fixed catalogue, incomplete inheritance annotations, local per-record consequences and absence of structural/repeat/RNA analysis. '
        'Public participant approaches were reviewed, so this work is not a blinded discovery.','',
        '## Data and reference scope','',
        f'Run UTC: {manifest["run_utc"]}. Main analysis runtime: {manifest["elapsed_seconds"]:.1f} seconds, excluding reference downloads and initial genome unpacking. '
        'Local CPU execution; no cloud analysis job or paid inference API was invoked by the pipeline. Subscription and total development cost are not measured.',
        'The supplied single-sample VCF, its index and clinical DOCX are the inputs. FASTQ realignment was not performed. '
        'No train/test split was claimed from one patient or its technical sequencing lanes. '
        'The public MANE annotation is pinned to release 1.5. ClinVar/HPO snapshots and input/reference/source hashes are recorded in the local run manifest. '
        'gnomad_r4 is the API dataset label; no more specific minor version is invented. Public catalogue downloads covered eleven predefined genes.','',
        '## Approach and scoring','',
        (ROOT/'docs/TRACK1_V2.md').read_text().split('## Model definition\n',1)[1].split('## Limitations and prior work',1)[0],
        '## Leading hypothesis','',
        f'Gene: **{top["gene"]}**. Disease profile: {top["disease_id"]}. Inheritance model: {top["inheritance"]}. '
        f'Utility score: {top["score"]:.6f}; phase: **{top["phase"]}**. {top["phase_evidence"]}.',
        'This utility score is not a probability of diagnosis. The CSV uses a separate provisional subjective EPCR of 0.60 for the leading hypothesis. '
        'It is a conservative draft value, not a calibrated estimate or a value optimized against the hidden answer. The participant must review it before submission.','',
        '| Allele in GRCh38 | Called genotype | DP | GQ | Predicted effects | Population evidence |',
        '|---|---|---:|---:|---|---|']
    evidence_lines=[]
    for v in variants:
        ann=[a for a in v['annotations'] if a['gene']==top['gene']]
        effects='; '.join(f'{a["transcript"]} {"&".join(a["consequences"])} {a["amino_acid_change"]}' for a in ann)
        pop=v['population'];freq=str(pop.get('af')) if pop['known'] else 'unknown / unusable'
        lines.append('| '+' | '.join(map(cell,[v['id'],v['gt'],v['dp'],v['gq'],effects,f'{freq}; {pop["source"]}']))+' |')
        evidence_lines.extend(['',f'**Evidence for {v["id"]}**',''])
        for cv in v['clinvar_records']:
            evidence_lines.append(f'ClinVar {cv.get("id")}: {cell(cv.get("CLNSIG"))}; review {cell(cv.get("CLNREVSTAT"))}; conditions {cell(cv.get("CLNDN"))}; condition identifiers {cell(cv.get("CLNDISDB"))}.')
        public=v.get('public_annotation',{})
        if public.get('alphamissense'):
            am=public['alphamissense'];evidence_lines.append(f'Public AlphaMissense entry: {am.get("protein_variant")}; score {am.get("am_pathogenicity")}; class {am.get("am_class")}; transcript {am.get("transcript_id")}. Only a matching transcript and substitution contribute to model support.')
        if public.get('gnomad'):
            for platform in ['exome','genome','joint']:
                g=public['gnomad'].get(platform)
                if g:evidence_lines.append(f'gnomAD {platform}: AC {g.get("ac")}, AN {g.get("an")}, filters {g.get("filters")}; AF {g.get("af", "see AC/AN")}.')
    lines.extend(evidence_lines)
    lines.extend(['','Unresolved evidence flags: '+', '.join(top['flags'])+'.','',
        'Both technical validity and biological causality require review. For a recessive pair, trans phase and pathogenicity of both alleles remain separate questions. '
        'No chromosome mosaicism diagnosis is inferred from approximately balanced read counts. No ACMG classification is assigned.','',
        '## Alternative hypotheses and sensitivity','',
        'Alternatives are documented for review and are not automatically added as incidental findings to the compact CSV. '
        'A high utility with missing population data or nonspecific phenotype evidence is not sufficient for a secondary diagnosis.','',
        '| Rank | Gene | Disease profile | Hypothesis | Score | Phase |',
        '|---:|---|---|---|---:|---|'])
    for r in rows[:10]:lines.append('| '+' | '.join(map(cell,[r['rank'],r['gene'],r['disease_id'],' + '.join(r['members']),f'{r["score"]:.5f}',r['phase']]))+' |')
    lines.extend(['','| Ablation | Leading gene | Original leading hypothesis rank within top 20 |','|---|---|---:|'])
    identity=set(top['members'])
    for name,rs in scenarios.items():
        rank=next((i+1 for i,r in enumerate(rs) if set(r['members'])==identity),'not within top 20')
        lines.append(f'| {name} | {rs[0]["gene"] if rs else "none"} | {rank} |')
    baseline_path=ROOT/'outputs/track1_baseline/ranked_hypotheses.json'
    if baseline_path.exists():
        baseline=json.loads(baseline_path.read_text())
        if baseline:
            old=baseline[0]
            new_rank=next((r['rank'] for r in rows if set(r['members'])==set(old['members'])),None)
            lines.extend(['','### Comparison with the first baseline','',
                f'The first six-gene baseline prioritized {old["gene"]}: {" + ".join(old["members"])}. '
                f'That hypothesis has rank {new_rank if new_rank is not None else "not eligible"} in the broader v2 run. '
                'The two numerical utility scales differ and should not be compared as diagnostic probabilities. '
                'A changed rank can reflect broader competition, disease-specific HPO profiles or additional public functional/population evidence; it is not proof that either model is clinically correct.'])
    lines.extend(['','These are sensitivity analyses, not independent validation cohorts. We did not choose weights from the clinical key or challenge scores.','',
        '## Independent ClinVar route',''])
    if screen:
        hits=screen.get('prioritized_hits',[])
        lines.append(f'The independent screen retained {len(hits)} observed calls with pathogenic evidence or conflicts for review. '
                     'It checks exact alleles and genotypes, not whether a ClinVar label explains the full phenotype. '
                     'Its full local findings remain in outputs/track1_v2/clinvar_screen.md and are not distributed as a genotype table.')
    else:lines.append('Independent ClinVar screen output was unavailable at packaging; do not claim it ran.')
    lines.extend(['','## Validation and limitations','',
        'All proposed top-ten alleles were programmatically checked against the original called genotypes and GRCh38 reference bases. '
        'Synthetic regressions cover allele-specific annotation, cis/trans handling, exact ClinVar condition identifiers, review stars, population filtering, predictor transcript matching and the official CSV schema. '
        'These verify implementation behavior, not clinical accuracy. No challenge submission or score has been produced.',
        *['- '+s for s in manifest['limitations']],
        '- Manual case review and account/disclosure confirmation remain required.',
        '- The local reference-based SNV crosscheck is retained in validation_evidence.json; transcript differences must be reviewed before treating it as full concordance.','',
        '## Reproducibility and prior work','',
        'The accompanying repository contains pinned Python requirements, executable scripts, configurations, synthetic tests and commands. '
        'No raw genome, index, aligned reads, clinical source document or genotype-scale evidence table belongs in the repository. '
        'Public reference receipts and local run manifests preserve provenance. A small final candidate list and methods report can be retained under organizer guidance once publication is authorized.',
        'Public literature and participant repositories were consulted, including K-Dense, swoopeagle and ivillafuerte approaches documented in the local research memo. '
        'An incidental public performance claim was encountered but not used for tuning. Exomiser/Genomiser was evaluated as an approach but was not executed here. '
        'No comparator score is claimed as our own. The literature review found no direct functional validation for every possible missense candidate; predictions cannot fill that evidence gap.','',
        '## Sources','',
        (ROOT/'docs/TRACK1_V2.md').read_text().split('## Sources\n',1)[1],
        '- [MVA literature and functional checkpoint study](https://doi.org/10.1158/0008-5472.CAN-09-4319).',
        '- [MAD2L1BP primary study](https://insight.jci.org/articles/view/170079).',
        '- [CEP192 primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10716027/).',
        '- [SLF2 and SMC5 primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9636423/).',
        '- [Public K-Dense approach](https://github.com/K-Dense-AI/hf-mva-hackathon-2026).',
        '- [Public swoopeagle approach](https://github.com/swoopeagle/mva-hackathon-methods).',
        '- [Official submission and scoring source](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/tree/main).','',
        '## Required acknowledgement','',ACK,'',
        'Dataset citation for any publication must be completed from the organizer’s current Synapse citation: [TO CONFIRM BEFORE PUBLICATION].'])
    (out/'methods.md').write_text('\n'.join(lines)+'\n')


def prepare_repository(out):
    repo=out/'repository';repo.mkdir(exist_ok=True)
    for directory in ['scripts','config','tests']:
        dest=repo/directory;dest.mkdir(exist_ok=True)
        for src in (ROOT/directory).iterdir():
            if src.is_file() and src.suffix in {'.py','.json'}:shutil.copy2(src,dest/src.name)
    (repo/'docs').mkdir(exist_ok=True)
    shutil.copy2(ROOT/'docs/TRACK1_V2.md',repo/'docs/TRACK1_V2.md')
    for name in ['requirements-track1.txt','.gitignore']:shutil.copy2(ROOT/name,repo/name)
    (repo/'README.md').write_text('# MVA Track 1 local variant prioritizer\n\n'
        'See [the model and reproduction guide](docs/TRACK1_V2.md). Public reference download is separate from offline patient analysis. '
        'The submission report and findings are retained locally pending participant approval. No private inputs or results are included in this source package.\n\n'
        'Run synthetic checks with `.venv/bin/python -m unittest discover -s tests -v`.\n\n'
        'Python dependencies: `requirements-track1.txt`. Scientific weights are manually specified; the model is not clinically calibrated. '
        'No clinical answer key or leaderboard tuning was used.\n')
    (repo/'LICENSE').write_text('MIT License\n\nCopyright (c) 2026 MVA Track 1 contributors\n\n'
        'Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), '
        'to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, '
        'and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:\n\n'
        'The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.\n\n'
        'THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, '
        'FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER '
        'LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.\n\n'
        'This license applies to original source code only. External references, data and any challenge submission retain their own terms.\n')
    archive=out/'repository.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for path in sorted(repo.rglob('*')):
            if path.is_file():z.write(path,Path('mva-track1-v2')/path.relative_to(repo))
    return repo


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis',type=Path,default=ROOT/'outputs/track1_v2')
    parser.add_argument('--output',type=Path,default=ROOT/'outputs/submission_track1_v2')
    parser.add_argument('--vcf',type=Path,default=input_path('vcf'))
    parser.add_argument('--fasta',type=Path,default=ROOT/'references/track1_v2/hg38.fa')
    args=parser.parse_args();install_network_guard();out=args.output;out.mkdir(parents=True,exist_ok=True)
    rows=json.loads((args.analysis/'ranked_hypotheses.json').read_text())
    by_id={v['id']:v for v in json.loads((args.analysis/'candidate_evidence.json').read_text())}
    manifest=json.loads((args.analysis/'run_manifest.json').read_text())
    scenarios=json.loads((args.analysis/'sensitivity.json').read_text())
    checks=check_findings(rows,by_id,args.vcf,args.fasta)
    top=rows[0];line={k:'' for k in FIELDS};line.update(proband_id='PROBAND01',epcr='0.60',finding_type='primary',
        notes=f'{top["gene"]}; {top["inheritance"]}; phase {top["phase"]}; provisional subjective EPCR; see methods; functional uncertainty remains')
    for i,id in enumerate(top['members'],1):
        v=by_id[id]
        for name in ['chrom','pos','ref','alt']:line[f'{name}_{i}']=v[name]
    with (out/'predictions.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=FIELDS);writer.writeheader();writer.writerow(line)
    validation=validate_csv(out/'predictions.csv')
    if not validation['valid']:raise ValueError('Submission CSV validation failed; inspect locally')
    write_json(out/'validation_evidence.json',checks)
    screen_path=args.analysis/'clinvar_screen.json';screen=json.loads(screen_path.read_text()) if screen_path.exists() else None
    report(out,rows,by_id,manifest,scenarios,checks,screen)
    write_json(out/'package_validation.json',validate_package(out/'predictions.csv',out/'methods.md',require_package=True))
    shutil.copy2(args.analysis/'run_manifest.json',out/'run_manifest.json')
    prepare_repository(out)
    (out/'START_HERE.md').write_text('# Track 1 local submission package\n\n'
        'The analysis and local artifact construction are complete. No challenge entry has been uploaded and no submission slot was consumed.\n\n'
        '1. Review `methods.md` locally. It contains the leading prediction, supporting evidence, alternatives, ablations and limitations.\n'
        '2. Review `predictions.csv`. It contains one complete hypothesis; the EPCR value of 0.60 is provisional subjective confidence and requires participant review.\n'
        '3. Supply team name, actual GitHub account/repository, exact AI model/product tier, feedback/credit-term disclosure, and confirm the reported training opt-out. Complete placeholders in methods.\n'
        '4. Review the source-only `repository/` or `repository.zip`; choose repository name and license before publishing. Do not add the data, references or private analysis tables.\n'
        '5. Authorize release of the small candidate list and report if you want an external submission. The current local-only instruction remains in effect.\n'
        '6. Upload the CSV, methods Markdown/PDF and real GitHub URL using the official challenge Space. Confirm available quota first. There are at most six submissions.\n\n'
        'The local `package_validation.json` intentionally reports missing participant fields and GitHub URL. A valid CSV is not a completed disclosure or clinical validation. '
        'The source-only archive contains no genotype-scale data. `validation_evidence.json` is a private review artifact; do not upload it by default.\n\n'
        'Official interface: https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026\n\n'
        'Keep raw genomes, indexes, annotated BCF and genotype-scale JSON tables local and include them in the organizer-required deletion inventory. '
        'No deletion action has been performed.\n')
    write_json(out/'checksums.json',{p.name:digest(p) for p in out.iterdir() if p.is_file() and p.name!='checksums.json'})
    print('Local submission files created; participant metadata and authorization remain pending.')

if __name__=='__main__':main()
