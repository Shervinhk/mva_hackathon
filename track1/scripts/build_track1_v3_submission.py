"""Assemble a private v3 package with automated or explicitly reviewed export.

This program is offline. Never return its patient-derived files to a hosted model.
Automated EPCR values encode ordinal rank and are explicitly not probabilities.
"""
from __future__ import annotations

if __package__:
    from .local_inputs import input_path
else:
    from local_inputs import input_path

import argparse
import csv
import hashlib
import html
import json
import math
import shutil
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import pysam

from build_track1_submission import ACK
from track1_baseline import digest, write_json
from track1_v2 import install_network_guard
from track1_automated_export import automated_choices, POLICY as AUTOMATED_POLICY
from track1_entry_identity import validate_identity, participant_lines
from validate_track1_submission import FIELDS, validate_csv, validate_package

ROOT = Path(__file__).resolve().parents[1]
LIMITED_REFERENCE_NOTICE = ('Intermediate reference-limited run: the full Exomiser population database was not used. '
                            'This is a local workflow check; the full-reference analysis and comparison remain pending.')


def hypothesis_id(row):
    value = [row['gene'], row['disease_id'], row['inheritance'], sorted(row['members'])]
    return hashlib.sha256(json.dumps(value, separators=(',', ':')).encode()).hexdigest()[:24]


def cell(value):
    return str(value).replace('|', ' / ').replace('\n', ' ')


def safe_script_json(value):
    return json.dumps(value, allow_nan=False).replace('<', '\\u003c').replace('>', '\\u003e').replace('&', '\\u0026')


def review_manifest(rows, analysis_hash):
    return {
        'schema_version': 1,
        'analysis_sha256': analysis_hash,
        'reviewed_locally': False,
        'confidence_policy': 'EPCR must be supplied with reasons after local review; utility is not probability.',
        'identity': {'team': '', 'github_url': '', 'ai_disclosure': ''},
        'choices': [{'hypothesis_id': hypothesis_id(row), 'selected': False,
                     'epcr': None, 'finding_type': 'primary', 'reason': '',
                     'evidence_reviewed': False} for row in rows],
    }


def validate_choices(selection, rows, analysis_hash):
    """Fail closed on stale evidence, unreviewed choices and invented confidence."""
    if selection.get('schema_version') != 1 or selection.get('analysis_sha256') != analysis_hash:
        raise ValueError('Selection is from a different analysis or schema.')
    if selection.get('reviewed_locally') is not True:
        raise ValueError('Local review has not been confirmed.')
    lookup = {hypothesis_id(row): row for row in rows}
    selected = []
    seen = set()
    allele_sets = set()
    for choice in selection.get('choices', []):
        if choice.get('selected') is not True:
            continue
        key = choice.get('hypothesis_id')
        if key not in lookup or key in seen:
            raise ValueError('Unknown or duplicate selected hypothesis.')
        seen.add(key)
        if lookup[key].get('eligible_for_manual_selection') is False:
            raise ValueError('This record is an incomplete genotype or unresolved-inheritance review finding, not a complete submission hypothesis.')
        alleles = tuple(sorted(lookup[key]['members']))
        if alleles in allele_sets:
            raise ValueError('The same allele hypothesis cannot be submitted twice under different disease profiles.')
        allele_sets.add(alleles)
        confidence = choice.get('epcr')
        if isinstance(confidence, bool) or not isinstance(confidence, (int, float)):
            raise ValueError('Supply an explicit numeric EPCR.')
        if not math.isfinite(confidence) or not 0 < confidence <= 1:
            raise ValueError('EPCR must be in (0, 1].')
        if choice.get('evidence_reviewed') is not True or len(str(choice.get('reason', '')).strip()) < 15:
            raise ValueError('Each selected hypothesis needs evidence review and a confidence rationale.')
        if choice.get('finding_type') not in {'primary', 'secondary'}:
            raise ValueError('Finding type must be primary or secondary.')
        selected.append((lookup[key], choice))
    if not 1 <= len(selected) <= 10:
        raise ValueError('Select between one and ten justified hypotheses.')
    return sorted(selected, key=lambda item: -item[1]['epcr'])


def csv_row(row, choice, evidence):
    result = {field: '' for field in FIELDS}
    result.update(proband_id='PROBAND01', epcr=choice['epcr'], finding_type=choice['finding_type'],
                  notes=f"{row['gene']}; {row['inheritance']}; phase {row['phase']}; " + str(choice['reason']).strip())
    if not 1 <= len(row['members']) <= 2:
        raise ValueError('Official schema supports one or two alleles per hypothesis.')
    for i, member in enumerate(row['members'], 1):
        variant = evidence[member]
        for name in ('chrom', 'pos', 'ref', 'alt'):
            result[f'{name}_{i}'] = variant[name]
    return result


def verify_calls(rows, evidence, vcf_path, fasta_path):
    """Check original observed GT separately from normalized submission allele."""
    checked = set()
    with pysam.VariantFile(str(vcf_path)) as original, pysam.FastaFile(str(fasta_path)) as fasta:
        for row in rows:
            for member in row['members']:
                if member in checked:
                    continue
                variant = evidence[member]
                if not variant.get('quality_pass'):
                    raise ValueError('Review candidate failed call-quality checks.')
                chrom, pos, ref, alt = [variant[name] for name in ('chrom', 'pos', 'ref', 'alt')]
                if fasta.fetch(chrom, pos - 1, pos - 1 + len(ref)).upper() != ref:
                    raise ValueError('Normalized allele disagrees with reference sequence.')
                source_chrom = variant.get('source_contig', variant.get('source_chrom', chrom))
                source_pos = variant.get('source_pos', pos)
                source_ref = variant.get('source_ref', ref)
                source_alt = variant.get('source_alt', alt)
                if source_chrom not in original.header.contigs:
                    source_chrom = next(c for c in original.header.contigs
                                        if c.removeprefix('chr') == source_chrom.removeprefix('chr'))
                found = False
                for record in original.fetch(source_chrom, source_pos - 1, source_pos):
                    if record.pos != source_pos or record.ref != source_ref or source_alt not in (record.alts or []):
                        continue
                    sample = next(iter(record.samples.values()))
                    index = record.alts.index(source_alt) + 1
                    if index in (sample.get('GT') or ()):
                        found = True
                if not found:
                    raise ValueError('Source allele is not present in the original called genotype.')
                checked.add(member)
    return {'original_genotypes_checked': True, 'normalized_reference_checked': True,
            'review_alleles_checked': len(checked), 'clinical_causality_verified': False}


def create_review_html(path, rows, evidence, selection, notice=''):
    payload = {'selection': selection, 'fields': FIELDS, 'rows': [
        {'hypothesis_id': hypothesis_id(row), 'row': row,
         'variants': [evidence[member] for member in row['members']]} for row in rows]}
    template = '''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta http-equiv="Content-Security-Policy" content="default-src 'none'; script-src 'unsafe-inline'; style-src 'unsafe-inline'; connect-src 'none'; img-src data:; base-uri 'none'; form-action 'none'">
<title>Track 1 — local evidence review</title><style>
body{font:16px/1.55 system-ui,sans-serif;color:#183445;background:#f3f6f7;margin:0}main{max-width:1160px;margin:auto;padding:30px}h1{font-size:32px;margin-bottom:6px}h2{font-size:21px}p{max-width:900px}.notice{background:#fff2d8;padding:18px;border-left:4px solid #b76b00;border-radius:6px}.card{background:white;padding:22px;margin:20px 0;border:1px solid #cbd8df;border-radius:10px}.meta{color:#456270;font-size:14px}button{background:#096477;color:white;border:0;border-radius:6px;padding:12px 18px;font:inherit;cursor:pointer;margin:6px 8px 6px 0}input,select,textarea{font:inherit;padding:8px;border:1px solid #829ba8;border-radius:4px}input[type=number]{width:100px}input[type=checkbox]{width:18px;height:18px}textarea{display:block;width:95%;min-height:64px}label{display:inline-block;margin:8px 14px 8px 0}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f1f5f7;padding:12px;font-size:12px}summary{cursor:pointer;font-weight:600}.error{color:#a32121;font-weight:600}.toolbar{position:sticky;bottom:0;background:#e1edf0;padding:14px;border-top:1px solid #9cb6c1}.mono{font-family:ui-monospace,monospace;overflow-wrap:anywhere}.badge{display:inline-block;border-radius:4px;background:#e4eef2;padding:2px 7px;margin-right:8px}a{color:#086477}
</style></head><body><main><p class="meta">MVA Track 1 · v3 · Private local document</p><h1>Review the evidence, then choose the submission</h1>
<p>This page contains patient-derived findings. Open it locally. Its code makes no network requests; do not upload it or share it with an AI assistant.</p>
__RUN_NOTICE__
<div class="notice"><strong>No confidence has been invented.</strong> The ranking score orders hypotheses; it is not a probability. Select only explanations you can justify, supply the challenge EPCR confidence and explain your judgment. Unresolved phase and computational predictions remain limitations.</div>
<p>Read <a href="methods.md">methods.md</a> and <a href="comparison.md">comparison.md</a> alongside the evidence. The list below is a review queue, not a set of diagnoses. Distinguish candidate explanations from supported secondary findings.</p>
<div id="candidates"></div><div class="toolbar"><label><input id="reviewed" type="checkbox"> I reviewed the selected hypotheses and their uncertainty locally.</label><br><button id="manifest">Save reviewed selection</button><button id="csv">Save predictions CSV</button><div id="message" role="status"></div></div>
<p>Saving a file does not upload it or consume a submission slot. The Python finalizer can independently validate a saved selection and regenerate the CSV. Identity, GitHub and AI-use disclosure still need completion in the methods document.</p></main>
<script id="payload" type="application/json">__PAYLOAD__</script><script>
'use strict';
const data=JSON.parse(document.getElementById('payload').textContent);
const elements=[];
function node(tag,text,cls){const e=document.createElement(tag);if(text!==undefined)e.textContent=text;if(cls)e.className=cls;return e;}
for(const item of data.rows){
 const r=item.row,c=node('section',undefined,'card');
 const heading=node('h2',String(r.rank)+'. '+r.gene+' — '+r.disease_id);c.append(heading);
 c.append(node('p','Inheritance: '+r.inheritance+' · Phase: '+r.phase+' · Ranking utility: '+Number(r.score).toFixed(4),'meta'));
 c.append(node('p',r.members.join(' + '),'mono'));
 c.append(node('p','Review tier: '+(r.review_tier||'unresolved')+' · '+(r.flags||[]).join('; ')));
 const include=node('input');include.type='checkbox';const includeLabel=node('label');includeLabel.append(include,document.createTextNode(' Select hypothesis'));c.append(includeLabel);
 const confidence=node('input');confidence.type='number';confidence.min='0.000001';confidence.max='1';confidence.step='any';confidence.placeholder='EPCR';confidence.setAttribute('aria-label','EPCR confidence for hypothesis '+r.rank);c.append(confidence);
 const type=node('select');for(const value of ['primary','secondary']){const o=node('option',value);o.value=value;type.append(o);}type.setAttribute('aria-label','Finding type');c.append(type);
 const reason=node('textarea');reason.placeholder='Why this hypothesis, and why this confidence? Address both alleles, disease mechanism, population evidence and phase where relevant.';reason.setAttribute('aria-label','Evidence and confidence rationale');c.append(reason);
 const reviewed=node('input');reviewed.type='checkbox';const label=node('label');label.append(reviewed,document.createTextNode(' I reviewed this hypothesis’s supporting and conflicting evidence.'));c.append(label);
 const details=node('details');details.append(node('summary','Inspect the local evidence record'),node('pre',JSON.stringify({hypothesis:r,variants:item.variants},null,2)));c.append(details);
 document.getElementById('candidates').append(c);elements.push({item,include,confidence,type,reason,reviewed});
}
function collect(){
 if(!document.getElementById('reviewed').checked)throw Error('Confirm local review before exporting.');
 const chosen=elements.filter(e=>e.include.checked);if(chosen.length<1||chosen.length>10)throw Error('Choose one to ten justified hypotheses.');
 const selection=structuredClone(data.selection);selection.reviewed_locally=true;selection.choices=[];const seen=new Set();
 for(const e of chosen){const key=JSON.stringify([...e.item.row.members].sort());if(seen.has(key))throw Error('The same allele hypothesis cannot be submitted twice under different disease profiles.');seen.add(key);const p=Number(e.confidence.value);if(!e.confidence.value||!Number.isFinite(p)||p<=0||p>1)throw Error('Supply EPCR in (0, 1] for each selected hypothesis.');if(!e.reviewed.checked||e.reason.value.trim().length<15)throw Error('Review each selected hypothesis and give an evidence/confidence rationale.');selection.choices.push({hypothesis_id:e.item.hypothesis_id,selected:true,epcr:p,finding_type:e.type.value,reason:e.reason.value.trim(),evidence_reviewed:true});}
 selection.choices.sort((a,b)=>b.epcr-a.epcr);return selection;
}
function download(name,text,mime){const u=URL.createObjectURL(new Blob([text],{type:mime}));const a=node('a');a.href=u;a.download=name;document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(u),2000);}
function protect(fn){try{fn();document.getElementById('message').textContent='File prepared locally. No upload was made.';document.getElementById('message').className='';}catch(e){document.getElementById('message').textContent=e.message;document.getElementById('message').className='error';}}
document.getElementById('manifest').onclick=()=>protect(()=>download('selection.reviewed.json',JSON.stringify(collect(),null,2)+'\\n','application/json'));
document.getElementById('csv').onclick=()=>protect(()=>{const s=collect();const rows=[];for(const c of s.choices){const item=data.rows.find(x=>x.hypothesis_id===c.hypothesis_id);const r=item.row;const out=Object.fromEntries(data.fields.map(f=>[f,'']));Object.assign(out,{proband_id:'PROBAND01',epcr:c.epcr,finding_type:c.finding_type,notes:r.gene+'; '+r.inheritance+'; phase '+r.phase+'; '+c.reason});item.variants.forEach((v,i)=>{for(const k of ['chrom','pos','ref','alt'])out[k+'_'+(i+1)]=v[k];});rows.push(out);}const quote=v=>'"'+String(v).replaceAll('"','""')+'"';const csv=[data.fields,...rows.map(r=>data.fields.map(f=>r[f]))].map(r=>r.map(quote).join(',')).join('\\r\\n')+'\\r\\n';download('predictions.csv',csv,'text/csv');});
</script></body></html>'''
    run_notice = '<p class="notice"><strong>' + html.escape(notice) + '</strong></p>' if notice else ''
    path.write_text(template.replace('__PAYLOAD__', safe_script_json(payload)).replace('__RUN_NOTICE__', run_notice))


def allele_identity(members, fasta=None):
    from track1_v3 import normalize_allele
    keys = []
    for member in members:
        chrom, pos, ref, alt = member.split(':')
        keys.append(normalize_allele(chrom, int(pos), ref, alt, fasta))
    return tuple(sorted(set(keys)))


def write_comparison(path, rows, previous, comparator=None, fasta_path=None):
    fasta = pysam.FastaFile(str(fasta_path)) if fasta_path else None
    try:
        current_identities = {hypothesis_id(row): allele_identity(row['members'], fasta) for row in rows}
        previous_identities = {hypothesis_id(row): allele_identity(row['members'], fasta) for row in previous[:10]}
        compared = []
        for row in (comparator or {}).get('hypotheses', []):
            try:
                key = allele_identity(row['variant_ids'], fasta)
            except (ValueError, KeyError):
                continue
            compared.append((row, key))
    finally:
        if fasta:
            fasta.close()
    lines = ['# Comparison with v2', '',
             'Local patient-derived findings. Rank changes do not establish clinical correctness; the utility scales are not calibrated probabilities.', '',
             '| Version | Leading hypothesis | Phase |', '|---|---|---|']
    for label, hypotheses in [('v2', previous), ('v3', rows)]:
        if hypotheses:
            row = hypotheses[0]
            lines.append('| ' + ' | '.join(map(cell, [label, row['gene'] + ': ' + ' + '.join(row['members']), row['phase']])) + ' |')
    lines += ['', '| v2 rank | Hypothesis | v3 rank (normalized allele identity) |', '|---:|---|---:|']
    for old in previous[:10]:
        new = next((r['rank'] for r in rows if current_identities[hypothesis_id(r)] == previous_identities[hypothesis_id(old)]), None)
        lines.append('| ' + ' | '.join(map(cell, [old['rank'], old['gene'] + ': ' + ' + '.join(old['members']), new or 'not present in complete ranked hypotheses'])) + ' |')
    lines += ['', 'Equivalent sequence indels are compared after reference-aware normalization when the FASTA is supplied. '
              'Expanded public evidence, condition-specific assertions and changed eligibility can all change the ranking. See the evidence records for reasons.']
    lines += ['', '## Independent Exomiser comparison', '',
              'No numerical score averaging is performed. Agreement on a gene is weaker evidence than agreement on the complete allele hypothesis; neither establishes a diagnosis.', '']
    if comparator is None:
        lines.append('No completed comparator output was available for this package build. No agreement is claimed.')
    else:
        lines += ['| v3 rank | Gene | Exomiser exact-allele rank | Exomiser same-gene rank |', '|---:|---|---:|---:|']
        for row in distinct_review_rows(rows, 20):
            key = current_identities[hypothesis_id(row)]
            exact = [other['rank'] for other, other_key in compared if key == other_key]
            gene = [other['rank'] for other, _ in compared if row['gene'] == other['gene']]
            lines.append('| ' + ' | '.join(map(cell, [row['rank'], row['gene'], min(exact) if exact else 'no exact match', min(gene) if gene else 'not exported'])) + ' |')
        lines += ['', 'A comparator record with more than two contributing alleles is not a confirmed recessive pair. '
                  'Allele agreement does not establish phase, and inheritance explanations require separate review. '
                  'Both pipelines share some public reference evidence and are not statistically independent.']
    path.write_text('\n'.join(lines) + '\n')


def distinct_review_rows(rows, limit=30):
    """One choice per allele set; retain other disease interpretations as context."""
    groups = {}
    for row in rows:
        if row.get('eligible_for_manual_selection') is False:
            continue
        key = tuple(sorted(row['members']))
        if key not in groups:
            groups[key] = {**row, 'alternative_disease_interpretations': []}
        else:
            groups[key]['alternative_disease_interpretations'].append({
                field: row.get(field) for field in ['gene', 'disease_id', 'inheritance', 'score', 'review_tier']})
    return list(groups.values())[:limit]


def source_archive(out):
    repo = out / 'repository'
    repo.mkdir(exist_ok=True)
    # This previously generated copy is no longer in the public allowlist.
    # Keep the original project document; remove only its generated duplicate.
    legacy_copy = repo / 'docs/TRACK1_V2.md'
    if legacy_copy.is_file():
        legacy_copy.unlink()
    allow = {'scripts': {'.py'}, 'tests': {'.py', '.cjs'}, 'config': {'.json'}}
    files = []
    for folder, extensions in allow.items():
        (repo / folder).mkdir(exist_ok=True)
        for source in sorted((ROOT / folder).iterdir()):
            if source.is_file() and source.suffix in extensions:
                target = repo / folder / source.name
                shutil.copy2(source, target)
                files.append(target)
    (repo / 'docs').mkdir(exist_ok=True)
    for name in ['TRACK1_V3.md', 'AI_USE_AND_PRIVACY.md']:
        source = ROOT / 'docs' / name
        if source.exists():
            target = repo / 'docs' / name
            shutil.copy2(source, target)
            files.append(target)
    (repo / 'research').mkdir(exist_ok=True)
    for name in ['TRACK1_V3_PUBLIC_REFERENCES.md', 'TRACK1_V3_COMPARATOR.md']:
        source = ROOT / 'research' / name
        if source.exists():
            target = repo / 'research' / name
            shutil.copy2(source, target)
            files.append(target)
    for name in ['requirements-track1.txt', '.gitignore']:
        target = repo / name
        shutil.copy2(ROOT / name, target)
        files.append(target)
    # The legacy ranker tests depend on this explicitly fabricated fixture.
    # Do not broaden this to all example files: only known public inputs belong
    # in the shareable repository.
    fixture = Path('examples/synthetic_candidates.json')
    (repo / fixture.parent).mkdir(exist_ok=True)
    shutil.copy2(ROOT / fixture, repo / fixture)
    files.append(repo / fixture)
    readme = repo / 'README.md'
    readme.write_text('# MVA Track 1 v3\n\nSee [the complete approach](docs/TRACK1_V3.md). '
                      'This archive includes original source, configurations and fabricated tests only. '
                      'Public references must be downloaded separately under their own terms. '
                      'Private input data and results are deliberately excluded.\n\n'
                      'Set `MVA_VCF_PATH` and `MVA_PHENOTYPE_PATH` to your authorized local inputs, '
                      'or configure `work/private_analysis.nosync/local_inputs.json` with `vcf` and `phenotype` paths. '
                      'This private configuration is excluded from the archive. Without overrides, defaults are '
                      '`data/input.vcf.gz` and `data/phenotype.docx`.\n\n'
                      'Run synthetic tests with `python -m unittest discover -s tests -v`. '
                      'Checks requiring the separately downloaded Java adapter or retained official scorer '
                      'snapshot skip when those public resources are absent. '
                      'A research ranking is not a clinical diagnosis or a calibrated probability.\n')
    files.append(readme)
    with zipfile.ZipFile(out / 'repository.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
        for target in sorted(files):
            archive.write(target, Path('mva-track1-v3') / target.relative_to(repo))
    return {'files': len(files), 'sha256': digest(out / 'repository.zip'),
            'private_inputs_or_outputs_included': False}


def submission_report(out, selected, manifest, comparator, selection):
    """Brief organizer-facing report; detailed patient evidence stays separate."""
    automated = selection.get('selection_mode') == 'automated_ordinal'
    cfg = manifest.get('config', {})
    elapsed = manifest.get('elapsed_seconds')
    runtime = f'{elapsed / 60:.1f} minutes' if isinstance(elapsed, (int, float)) else 'not recorded'
    lines = ['# MVA Track 1 — disease-aware variant ranking v3', '',
             *participant_lines(selection.get('identity', {})), '',
             '## Abstract', '',
             'We prioritize small-variant hypotheses from the supplied GRCh38 VCF and clinical document. '
             'Local BCFtools csq annotation with MANE transcripts, reference-aware allele normalization, '
             'call-quality filtering, public population frequencies, condition-specific ClinVar assertions '
             'and disease-specific HPO similarity contribute to an interpretable, untrained ranking. '
             'Inheritance-compatible single alleles and recessive pairs are retained; unresolved phase '
             'is not presented as confirmed compound heterozygosity. An independent Exomiser run provides '
             'comparison without score averaging. The method preserves evidence provenance and missingness. '
             'It has not been calibrated or evaluated on an independent clinical cohort and does not cover '
             'all structural, repeat, mosaic or regulatory causes.', '',
             '## Inputs and public resources', '',
             'The authorized gated challenge inputs are the existing VCF and clinical document. FASTQ reads '
             'were not realigned and no new variant caller was run. Positive proband HPO observations are '
             'extracted locally; family observations and explicitly absent features are not treated as '
             'positive proband terms. No proprietary reference cohort or commercial annotation service was used.', '',
             'Auxiliary resources include GRCh38, MANE 1.5, the Human Phenotype Ontology and gene–disease '
             'annotations, complete ClinVar variant/submission summaries, NCBI disease-identifier mappings, '
             'and Exomiser 2602 GRCh38/phenotype releases. Earlier targeted public predictor annotations '
             'remain identifiable as a separate source. Exact versions, hashes and effective configuration '
             'are retained in the private run manifest; reference acquisition and source links are documented '
             'in the source repository.', '',
             '## Method and biology', '',
             'For most autosomal genes, a child has two copies. A recessive explanation can require variants '
             'affecting both copies; finding two variants in one gene does not establish that they are on '
             'opposite parental copies. The model preserves that distinction and does not infer de novo status.', '',
             f'Calls require PASS, depth at least {cfg.get("minimum_dp", 10)}, genotype quality at least '
             f'{cfg.get("minimum_gq", 20)}, and the configured heterozygous allele-balance interval. '
             'Original called alleles are retained alongside normalized GRCh38 representations. Missing '
             'population annotation is not interpreted as frequency zero. Population sources and their '
             'available quality information are retained. Condition-specific clinical support requires '
             'an actual matching SCV/RCV assertion; aggregate gene or ClinVar labels are not transferred '
             'to unrelated diseases.', '',
             'The default ranking combines genetic evidence (60%) and disease-specific phenotype similarity '
             '(40%). For a recessive pair the weaker allele contributes 85% of the genetic component and '
             'the stronger allele 15%. Predictor and literature-mechanism bonuses are zero by default; '
             'predictors and mechanism caveats remain visible, and sensitivity analyses explore alternative '
             'weights. Novel dominant and X-linked hypotheses are not excluded merely for lacking a ClinVar '
             'assertion. Incomplete recessive and unknown-inheritance findings remain outside the complete '
             'hypothesis ranking. The supplied configuration and source specify the full scoring rules.', '',
             '## Automated output or downstream review', '',
             (AUTOMATED_POLICY if automated else
              'A participant-supplied reviewed manifest determines the exported choices and confidence. '
              'See its recorded rationales and analysis-hash binding.'), '',
             'The automated shortlist preserves the existing v3 order and deduplicates exact normalized '
             'allele sets. It does not tune to leaderboard feedback, the hidden answer, or comparator '
             'agreement. These are alternative primary hypotheses, not an assertion that every row is causal. '
             'No secondary/incidental diagnosis is asserted by this automated export.', '',
             'The submission instructions say relative ranking controls scoring. The formal EPCR field '
             'is labelled as a probability, but the automated values here explicitly carry only ordinal '
             'ordering information. Absolute causal probabilities have not been estimated or calibrated; '
             'these values must not be used for clinical risk interpretation.', '',
             '## Ranked submission', '',
             '| Submission rank | Gene / disease | Allele hypothesis (GRCh38) | Inheritance | Phase |',
             '|---:|---|---|---|---|']
    for index, (row, _) in enumerate(selected, 1):
        values = [index, row['gene'] + ' / ' + row['disease_id'], ' + '.join(row['members']),
                  row['inheritance'], row['phase']]
        lines.append('| ' + ' | '.join(map(cell, values)) + ' |')
    lines += ['', '## Validation, runtime and limitations', '',
              'Each exported allele was checked against the original called genotype and the GRCh38 '
              'reference. CSV schema, unique hypotheses, complete pairs and strict rank ordering were '
              'checked locally. Synthetic tests assess implementation behavior, not clinical accuracy. '
              'The public Exomiser demonstration is an integration check, not a held-out accuracy estimate.', '',
              ('The production Exomiser comparison completed and matched the current input provenance. '
               'Its allele-level agreement is available in the local comparison report. Shared public '
               'evidence means the two methods are not statistically independent.'
               if comparator.get('eligible_for_current_comparison') else
               'No current verified comparator is claimed for this report.'), '',
              f'The main ranking run recorded {runtime}; this excludes initial reference downloads and '
              'preparation, prior reusable annotation, comparator time and packaging. Execution used local '
              'CPU/storage. No total financial cost is estimated. No model was trained or fine-tuned on this child.', '',
              'Limitations include incomplete population callability, transcript coverage, unknown phase, '
              'limited mechanism validation and lack of an independent calibration cohort. Public participant '
              'methods were encountered during project research; this is not described as blinded discovery. '
              'The v3 ranking was not changed during submission-file preparation.', '',
              '## AI use and data handling', '',
              'OpenAI Codex assisted with public-source research, software and documentation. The participant '
              'reports training opt-out was enabled before this work; the account setting and plan/tier '
              'were not independently verified. The participant disclosure above must specify actual product '
              'terms and controls. Earlier conversation contained patient-specific material; later source '
              'review also encountered case-specific filenames that have now been removed from shareable '
              'code. Current v3 genomic processing and findings stayed in deterministic local programs and '
              'excluded local files. No retrospective no-disclosure or no-training guarantee is made.', '',
              '## References and acknowledgement', '',
              '- SageBio. MVA Hackathon 2026 Data. Authorized challenge dataset: '
              'https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data . '
              'Organizer Synapse project: https://www.synapse.org/Synapse:syn76251147/wiki/642892 . '
              'Use the organizer-current formal dataset citation at external publication; no distinct formal '
              'citation was located in the public landing wiki during preparation.',
              '- Submission instructions: https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/submit_track1.py',
              '- BCFtools csq: https://samtools.github.io/bcftools/howtos/csq-calling.html',
              '- MANE: https://www.ncbi.nlm.nih.gov/refseq/MANE/',
              '- ClinVar assertions: https://www.ncbi.nlm.nih.gov/clinvar/docs/details/',
              '- Exomiser: https://github.com/exomiser/Exomiser', '', ACK, '']
    (out / 'methods_submission.md').write_text('\n'.join(lines))


def methods(out, rows, evidence, manifest, checks, selection, comparator):
    automated = selection.get('selection_mode') == 'automated_ordinal'
    lines = ['# Track 1 v3 — methods and local findings', '',
             ('**Automated local submission draft.** ' + AUTOMATED_POLICY if automated else
              '**Local review package.** The pipeline has not assigned EPCR probabilities or uploaded a submission. '
              'Review the evidence and complete participant disclosure before release.'), '',
             '## Participant and AI disclosure', '',
             *participant_lines(selection.get('identity', {})),
             '- OpenAI Codex assisted with code, public literature and documentation. The participant-supplied disclosure must state the actual plan/tier and relevant controls; account settings have not been independently verified.',
             '- The participant reports training opt-out was enabled long before this work; this was not independently inspected.',
             '- Earlier conversation content included patient-specific information before the local-only restriction. Subsequent v3 processing and results remained in local programs/files; the hosted assistant reviewed source and fabricated tests, not these findings.',
             '- Earlier source defaults included case-specific input filenames visible during source review. These defaults were replaced with a generic resolver and excluded local configuration; the change does not remove earlier conversation records.',
             '- No model was trained or fine-tuned on this child by this pipeline. Published pretrained annotations are evidence inputs, not functional experiments.',
             ('- Case-level human review: not performed for this automated selection. '
              'No retrospective guarantee of earlier data handling or provider settings is made.' if automated else
              '- Case-level human review: [TO CONFIRM]. No retrospective guarantee of earlier data handling or provider settings is made.'), '',
             '## Approach', '', (ROOT / 'docs/TRACK1_V3.md').read_text(), '',
             '## Execution record', '', '```json', json.dumps(manifest, indent=2), '```', '',
             '## Candidate review queue', '',
             'This table is broader than the final submission. Coverage, mechanism and phase limitations must be considered alongside rank.', '',
             '| Rank | Gene / disease | Alleles | Utility | Phase | Review tier |', '|---:|---|---|---:|---|---|']
    if not manifest.get('exomiser_mvstore_sha256'):
        lines.insert(4, '**' + LIMITED_REFERENCE_NOTICE + '**\n')
    for row in rows[:20]:
        values = [row['rank'], row['gene'] + ' / ' + row['disease_id'], ' + '.join(row['members']),
                  f"{row['score']:.5f}", row['phase'], row.get('review_tier', 'review required')]
        lines.append('| ' + ' | '.join(map(cell, values)) + ' |')
    lines += ['', '## Evidence for the leading review hypotheses', '']
    for row in rows[:10]:
        lines += [f"### Rank {row['rank']}: {row['gene']}", '', '```json',
                  json.dumps({'hypothesis': row, 'alleles': [evidence[x] for x in row['members']]}, indent=2), '```', '']
    lines += ['## Technical validation', '', '```json', json.dumps(checks, indent=2), '```', '',
              'These checks establish allele/reference consistency, not biological causality. '
              'Synthetic tests establish implementation behavior, not patient-level accuracy.', '',
              '## Independent comparator and benchmark status', '', '```json',
              json.dumps(comparator, indent=2), '```', '',
              'Comparator agreement is supporting context and is not an independent diagnosis: tools can share reference evidence. '
              'Any synthetic demonstration is not an external clinical accuracy estimate.', '',
              '## Selection and confidence', '',
              (AUTOMATED_POLICY + ' The organizer states that relative ranking determines scoring; '
               'the formal EPCR label remains a probability label, so these proxy values must not be '
               'interpreted as diagnostic risk estimates. This is a transparent ranking submission, '
               'not a calibrated probability model. All selected rows are primary hypotheses; no '
               'unreviewed secondary finding is asserted. Exomiser agreement is reported separately '
               'and does not alter the v3 order.' if automated else
               'No constant confidence is inserted. The local reviewer chooses hypotheses and explains each EPCR value. '
               'Uncalibrated utilities must not be interpreted as probabilities.') +
              ' Selection is bound to this analysis by its file hash.', '',
              '```json', json.dumps(selection, indent=2), '```', '',
              '## Required acknowledgement', '', ACK, '',
              'Dataset citation: [TO CONFIRM FROM CURRENT ORGANIZER CITATION BEFORE PUBLICATION].']
    (out / 'methods.md').write_text('\n'.join(lines) + '\n')


def checked_comparator(status, hypotheses_path, expected_inputs):
    """Exclude stale runs or altered exports without exposing private values."""
    assessed = dict(status)
    reasons = []
    if not status.get('completed'):
        reasons.append('comparator_not_completed')
    if status.get('provenance_schema_version') != 1:
        reasons.append('missing_current_provenance')
    recorded = status.get('analysis_input_hashes', {})
    if any(not value or recorded.get(key) != value for key, value in expected_inputs.items()):
        reasons.append('input_hash_mismatch')
    if status.get('phenotype_input_kind') != 'docx':
        reasons.append('phenotype_source_not_directly_comparable')
    if not status.get('inputs_unchanged_during_run') or not status.get('reference_size_and_hash_verified'):
        reasons.append('incomplete_input_or_reference_verification')
    if not hypotheses_path.is_file() or digest(hypotheses_path) != status.get('output_hashes', {}).get('hypotheses'):
        reasons.append('missing_or_changed_comparator_export')
    assessed['eligible_for_current_comparison'] = not reasons
    assessed['comparison_exclusion_reasons'] = reasons
    return assessed, json.loads(hypotheses_path.read_text()) if not reasons else None


def validate_analysis_artifacts(analysis, current_inputs, fasta_path):
    """Only package a completed run whose complete evidence set is unchanged."""
    manifest_path = analysis / 'run_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('provenance_schema_version') != 1 or manifest.get('completed') is not True:
        raise ValueError('Main analysis is incomplete or lacks current provenance; rerun before packaging.')
    if (manifest.get('input_hashes', {}).get('vcf') != current_inputs['vcf']
            or manifest.get('phenotype_sha256') != current_inputs['phenotype']
            or manifest.get('input_hashes', {}).get('fasta') != digest(fasta_path)):
        raise ValueError('Main analysis input provenance is stale; rerun before packaging.')
    expected_files = {name + '.json' for name in (
        'candidate_evidence', 'ranked_hypotheses', 'review_hypotheses',
        'disease_phenotype', 'sensitivity', 'review_counts', 'annotation_counts')}
    output_hashes = manifest.get('output_hashes')
    if not isinstance(output_hashes, dict) or set(output_hashes) != expected_files:
        raise ValueError('Main analysis output provenance is missing or incomplete.')
    for name, expected in output_hashes.items():
        path = analysis / name
        if not expected or not path.is_file() or digest(path) != expected:
            raise ValueError('Main analysis evidence changed or is incomplete; rerun before packaging.')
    return manifest


def archive_previous_exports(out):
    """Keep prior reviewed files locally, outside the active submission paths."""
    previous = [out / name for name in ('predictions.csv', 'selection.reviewed.json', 'selection.automated.json',
                                       'methods_submission.md', 'submission_bundle.zip', 'submission_readiness.json')
                if (out / name).exists()]
    archived = []
    if previous:
        directory = out / 'superseded_exports' / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        directory.mkdir(parents=True)
        for path in previous:
            target = directory / path.name
            path.replace(target)
            archived.append(str(target.relative_to(out)))
    return archived


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis', type=Path, default=ROOT / 'work/private_analysis.nosync/track1_v3')
    parser.add_argument('--output', type=Path, default=ROOT / 'work/private_analysis.nosync/submission_track1_v3')
    parser.add_argument('--comparator-output', type=Path, default=ROOT / 'work/private_analysis.nosync/track1_exomiser')
    parser.add_argument('--vcf', type=Path, default=input_path('vcf'))
    parser.add_argument('--fasta', type=Path, default=ROOT / 'references/track1_v2/hg38.fa')
    parser.add_argument('--phenotype', type=Path, default=input_path('phenotype'))
    parser.add_argument('--identity', type=Path, help='Local JSON containing team, github_url and accurate ai_disclosure.')
    export = parser.add_mutually_exclusive_group()
    export.add_argument('--selection', type=Path, help='Explicit locally reviewed manifest.')
    export.add_argument('--automated', action='store_true',
                        help='Export model-ranked hypotheses with disclosed ordinal EPCR weights, not probabilities.')
    args = parser.parse_args()
    install_network_guard()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / 'package_status.json', {
        'completed': False, 'review_package_created': False,
        'reviewed_csv_exported_this_run': False, 'submission_status': 'not_submitted',
    })
    # Read a supplied selection before moving a prior active export: callers
    # may intentionally finalize using this directory's selection.reviewed.json.
    supplied_selection = args.selection.read_text() if args.selection else None
    identity_path = args.identity or out / 'entry_identity.json'
    supplied_identity = json.loads(identity_path.read_text()) if identity_path.exists() else {}
    identity, identity_missing = validate_identity(supplied_identity)
    archived_exports = archive_previous_exports(out)
    current_inputs = {'vcf': digest(args.vcf), 'phenotype': digest(args.phenotype)}
    manifest = validate_analysis_artifacts(args.analysis, current_inputs, args.fasta)
    ranked_path = args.analysis / 'ranked_hypotheses.json'
    rows = json.loads(ranked_path.read_text())
    if not rows:
        raise ValueError('No hypotheses were ranked. Do not invent a submission.')
    by_id = {v['id']: v for v in json.loads((args.analysis / 'candidate_evidence.json').read_text())}
    # Bind review approval to all recorded evidence and provenance, including
    # evidence updates that happen to leave the ranked allele list unchanged.
    analysis_hash = digest(args.analysis / 'run_manifest.json')
    review_rows = distinct_review_rows(rows)
    checks = verify_calls(review_rows, by_id, args.vcf, args.fasta)
    selection = review_manifest(review_rows, analysis_hash)
    write_json(out / 'selection.template.json', selection)
    notice = LIMITED_REFERENCE_NOTICE if not manifest.get('exomiser_mvstore_sha256') else ''
    if args.automated:
        notice += (' Automated predictions.csv is prepared by the model; this page is an optional '
                   'human-reviewed alternative. Automated EPCR numbers encode rank only, not probability.')
    create_review_html(out / 'index.html', review_rows, by_id, selection, notice)
    previous_path = ROOT / 'outputs/track1_v2/ranked_hypotheses.json'
    previous = json.loads(previous_path.read_text()) if previous_path.exists() else []
    with (out / 'ranked_review.csv').open('w', newline='') as handle:
        fields = ['rank', 'gene', 'disease_id', 'inheritance', 'members', 'score', 'phase', 'review_tier', 'flags']
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in review_rows:
            writer.writerow({key: '; '.join(row[key]) if isinstance(row.get(key), list) else row.get(key, '') for key in fields})
    comparator_path = args.comparator_output / 'comparator_status.json'
    comparator = json.loads(comparator_path.read_text()) if comparator_path.exists() else {
        'status': 'not_completed', 'accuracy_estimate': None,
        'explanation': 'See public comparator setup report; no completed comparison is claimed.'}
    comparator_results_path = comparator_path.with_name('comparator_hypotheses.json')
    comparator, comparator_hypotheses = checked_comparator(comparator, comparator_results_path, current_inputs)
    write_comparison(out / 'comparison.md', rows, previous, comparator_hypotheses, args.fasta)
    exported = False
    selected = []
    if args.automated:
        selected, selection = automated_choices(rows, analysis_hash, hypothesis_id)
    if args.selection:
        selection = json.loads(supplied_selection)
        selected = validate_choices(selection, rows, analysis_hash)
        if not supplied_identity:
            identity, identity_missing = validate_identity(selection.get('identity', {}))
    selection['identity'] = identity
    write_json(out / 'entry_identity.json', identity)
    if args.automated or args.selection:
        verify_calls([row for row, _ in selected], by_id, args.vcf, args.fasta)
        target = out / 'predictions.csv'
        with target.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS)
            writer.writeheader()
            for row, choice in selected:
                writer.writerow(csv_row(row, choice, by_id))
        validation = validate_csv(target)
        if not validation['valid']:
            target.unlink()
            raise ValueError('Prediction CSV failed schema validation.')
        write_json(out / ('selection.automated.json' if args.automated else 'selection.reviewed.json'), selection)
        checks['csv_validation'] = validation
        exported = True
    methods(out, rows, by_id, manifest, checks, selection, comparator)
    if exported:
        submission_report(out, selected, manifest, comparator, selection)
        readiness = validate_package(out / 'predictions.csv', out / 'methods_submission.md',
                                     identity.get('github_url') or None, require_package=True)
        write_json(out / 'submission_readiness.json', readiness)
    write_json(out / 'validation_evidence.json', checks)
    package = source_archive(out)
    if exported:
        with zipfile.ZipFile(out / 'submission_bundle.zip', 'w', zipfile.ZIP_DEFLATED) as bundle:
            for name in ('predictions.csv', 'methods_submission.md', 'repository.zip'):
                bundle.write(out / name, name)
    (out / 'HOW_TO_SUBMIT.md').write_text(
        '# Submit Track 1\n\n'
        'The completed local analysis is already available. Do not download or rerun the genome to finish the entry.\n\n'
        '1. Fill `entry_identity.json` with the team name, real GitHub repository URL, and accurate AI-use disclosure '
        '(provider, actual plan/tier, relevant data-handling control and known limitations). Do not invent unknown facts.\n'
        '2. Publish only the prepared `repository/` source folder or the contents of `repository.zip` to that GitHub repository. '
        'Do not publish this private package, the whole workspace, raw genomic files, logs, or evidence tables. '
        'GitHub authentication and the destination must be supplied by the participant.\n'
        '3. Rebuild the entry files to apply those details:\n\n'
        '```bash\n.venv/bin/python scripts/build_track1_v3_submission.py --automated '
        '> work/private_analysis.nosync/track1_v3/finalize.log 2>&1\n```\n\n'
        '4. Read `submission_readiness.json` locally. The automated CSV preserves model order; its EPCR values are '
        'ordinal weights, not probabilities. A human-reviewed alternative remains available through `index.html`.\n'
        '5. If releasing these limited findings is authorized, sign in at '
        'https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026 and open Submit Track 1. '
        'Rename a copy of `predictions.csv` to include your HF username or team plus `v3-automated`. '
        'Upload that CSV and `methods_submission.md`, and enter the GitHub link. The ZIP is a local convenience bundle, '
        'not a substitute for the separate fields in the submission form.\n'
        '6. Check the files and submit once. The challenge allows six Track 1 submissions; none has been used by this workflow. '
        'Keep the returned score/receipt locally. Do not paste patient findings or upload these files to an AI assistant.\n\n'
        'The methods include required acknowledgements and resource links. At external publication, verify the organizer-current '
        'formal dataset citation on the Synapse project. No remote repository or challenge action is performed by these scripts.\n')
    write_json(out / 'package_status.json', {'created_utc': datetime.now(timezone.utc).isoformat(),
               'completed': True, 'archived_previous_exports': archived_exports,
               'analysis_sha256': analysis_hash, 'review_package_created': True,
               'reviewed_csv_exported_this_run': exported and not args.automated,
               'automated_csv_exported_this_run': exported and args.automated,
               'prediction_csv_exported_this_run': exported,
               'selection_mode': 'automated_ordinal' if args.automated else ('reviewed' if args.selection else 'review_only'),
               'submission_status': 'not_submitted',
               'full_population_reference_used': bool(manifest.get('exomiser_mvstore_sha256')),
               'current_comparator_eligible': comparator.get('eligible_for_current_comparison', False),
               'identity_and_disclosure_complete': not identity_missing,
               'participant_fields_missing': identity_missing,
               'provider_settings_independently_verified': False,
               'github_remote_verified': False, 'source_archive': package})
    automated_intro = ('The automated [predictions CSV](predictions.csv) is ready locally. '
                       'It contains up to ten distinct complete allele hypotheses in model order. '
                       'Its EPCR numbers are ordinal weights, not causal probabilities. '
                       'Use [the concise submission methods](methods_submission.md); '
                       '[the bundle](submission_bundle.zip) contains the CSV, concise report and source archive. '
                       'See [the automated selection record](selection.automated.json) for provenance. '
                       'No case-level human review or upload is claimed.\n\n') if args.automated else ''
    (out / 'START_HERE.md').write_text('# Track 1 v3 local package\n\n' + (notice + '\n\n' if notice else '') + automated_intro +
        '1. Open [the private review page](index.html) locally. Read [methods](methods.md) and [the v2 comparison](comparison.md).\n'
        '2. An optional human-reviewed alternative can be exported from the review page. Review both alleles, disease mechanism, population evidence and unresolved phase, and enter your confidence with reasons.\n'
        '3. Save the reviewed selection and CSV. For independent Python validation, run the command below with the path to the saved selection.\n'
        '4. Complete team, GitHub and accurate AI-use disclosure in `entry_identity.json`, then rerun the finalizer with `--automated` (or the reviewed selection) to apply them to both methods reports. The source-only [repository archive](repository.zip) contains no patient results.\n'
        '5. Follow [the submission instructions](HOW_TO_SUBMIT.md) when releasing the small findings and methods is authorized. Nothing has been published or submitted by this workflow.\n\n'
        '```bash\n.venv/bin/python scripts/build_track1_v3_submission.py --selection /path/to/selection.reviewed.json > work/private_analysis.nosync/track1_v3/finalize.log 2>&1\n```\n\n'
        'The ranked_review.csv file is a review queue, not an official predictions CSV. The original v2 draft remains a separate historical artifact. '
        'Do not upload this HTML, private evidence tables, logs, genomic files or the whole workspace to an AI tool or a repository.\n')
    print('Private v3 submission package prepared. No upload was made.' if exported else
          'Private v3 review package prepared. No upload was made.')


if __name__ == '__main__':
    main()
