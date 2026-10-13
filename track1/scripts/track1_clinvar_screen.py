"""Offline genome-wide observed-allele intersection with the public ClinVar VCF.

This is an independent safety-net screen, not an ACMG classifier. It does not
infer that an exact ClinVar match explains this patient. No network is used.
"""
from __future__ import annotations

if __package__:
    from .local_inputs import input_path
else:
    from local_inputs import input_path

import argparse
from collections import Counter
from datetime import datetime, timezone
import heapq
import itertools
import json
from pathlib import Path
import re
import sys
import time

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
from track1_baseline import (
    clinical_support, digest, minimal_allele, PhenotypeMatcher,
    rarity_evidence, serial, text_field,
)

ROOT = Path(__file__).resolve().parents[1]
CV_FIELDS = ['CLNSIG', 'CLNSIGCONF', 'CLNREVSTAT', 'CLNDN', 'CLNDISDB',
             'CLNHGVS', 'GENEINFO', 'MC', 'AF_EXAC', 'AF_ESP', 'AF_TGP',
             'ALLELEID', 'CLNSIGINCL', 'CLNDNINCL']


def canonical_contig(name):
    name = name.removeprefix('chr')
    return 'MT' if name in {'M', 'MT'} else name


def normalized_alleles(records, sample_name=None, counts=None):
    """Yield minimally trimmed sequence alleles in normalized coordinate order.

    Prefix trimming may move an allele beyond the next VCF record. A small heap
    delays emission until the next original position guarantees ordering. This
    is minimal representation only, not repeat left-alignment. Patient ALTs
    must actually be called in GT; ClinVar ALTs are reference annotations.
    """
    counts = counts if counts is not None else Counter()
    pending = []
    sequence_number = 0
    previous_pos = 0
    for rec in records:
        counts['records'] += 1
        if rec.pos < previous_pos:
            raise ValueError('Input records are not sorted within a contig')
        previous_pos = rec.pos
        while pending and pending[0][0][0] < rec.pos:
            key, _, saved, alt_index = heapq.heappop(pending)
            yield key, saved, alt_index
        if sample_name is not None:
            gt = rec.samples[sample_name].get('GT')
            if not gt or any(x is None for x in gt):
                counts['missing_or_partial_GT_records'] += 1
                continue
            indices = sorted(set(gt) - {0})
        else:
            indices = range(1, len(rec.alts or ()) + 1)
        for alt_index in indices:
            if not rec.alts or alt_index < 1 or alt_index > len(rec.alts):
                counts['invalid_ALT_index'] += 1
                continue
            ref, alt = rec.ref.upper(), rec.alts[alt_index - 1].upper()
            if not ref or not alt or not set(ref + alt) <= set('ACGT'):
                counts['nonsequence_alleles_excluded'] += 1
                continue
            key = minimal_allele(rec.pos, ref, alt)
            heapq.heappush(pending, (key, sequence_number, rec, alt_index))
            sequence_number += 1
            counts['sequence_alleles'] += 1
    while pending:
        key, _, rec, alt_index = heapq.heappop(pending)
        yield key, rec, alt_index


def grouped_alleles(alleles):
    for key, group in itertools.groupby(alleles, key=lambda item: item[0]):
        yield key, [(record, alt_index) for _, record, alt_index in group]


def intersect_alleles(patient, reference):
    """Merge two sorted allele iterators; match position AND REF AND ALT."""
    patient, reference = iter(grouped_alleles(patient)), iter(grouped_alleles(reference))
    a, b = next(patient, None), next(reference, None)
    while a is not None and b is not None:
        if a[0] < b[0]:
            a = next(patient, None)
        elif b[0] < a[0]:
            b = next(reference, None)
        else:
            yield a[0], a[1], b[1]
            a, b = next(patient, None), next(reference, None)
    # Fully exhaust both sides so all genome-wide input counters are complete.
    for _ in patient:
        pass
    for _ in reference:
        pass


def clinvar_payload(record):
    return {'id': record.id, **{k: serial(record.info[k]) for k in CV_FIELDS if k in record.info}}


def has_pathogenic_conflict(cv):
    # The word "pathogenicity" in the aggregate conflict label is not itself a
    # pathogenic assertion: many conflicts are benign versus uncertain.
    assertions = text_field(cv.get('CLNSIGCONF')).lower()
    return bool(re.search(r'(?<![a-z])(?:likely_)?pathogenic(?:\(|[|,/]|$)', assertions))


def patient_payload(rec, alt_index, sample_name, minimum_dp, minimum_gq):
    sample = rec.samples[sample_name]
    gt, dp, gq, ad = sample.get('GT'), sample.get('DP'), sample.get('GQ'), sample.get('AD')
    flags = []
    if list(rec.filter) != ['PASS']:
        flags.append('FILTER_not_PASS')
    if dp is None or dp < minimum_dp:
        flags.append('DP_missing_or_low')
    if gq is None or gq < minimum_gq:
        flags.append('GQ_missing_or_low')
    allele_balance = (ad[alt_index] / sum(ad)
                     if ad and len(ad) > alt_index and all(v is not None for v in ad) and sum(ad) > 0 else None)
    dosage = gt.count(alt_index)
    genotype = 'hom_alt' if len(gt) == 2 and dosage == 2 else 'het' if len(gt) == 2 else 'haploid_or_other'
    return {'chrom': 'chr' + canonical_contig(rec.contig).replace('MT', 'M'),
            'pos': rec.pos, 'ref': rec.ref, 'alt': rec.alts[alt_index - 1],
            'alt_index': alt_index, 'gt': list(gt), 'ploidy': len(gt),
            'alt_dosage': dosage, 'genotype': genotype, 'dp': dp, 'gq': gq,
            'ad': serial(ad), 'allele_balance': allele_balance, 'filter': list(rec.filter),
            'quality_pass': not flags, 'quality_flags': flags}


def screen(vcf_path, clinvar_path, sample_name=None, minimum_dp=10, minimum_gq=20):
    patient_counts, reference_counts, matches = Counter(), Counter(), Counter()
    by_contig, prioritized, unrepresented = {}, [], []
    with pysam.VariantFile(str(vcf_path)) as patient, pysam.VariantFile(str(clinvar_path)) as reference:
        if 'GRCh38' not in str(patient.header):
            raise ValueError('Patient VCF header does not establish GRCh38')
        samples = list(patient.header.samples)
        if sample_name is None:
            if len(samples) != 1:
                raise ValueError('Select one patient sample explicitly')
            sample_name = samples[0]
        if sample_name not in samples:
            raise ValueError('Requested sample is absent from VCF')
        reference_contigs = {canonical_contig(c): c for c in reference.header.contigs}
        for contig in patient.header.contigs:
            canonical = canonical_contig(contig)
            pc, rc, mc = Counter(), Counter(), Counter()
            called = normalized_alleles(patient.fetch(contig), sample_name, pc)
            if canonical not in reference_contigs:
                for _ in called:
                    pass
                if pc['records']:
                    unrepresented.append({'contig': contig, 'records': pc['records']})
            else:
                cv = normalized_alleles(reference.fetch(reference_contigs[canonical]), None, rc)
                for key, observed, cv_records in intersect_alleles(called, cv):
                    mc['distinct_minimal_alleles_matched'] += 1
                    cvs = [clinvar_payload(r) for r, _ in cv_records]
                    supports = [clinical_support(c) for c in cvs]
                    any_plp = any(c['pathogenic_assertion'] for c in supports)
                    any_conflict = any(c['conflicting'] for c in supports)
                    cross_record_conflict = any_plp and any(c['benign'] for c in supports)
                    conflicting_plp = any(has_pathogenic_conflict(c) for c in cvs) or cross_record_conflict
                    for label in {c['classification'] for c in supports}:
                        mc['classification:' + label] += 1
                    if not any_plp and not conflicting_plp:
                        continue
                    reviewed_plp = any(c['pathogenic_assertion'] and c['review_stars'] >= 1 for c in supports)
                    unconflicted_plp = any_plp and not any_conflict and not cross_record_conflict
                    priority_class = ('unconflicted_PLP_at_least_one_star' if unconflicted_plp and reviewed_plp
                                      else 'unconflicted_PLP_zero_star' if unconflicted_plp
                                      else 'conflicting_pathogenic_evidence_review_only')
                    mc[priority_class] += 1
                    for rec, alt_index in observed:
                        p = patient_payload(rec, alt_index, sample_name, minimum_dp, minimum_gq)
                        genes = sorted({part.split(':')[0] for c in cvs
                                        for part in re.split(r'[|,]', text_field(c.get('GENEINFO'))) if ':' in part})
                        frequencies = [rarity_evidence(c) for c in cvs]
                        known_afs = [r['af'] for r in frequencies if r['known']]
                        flags = list(p['quality_flags'])
                        if any_conflict or cross_record_conflict:
                            flags.append('ClinVar_conflicting_evidence')
                        if not known_afs:
                            flags.append('population_frequency_unknown')
                        elif max(known_afs) > .01:
                            flags.append('legacy_population_AF_over_one_percent')
                        if p['allele_balance'] is not None and p['genotype'] == 'het' and not .2 <= p['allele_balance'] <= .8:
                            flags.append('heterozygous_allele_balance_outside_0.2_to_0.8')
                        if any('low_penetrance' in text_field(c.get('CLNSIG')).lower() or 'risk' in text_field(c.get('CLNSIG')).lower() for c in cvs):
                            flags.append('risk_or_low_penetrance_assertion_requires_context')
                        p.update(id=f"{p['chrom']}:{p['pos']}:{p['ref']}:{p['alt']}",
                                 minimal_allele={'pos': key[0], 'ref': key[1], 'alt': key[2]},
                                 genes=genes, clinical_priority_class=priority_class,
                                 clinvar_records=cvs, clinical_evidence=supports,
                                 maximum_known_legacy_population_af=max(known_afs) if known_afs else None,
                                 flags=sorted(set(flags)),
                                 interpretation='Exact public annotation match; condition, inheritance and mechanism require review; no patient-level causal claim')
                        prioritized.append(p)
                        if p['quality_pass']:
                            mc[priority_class + '_quality_pass_observed_calls'] += 1
            patient_counts.update(pc)
            reference_counts.update(rc)
            matches.update(mc)
            if pc['records']:
                by_contig[contig] = {'patient': dict(pc), 'clinvar': dict(rc), 'matches': dict(mc)}
    return {'sample': sample_name, 'patient_counts': dict(patient_counts),
            'clinvar_counts_in_intersecting_contigs': dict(reference_counts),
            'match_counts': dict(matches), 'by_contig': by_contig,
            'patient_contigs_absent_from_clinvar': unrepresented,
            'prioritized_hits': prioritized}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--vcf', type=Path, default=input_path('vcf'))
    p.add_argument('--references', type=Path, default=ROOT / 'references/track1')
    p.add_argument('--phenotypes', type=Path, default=ROOT / 'outputs/track1_baseline/phenotypes.json')
    p.add_argument('--output', type=Path, default=ROOT / 'outputs/track1_v2/clinvar_screen.json')
    p.add_argument('--sample')
    p.add_argument('--minimum-dp', type=int, default=10)
    p.add_argument('--minimum-gq', type=int, default=20)
    a = p.parse_args()
    start = time.perf_counter()
    cv_path = a.references / 'clinvar.vcf.gz'
    manifest = json.loads((a.references / 'manifest.json').read_text())
    hashes = {}
    for name in ['clinvar.vcf.gz', 'clinvar.vcf.gz.tbi', 'hp.obo', 'genes_to_phenotype.txt']:
        hashes[name] = digest(a.references / name)
        if hashes[name] != manifest[name]:
            raise ValueError('Reference hash mismatch: ' + name)
    result = screen(a.vcf, cv_path, a.sample, a.minimum_dp, a.minimum_gq)
    matcher = PhenotypeMatcher(a.references / 'hp.obo', a.references / 'genes_to_phenotype.txt', {'CCDC84': 'CENATAC'})
    phenotypes = json.loads(a.phenotypes.read_text())
    query = [r['hpo_id'] for r in phenotypes if r['subject'] == 'proband']
    for hit in result['prioritized_hits']:
        hit['gene_phenotype_evidence'] = {g: matcher.score(g, query) for g in hit['genes']}
        hit['maximum_gene_phenotype_score'] = max((x['score'] for x in hit['gene_phenotype_evidence'].values()), default=0)
    result['prioritized_hits'].sort(key=lambda x: (not x['quality_pass'],
        x['clinical_priority_class'] != 'unconflicted_PLP_at_least_one_star',
        -x['maximum_gene_phenotype_score'], x['id']))
    result['provenance'] = {'run_utc': datetime.now(timezone.utc).isoformat(),
        'model': 'genomewide_exact_ClinVar_observed_allele_screen_v1',
        'code_sha256': digest(__file__), 'baseline_code_sha256': digest(ROOT / 'scripts/track1_baseline.py'),
        'vcf_sha256': digest(a.vcf),
        'vcf_index_sha256': digest(str(a.vcf) + '.tbi'),
        'reference_hashes': hashes, 'phenotypes_sha256': digest(a.phenotypes),
        'pysam_version': pysam.__version__, 'network_during_analysis': False,
        'minimum_dp': a.minimum_dp, 'minimum_gq': a.minimum_gq,
        'elapsed_seconds': time.perf_counter() - start}
    result['limitations'] = [
        'A clinically classified allele may be a recessive carrier finding or unrelated to the child; this is not a causal ranking.',
        'ClinVar aggregate labels mix conditions and review quality; condition-specific review remains essential.',
        'Only called sequence ALTs in complete GT are compared; missing or partial GT and symbolic structural alleles are excluded.',
        'Equivalent repeat-shifted indels can remain unmatched because minimal trimming is not full left-normalization.',
        'ClinVar absence is not benignity or rarity; novel causes and many population alleles are absent.',
        'Population AF uses legacy ClinVar fields only; no modern comprehensive ancestry-specific population data.',
        'HPO is gene-level annotation-union matching for triage, not a disease likelihood; inheritance is not scored.',
        'Input record counts should be compared with an independent sequential VCF audit to check indexed retrieval.',
    ]
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    write_report(a.output.with_suffix('.md'), result)
    print(json.dumps({'output': str(a.output), 'patient_records': result['patient_counts']['records'],
                      'matched_alleles': result['match_counts'].get('distinct_minimal_alleles_matched', 0),
                      'prioritized_hits': len(result['prioritized_hits']),
                      'elapsed_seconds': result['provenance']['elapsed_seconds']}, indent=2))


def write_report(path, result):
    """Write findings locally; callers can redirect execution logs as needed."""
    hits = result['prioritized_hits']
    reviewed = [x for x in hits if x['quality_pass'] and
                x['clinical_priority_class'] == 'unconflicted_PLP_at_least_one_star']
    text = ['# Local genome wide ClinVar screen', '',
        'This report and its JSON contain patient-derived findings and must remain local.', '',
        'This independent screen tests whether called patient alleles exactly match minimally trimmed public ClinVar alleles. '
        'It is a safety net for established annotated variants outside the selected transcript panel, not a genome-wide diagnostic model. '
        'Every hit still requires review of the condition, inheritance, mechanism, genotype and phenotype.', '',
        f"Patient VCF records visited: {result['patient_counts'].get('records', 0):,}.",
        f"Distinct minimally represented alleles matched: {result['match_counts'].get('distinct_minimal_alleles_matched', 0):,}.",
        f'Prioritized observed calls: {len(hits):,}.',
        f'Quality-passing unconflicted pathogenic or likely pathogenic calls with at least one review star: {len(reviewed):,}.', '',
        'The listed order is triage by quality, review tier and gene-level HPO compatibility. '
        'A high HPO score is not disease probability. A heterozygous recessive allele may be a carrier finding. '
        'No causal interpretation or official ACMG classification is assigned here.', '',
        '## Reviewed quality passing annotations', '',
        '| Observed allele | Genes | GT | DP | GQ | Maximum gene HPO score | ClinVar classification | Conditions |',
        '|---|---|---|---:|---:|---:|---|---|']
    def cell(value):
        return str(value).replace('|', '; ').replace('\n', ' ')
    for hit in reviewed:
        classifications = '; '.join(sorted({c['classification'] for c in hit['clinical_evidence']}))
        conditions = '; '.join(sorted({text_field(c.get('CLNDN')) for c in hit['clinvar_records']}))
        text.append('| ' + ' | '.join(cell(v) for v in [hit['id'], ', '.join(hit['genes']),
            hit['gt'], hit['dp'], hit['gq'], f"{hit['maximum_gene_phenotype_score']:.4f}", classifications, conditions]) + ' |')
    text.extend(['', '## Limitations', ''])
    text.extend('- ' + value for value in result['limitations'])
    text.extend(['', '## Provenance', '',
        'The companion JSON retains all prioritized calls, exact ClinVar records, quality and frequency flags, '
        'per-gene phenotype evidence, contig-level counts, reference hashes and execution provenance. '
        'Inputs were processed locally without network annotation.'])
    path.write_text('\n'.join(text) + '\n')


if __name__ == '__main__':
    main()
