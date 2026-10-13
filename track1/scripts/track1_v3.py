"""Offline disease-aware ranking with explicit evidence coverage and review gates.

No patient files are opened on import. Public downloads are a separate operation.
The command writes patient-derived results locally and prints generic status only.
"""
from __future__ import annotations

if __package__:
    from .local_inputs import input_path
else:
    from local_inputs import input_path

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import itertools
import json
import math
from pathlib import Path
import re
import sqlite3
import time

import pysam
import pysam.bcftools

from track1_baseline import clinical_support, digest, extract_phenotypes, minimal_allele, serial, write_json
from track1_v2 import DiseaseMatcher, read_csq, install_network_guard, local_path, CV_FIELDS

ROOT = Path(__file__).resolve().parents[1]
LOF = {'stop_gained', 'frameshift', 'splice_donor', 'splice_acceptor', 'start_lost'}


def canonical_chrom(chrom):
    name = str(chrom).removeprefix('chr')
    return 'M' if name == 'MT' else name


def positive_proband_terms(phenotypes):
    """Repeated mentions are not independent evidence; absent/family terms differ."""
    return sorted({entry['hpo_id'] for entry in phenotypes
                   if entry.get('subject') == 'proband' and entry.get('present', True)})


def normalize_condition_id(value):
    """Normalize identifier spelling without inferring cross-database equivalence."""
    value = str(value).strip()
    value = re.sub(r'^(?:Orphanet|ORPHANET):', 'ORPHA:', value)
    value = re.sub(r'^MONDO:MONDO:', 'MONDO:', value)
    return value


def condition_ids(value):
    if isinstance(value, (tuple, list, set)):
        return set().union(*(condition_ids(v) for v in value)) if value else set()
    text = str(value or '')
    text = re.sub(r'MONDO:MONDO:', 'MONDO:', text)
    text = text.replace('MedGen:', 'MEDGEN:')
    return {normalize_condition_id(v) for v in re.findall(r'(?:OMIM|ORPHA|Orphanet|MONDO|MEDGEN):[A-Za-z0-9_]+', text)}


def fasta_contig(fasta, chrom):
    wanted = canonical_chrom(chrom)
    candidates = [chrom, wanted, 'chr' + wanted]
    if wanted == 'M':
        candidates += ['MT', 'chrMT']
    names = set(fasta.references)
    return next((c for c in candidates if c in names), None)


def normalize_allele(chrom, pos, ref, alt, fasta=None):
    """Validate REF and left-align ordinary sequence indels when FASTA is given.

    Coordinates are one-based. This does not interpret symbolic or breakend alleles.
    It retains a VCF anchor and shifts repeat-equivalent alleles toward the left.
    """
    chrom, pos, ref, alt = canonical_chrom(chrom), int(pos), ref.upper(), alt.upper()
    if pos < 1 or not ref or not alt or ref == alt or not set(ref + alt) <= set('ACGT'):
        raise ValueError('Only positive coordinates and ordinary sequence alleles are supported')
    contig = fasta_contig(fasta, chrom) if fasta is not None else None
    if fasta is not None:
        if contig is None:
            raise ValueError('Reference chromosome unavailable')
        if fasta.fetch(contig, pos - 1, pos - 1 + len(ref)).upper() != ref:
            raise ValueError('Reference allele mismatch during normalization')
    pos, ref, alt = minimal_allele(pos, ref, alt)
    if fasta is not None and len(ref) != len(alt):
        while pos > 1 and ref[-1] == alt[-1]:
            previous = fasta.fetch(contig, pos - 2, pos - 1).upper()
            if len(previous) != 1 or previous not in 'ACGT':
                break
            ref, alt, pos = previous + ref[:-1], previous + alt[:-1], pos - 1
        pos, ref, alt = minimal_allele(pos, ref, alt)
    return chrom, pos, ref, alt


def extract_variants(path, config):
    """Preserve every original record until reference normalization can reconcile it.

    The v2 extractor coalesces minimally equivalent alleles before comparing
    genotypes, which loses duplicate-record conflicts and original representations.
    """
    variants, counts = [], Counter()
    with pysam.VariantFile(str(path)) as source:
        if len(source.header.samples) != 1:
            raise ValueError('Single explicitly selected sample required')
        if 'GRCh38' not in str(source.header):
            raise ValueError('GRCh38 header required')
        sample_name = next(iter(source.header.samples))
        for record_number, rec in enumerate(source):
            counts['records'] += 1
            sample = rec.samples[sample_name]
            gt = sample.get('GT')
            if not gt or None in gt:
                continue
            if len(gt) not in {1, 2}:
                counts['unsupported_ploidy_records_excluded'] += 1
                continue
            for alt_index in sorted(set(gt) - {0}):
                if alt_index < 1 or not rec.alts or alt_index > len(rec.alts):
                    raise ValueError('Genotype references an unavailable alternate allele')
                alt = rec.alts[alt_index - 1]
                if not set(rec.ref + alt) <= set('ACGT') or rec.ref == alt:
                    continue
                annotations = read_csq(rec, alt, config)
                if not annotations:
                    continue
                dp, gq, ad = sample.get('DP'), sample.get('GQ'), sample.get('AD')
                balance = ad[alt_index] / sum(ad) if ad and len(ad) > alt_index and all(x is not None for x in ad) and sum(ad) > 0 else None
                genotype = 'hom_alt' if len(gt) == 2 and gt.count(alt_index) == 2 else 'het' if len(gt) == 2 else 'hemizygous'
                flags = []
                if list(rec.filter) != ['PASS']:
                    flags.append('FILTER_not_PASS')
                if dp is None or dp < config['minimum_dp']:
                    flags.append('low_or_missing_DP')
                if gq is None or gq < config['minimum_gq']:
                    flags.append('low_or_missing_GQ')
                if genotype == 'het' and balance is not None and not config['minimum_het_balance'] <= balance <= config['maximum_het_balance']:
                    flags.append('heterozygous_allele_balance_outside_configured_range')
                phase_set = haplotype = None
                if genotype == 'het':
                    if sample.phased and sample.get('PS') is not None:
                        phase_set, haplotype = f'{canonical_chrom(rec.contig)}:PS:{sample["PS"]}', gt.index(alt_index)
                    elif sample.get('PID') and re.fullmatch(r'\d+\|\d+', sample.get('PGT') or ''):
                        pgt = tuple(map(int, sample['PGT'].split('|')))
                        if sorted(pgt) == sorted(gt):
                            phase_set, haplotype = f'{canonical_chrom(rec.contig)}:PID:{sample["PID"]}', pgt.index(alt_index)
                chrom = 'chr' + canonical_chrom(rec.contig)
                variants.append({'id': f'{chrom}:{rec.pos}:{rec.ref}:{alt}', 'chrom': chrom,
                    'source_contig': rec.contig, 'pos': rec.pos, 'ref': rec.ref, 'alt': alt, 'alt_index': alt_index,
                    'record_id': f'{rec.contig}:{rec.pos}:{rec.ref}:{",".join(rec.alts)}',
                    'source_record_number': record_number, 'gt': list(gt), 'genotype': genotype,
                    'dp': dp, 'gq': gq, 'ad': serial(ad), 'allele_balance': balance,
                    'phase_set': phase_set, 'haplotype': haplotype, 'quality_pass': not flags,
                    'quality_flags': flags, 'annotations': annotations, 'clinvar_records': []})
    counts['functional_source_alleles'] = len(variants)
    return variants, dict(counts)


def allele_dosage(v):
    return len(v['gt']), v['gt'].count(v['alt_index'])


def phase_call(a, b):
    """Combine phase evidence from retained source records without choosing a conflict."""
    phases = set()
    for left in a.get('source_alleles', [a]):
        for right in b.get('source_alleles', [b]):
            same_record = left.get('source_record_number') == right.get('source_record_number')
            if (same_record and left.get('record_id') == right.get('record_id')
                    and left.get('alt_index') != right.get('alt_index')):
                expected = sorted([left['alt_index'], right['alt_index']])
                if sorted(left['gt']) == expected and sorted(right['gt']) == expected:
                    phases.add('trans')
            if (left.get('phase_set') and left.get('phase_set') == right.get('phase_set')
                    and left.get('haplotype') is not None and right.get('haplotype') is not None):
                phases.add('cis' if left['haplotype'] == right['haplotype'] else 'trans')
    if len(phases) > 1:
        return 'unknown', 'Conflicting source-record phase evidence; orthogonal phasing required'
    if phases:
        return phases.pop(), 'Consistent source-record genotype or shared VCF phase block; phase requires validation'
    return 'unknown', 'No shared informative VCF phase block; parental/orthogonal phasing required'


def normalized_variants(variants, fasta=None):
    result = {}
    counts = Counter()
    for v in variants.values() if isinstance(variants, dict) else variants:
        v = dict(v)
        for field in ['chrom', 'pos', 'ref', 'alt']:
            v['source_' + field] = v.get('source_' + field, v[field])
        v['source_contig'] = v.get('source_contig', v['source_chrom'])
        v['source_id'] = v.get('source_id', v['id'])
        try:
            key = normalize_allele(v['chrom'], v['pos'], v['ref'], v['alt'], fasta)
        except ValueError:
            counts['reference_mismatch_excluded'] += 1
            continue
        chrom, pos, ref, alt = key
        v.update(chrom='chr' + chrom, pos=pos, ref=ref, alt=alt,
                 id=f'chr{chrom}:{pos}:{ref}:{alt}', normalization='reference_left_aligned' if fasta else 'minimal_trim_only')
        source_fields = ['source_chrom', 'source_contig', 'source_pos', 'source_ref', 'source_alt', 'source_id',
                         'record_id', 'source_record_number', 'gt', 'alt_index', 'genotype', 'phase_set', 'haplotype',
                         'quality_pass', 'quality_flags', 'dp', 'gq', 'ad', 'allele_balance']
        v['source_alleles'] = [{k: v.get(k) for k in source_fields}]
        if key in result:
            counts['equivalent_source_alleles_merged'] += 1
            old = result[key]
            old['annotations'].extend(v['annotations'])
            old['source_alleles'].extend(v['source_alleles'])
            old['quality_pass'] = old['quality_pass'] and v['quality_pass']
            old['quality_flags'] = sorted(set(old['quality_flags'] + v['quality_flags']))
            if allele_dosage(old) != allele_dosage(v) or old.get('genotype') != v.get('genotype'):
                old['quality_pass'] = False
                old['quality_flags'] = sorted(set(old['quality_flags'] + ['equivalent_records_genotype_disagreement']))
            if old.get('phase_set') != v.get('phase_set') or old.get('haplotype') != v.get('haplotype'):
                old['phase_set'] = old['haplotype'] = None
        else:
            result[key] = v
    counts['normalized_functional_alleles'] = len(result)
    return result, dict(counts)


class PublicIndex:
    """Read-only SQLite adapter; no extension loading and no network operations."""
    TABLES = {'population': 'population', 'alphamissense': 'alphamissense_annotations', 'clinvar_rcv': 'rcv_assertions'}

    def __init__(self, path=None):
        self.connection = None
        self.tables = set()
        self.metadata = {}
        if path is not None:
            path = Path(path)
            if not path.is_file():
                raise FileNotFoundError('Public annotation index is absent')
            self.connection = sqlite3.connect(path.resolve().as_uri() + '?mode=ro', uri=True)
            self.connection.execute('PRAGMA query_only=ON')
            self.tables = {r[0] for r in self.connection.execute("SELECT name FROM sqlite_master WHERE type IN ('table','view')")}
            if 'metadata' in self.tables:
                for key, value in self.connection.execute('SELECT key,value FROM metadata'):
                    try:
                        self.metadata[key] = json.loads(value)
                    except (TypeError, ValueError):
                        self.metadata[key] = value
            # A partially built index must not silently look like absent evidence.
            # The population placeholder is deliberately incomplete when no
            # population channel has been imported and is not an interrupted build.
            for channel in ['alphamissense', 'clinvar_rcv', 'clinvar_scv', 'condition_aliases_provenance']:
                metadata = self.metadata.get(channel, {})
                if isinstance(metadata, dict) and metadata.get('complete') is False:
                    self.close()
                    raise ValueError('Public annotation index contains an incomplete source build')
                expected = 'clinvar_rcv' if channel == 'clinvar_scv' else channel
                if isinstance(metadata, dict) and metadata.get('complete') is True and channel != 'condition_aliases_provenance' and expected not in self.tables:
                    self.close()
                    raise ValueError('Public annotation index is missing a completed source view')

    def lookup(self, key):
        if self.connection is None:
            return {}
        result = {}
        if 'variants' in self.tables:
            row = self.connection.execute('SELECT payload FROM variants WHERE chrom=? AND pos=? AND ref=? AND alt=?', key).fetchone()
            if row:
                result.update(json.loads(row[0]))
        for table, channel in self.TABLES.items():
            if table in self.tables:
                row = self.connection.execute(f'SELECT payload FROM {table} WHERE chrom=? AND pos=? AND ref=? AND alt=?', key).fetchone()
                if row:
                    result[channel] = json.loads(row[0])
        return result

    def close(self):
        if self.connection:
            self.connection.close()


def population_evidence(v):
    public = v.get('public_annotation', {})
    rows = []
    population = public.get('population')
    if population:
        rows.extend(population if isinstance(population, list) else [population])
    for platform in ['exome', 'genome', 'joint']:
        record = (public.get('gnomad') or {}).get(platform)
        if record:
            rows.append({**record, 'source': 'gnomAD ' + platform})
    usable, rejected = [], []
    for row in rows:
        source = row.get('source', 'public population reference')
        filters = row.get('filters') or []
        if isinstance(filters, str):
            filters = [] if filters in {'PASS', '.'} else [filters]
        if filters or row.get('callable') is False:
            rejected.append(source)
            continue
        try:
            an = row.get('an', row.get('AN'))
            an = float(an) if an is not None else None
            if an is not None and (not math.isfinite(an) or an <= 0):
                continue
            ac = row.get('ac', row.get('AC'))
            af = row.get('af', row.get('AF'))
            if af is None and an and ac is not None:
                af = float(ac) / an
            # Unknown callability cannot establish absence: zero needs an AN denominator.
            af = float(af)
            if not math.isfinite(af) or not 0 <= af <= 1 or (af == 0 and not an):
                continue
        except (TypeError, ValueError):
            continue
        frequencies = [af]
        groupmax = row.get('groupmax_af', row.get('popmax_af'))
        try:
            groupmax = float(groupmax)
            if math.isfinite(groupmax) and 0 <= groupmax <= 1:
                frequencies.append(groupmax)
        except (TypeError, ValueError):
            pass
        for pop in row.get('populations', []) or []:
            try:
                pop_an, pop_ac = float(pop['an']), float(pop['ac'])
                if math.isfinite(pop_an) and math.isfinite(pop_ac) and pop_an > 0 and 0 <= pop_ac <= pop_an:
                    frequencies.append(pop_ac / pop_an)
            except (KeyError, TypeError, ValueError):
                continue
        usable.append({'af': max(frequencies), 'overall_af': af, 'an': an, 'source': source,
                       'quality_verified': row.get('filters') is not None,
                       'quality_scope': row.get('quality_scope', 'per-allele filters available' if row.get('filters') is not None else 'per-allele filter status unavailable')})
    legacy = []
    for cv in v.get('clinvar_records', []):
        for field in ['AF_EXAC', 'AF_ESP', 'AF_TGP']:
            values = cv.get(field, [])
            for value in values if isinstance(values, (list, tuple)) else [values]:
                try:
                    value = float(value)
                    if math.isfinite(value) and 0 <= value <= 1:
                        legacy.append(value)
                except (TypeError, ValueError):
                    pass
    return {'known': bool(usable), 'af': max((r['af'] for r in usable), default=None),
            'status': 'measured' if usable else 'unknown_or_unusable',
            'source': '; '.join(sorted({r['source'] for r in usable})) if usable else 'missing or unusable; no rarity credit',
            'measurements': usable, 'unreliable_callsets': rejected, 'rarity_credit': 0,
            'all_measurement_quality_verified': bool(usable) and all(r['quality_verified'] for r in usable),
            'legacy_af_for_common_variant_exclusion_only': max(legacy) if legacy else None}


def rcv_evidence(v, disease):
    target = normalize_condition_id(disease['disease_id'])
    rows = v.get('public_annotation', {}).get('rcv_assertions') or []
    rows = rows if isinstance(rows, list) else [rows]
    relevant, ignored = [], Counter()
    for assertion in rows:
        if assertion.get('contributes_to_aggregate') is False:
            ignored['explicit_noncontributing_assertion'] += 1
            continue
        context = assertion.get('classification_type', assertion.get('classification_context'))
        if context and str(context).lower() not in {'germline', 'germlineclassification', 'germline classification'}:
            ignored['non_germline_classification_context'] += 1
            continue
        assertion_type = assertion.get('assertion_type')
        if assertion_type not in {'RCV', 'SCV'}:
            ignored['not_a_condition_specific_assertion'] += 1
            continue
        ids = condition_ids(assertion.get('disease_ids', assertion.get('condition_ids', [])))
        ids |= condition_ids(assertion.get('mapped_disease_ids', []))
        if target not in ids:
            continue
        classification = assertion.get('classification', assertion.get('clinical_significance', ''))
        if isinstance(classification, dict):
            classification = classification.get('description', '')
        support = clinical_support({'CLNSIG': str(classification).replace(' ', '_'),
            'CLNREVSTAT': str(assertion.get('review_status', '')).replace(' ', '_')})
        # Substring matches would incorrectly promote e.g. "not pathogenic" or
        # "association with pathogenic variant" to a clinical classification.
        label = str(classification).lower().replace('_', ' ').strip()
        labels = {part.strip() for part in label.split('/')}
        pathogenic_labels = {'pathogenic', 'likely pathogenic', 'pathogenic, low penetrance', 'likely pathogenic, low penetrance'}
        benign_labels = {'benign', 'likely benign'}
        support['pathogenic_assertion'] = bool(labels) and labels <= pathogenic_labels and not support['conflicting']
        support['benign'] = bool(labels) and labels <= benign_labels and not support['conflicting']
        stars = assertion.get('review_stars', support['review_stars'])
        if not isinstance(stars, int) or not 0 <= stars <= 4:
            raise ValueError('Invalid public RCV review stars')
        support['review_stars'] = stars
        support['support'] = [.65, .75, .85, .95, .95][stars] if support['pathogenic_assertion'] else 0
        support['conflicting'] = support['conflicting'] or bool(assertion.get('conflict'))
        relevant.append({**assertion, **support, 'disease_ids': sorted(ids),
                         'assertion_type': assertion_type, 'classification': classification})
    conflict = any(r['conflicting'] for r in relevant)
    conflict |= any(r['benign'] for r in relevant) and any(r['pathogenic_assertion'] for r in relevant)
    classes = {('pathogenic' if r['pathogenic_assertion'] else 'benign' if r['benign']
                else 'uncertain' if 'uncertain' in str(r['classification']).lower() else 'other') for r in relevant}
    conflict |= len(classes & {'pathogenic', 'benign', 'uncertain'}) > 1
    eligible = [r for r in relevant if r['pathogenic_assertion'] and not conflict]
    return {'support': max((r['support'] for r in eligible), default=0),
            'benign': bool(relevant) and any(r['benign'] for r in relevant) and not conflict and not eligible,
            'conflict': bool(conflict), 'assertions': relevant, 'ignored': dict(ignored),
            'source': 'condition-specific ' + '/'.join(sorted({r['assertion_type'] for r in relevant})) if relevant else 'no condition-specific assertion match'}


def predictor_evidence(v, annotation):
    public = v.get('public_annotation', {})
    rows = public.get('alphamissense_annotations') or ([public['alphamissense']] if public.get('alphamissense') else [])
    rows = rows if isinstance(rows, list) else [rows]
    amino = re.fullmatch(r'(\d+)([A-Z*])>(\d+)([A-Z*])', annotation.get('amino_acid_change', ''))
    expected = f'{amino[2]}{amino[1]}{amino[4]}' if amino and amino[1] == amino[3] else None
    matches = []
    for row in rows:
        if row.get('gene') and row['gene'] != annotation['gene']:
            continue
        if row.get('transcript_id', '').split('.')[0] != annotation['transcript'].split('.')[0]:
            continue
        if not expected or row.get('protein_variant') != expected:
            continue
        score = float(row['am_pathogenicity'])
        if not math.isfinite(score) or not 0 <= score <= 1:
            raise ValueError('Invalid public missense predictor score')
        matches.append({**row, 'am_pathogenicity': score})
    if not matches:
        return {'available': False, 'score': None, 'reason': 'no matching genomic allele, transcript and protein substitution'}
    values = {r['am_pathogenicity'] for r in matches}
    if len(values) > 1:
        return {'available': False, 'score': None, 'reason': 'conflicting duplicate predictor annotations'}
    return {'available': True, 'score': matches[0]['am_pathogenicity'], 'annotation': matches[0],
            'version_match': matches[0].get('transcript_id') == annotation['transcript']}


def mechanism_evidence(annotation, disease, mode, config):
    matching = [r for r in config.get('mechanism_rules', [])
                if r.get('gene') == annotation['gene'] and normalize_condition_id(r.get('disease_id', '')) == normalize_condition_id(disease['disease_id'])
                and r.get('inheritance') == mode]
    valid = [r for r in matching if r.get('source') and r.get('mechanism')]
    if not valid:
        return {'status': 'unreviewed', 'mechanism': None, 'source': None,
                'flags': ['disease_mechanism_requires_review']}
    rule = valid[0]
    effects = set(annotation['consequences'])
    compatible = bool(effects & set(rule.get('supported_consequences', [])))
    flags = list(rule.get('caveat_flags', []))
    if not compatible:
        flags.append('consequence_not_established_for_curated_mechanism')
    if rule.get('transcripts') and annotation['transcript'].split('.')[0] not in {t.split('.')[0] for t in rule['transcripts']}:
        compatible = False
        flags.append('transcript_not_in_curated_mechanism')
    plausible = bool(effects & set(rule.get('plausible_consequences', [])))
    if rule.get('requires_allele_specific_validation'):
        compatible = False
        flags.append('mechanism_cannot_be_assigned_from_consequence_class_alone')
    return {'status': 'supported' if compatible else 'plausible' if plausible else 'uncertain', 'mechanism': rule['mechanism'],
            'source': rule['source'], 'notes': rule.get('notes'), 'flags': flags}


def evidence(v, annotation, disease, config, mode='unknown', no_clinvar=False, no_am=False, predictor_weight=None):
    rcv = rcv_evidence(v, disease) if not no_clinvar else {'support': 0, 'benign': False, 'conflict': False, 'assertions': [], 'source': 'ablated'}
    pred = predictor_evidence(v, annotation) if 'missense' in annotation['consequences'] and not no_am else {'available': False, 'score': None}
    weight = config.get('predictor_weight', 0) if predictor_weight is None else predictor_weight
    if not isinstance(weight, (int, float)) or not 0 <= weight <= 1:
        raise ValueError('Invalid predictor ranking weight')
    effect = annotation['effect_score']
    if pred['available'] and weight:
        effect = (1 - weight) * effect + weight * (.2 + .5 * pred['score'])
    mechanism = mechanism_evidence(annotation, disease, mode, config)
    flags = list(mechanism['flags'])
    genomic_predictors = v.get('public_annotation', {}).get('other_predictors', {})
    splice = genomic_predictors.get('SPLICE_AI', genomic_predictors.get('SPLICEAI', {}))
    splice_score = splice.get('score') if isinstance(splice, dict) else None
    splice_rescue = False
    if splice_score is not None and {'splice_region', 'synonymous'} & set(annotation['consequences']):
        splice_score = float(splice_score)
        if not math.isfinite(splice_score) or not 0 <= splice_score <= 1:
            raise ValueError('Invalid public splice predictor score')
        if splice_score >= config.get('splice_review_threshold', .5):
            effect = max(effect, config.get('splice_review_utility', .4))
            splice_rescue = True
            flags.append('splice_prediction_rescue_requires_transcript_and_RNA_review')
    if genomic_predictors:
        flags.append('genomic_predictor_summary_has_no_transcript_identifiers')
    if not v['population']['known']:
        flags.append('population_frequency_unknown_or_unusable')
    elif not v['population'].get('all_measurement_quality_verified', False):
        flags.append('population_source_preprocessing_used_per_allele_filters_unavailable')
    if not rcv['assertions']:
        flags.append('no_condition_specific_ClinVar_assertion_support')
    if v.get('clinvar_records'):
        flags.append('aggregate_ClinVar_retained_for_review_not_condition_specific_support')
    if rcv['conflict']:
        flags.append('conflicting_condition_specific_ClinVar_evidence')
    if any('low penetrance' in str(a.get('classification', '')).lower() for a in rcv['assertions']):
        flags.append('condition_assertion_low_penetrance_requires_review')
    if 'missense' in annotation['consequences']:
        flags.append('missense_function_unconfirmed')
        if not pred['available']:
            flags.append('missense_predictor_missing_or_incompatible')
        else:
            flags.append('AlphaMissense_computational_evidence_not_functional_assay')
            if not pred['version_match']:
                flags.append('predictor_transcript_version_differs_stable_id_and_substitution_match')
    if set(annotation['consequences']) & LOF:
        flags.append('putative_loss_of_function_requires_NMD_exon_and_mechanism_review')
    return {'score': max(effect, rcv['support']), 'effect_score': annotation['effect_score'],
            'clinical_support': rcv['support'], 'has_pathogenic_assertion_any_condition': bool(rcv['support']),
            'benign': rcv['benign'], 'clinical_conflict': rcv['conflict'],
            'alphamissense_used': bool(pred['available'] and weight), 'predictor_score': pred['score'],
            'predictor_evidence': pred, 'selected_transcript': annotation['transcript'],
            'genomic_predictors_for_review': genomic_predictors, 'splice_prediction_rescue': splice_rescue,
            'rcv_assertions_used': rcv['assertions'], 'condition_assertions_used': rcv['assertions'], 'condition_support_source': rcv['source'],
            'mechanism_status': mechanism['status'], 'mechanism_evidence': mechanism,
            'coverage': {'population': v['population']['status'], 'predictor': pred['available'], 'condition_RCV': any(a['assertion_type']=='RCV' for a in rcv['assertions']),
                         'condition_assertion': bool(rcv['assertions'])},
            'flags': sorted(set(flags))}


def best_transcript_evidence(v, gene, disease, mode, config, **kwargs):
    options = []
    for ann in v['annotations']:
        if ann['gene'] == gene:
            options.append((ann, evidence(v, ann, disease, config, mode, **kwargs)))
    if not options:
        return None
    # Score every transcript first. On equal common-evidence score, matched
    # predictor availability may select a transcript, but does not change score.
    options.sort(key=lambda item: (-item[1]['score'], -item[1]['clinical_support'],
        item[1]['mechanism_status'] != 'supported', not item[1]['predictor_evidence']['available'], item[0]['transcript']))
    ann, ev = options[0]
    if len({tuple(sorted(a['consequences'])) for a, _ in options}) > 1:
        ev['flags'] = sorted(set(ev['flags'] + ['consequences_differ_between_transcripts']))
    ev['transcripts_evaluated'] = [a['transcript'] for a, _ in options]
    return ann, ev


def effective_modes(disease, gene, config):
    modes = set(disease.get('modes', ['unknown']))
    for rule in config.get('mechanism_rules', []):
        if rule.get('source') and rule.get('gene') == gene and normalize_condition_id(rule.get('disease_id', '')) == normalize_condition_id(disease['disease_id']):
            modes.add(rule.get('inheritance', 'unknown'))
    return sorted(modes - {'unknown'}) or ['unknown']


def model(variants, profiles, config, *, pheno_weight=None, mechanism_weight=None, predictor_weight=None,
          no_clinvar=False, no_am=False, unknown_phase_factor=None):
    grouped = defaultdict(list)
    counts = Counter()
    rows = []
    for v in variants:
        if not v['quality_pass']:
            counts['quality_excluded_alleles'] += 1
            continue
        for gene in sorted({a['gene'] for a in v['annotations']}):
            grouped[gene].append(v)
    pw = config['phenotype_weight'] if pheno_weight is None else pheno_weight
    mw = config.get('mechanism_weight', 0) if mechanism_weight is None else mechanism_weight
    pf = config.get('unknown_phase_factor', 1) if unknown_phase_factor is None else unknown_phase_factor
    if not 0 <= pw <= 1 or not 0 <= mw <= 1 or pw + mw > 1 or not 0 <= pf <= 1:
        raise ValueError('Invalid ranking weights')
    for gene, values in sorted(grouped.items()):
        for disease in profiles.get(gene, [{'disease_id': 'unannotated', 'modes': ['unknown'], 'score': 0}]):
            for mode in effective_modes(disease, gene, config):
                candidates = []
                for v in values:
                    choice = best_transcript_evidence(v, gene, disease, mode, config,
                        no_clinvar=no_clinvar, no_am=no_am, predictor_weight=predictor_weight)
                    if choice is None:
                        continue
                    ann, ev = choice
                    if ev['benign']:
                        counts['condition_specific_benign_excluded'] += 1
                        continue
                    if ev['score'] < config.get('minimum_allele_utility', .3):
                        counts['low_utility_excluded'] += 1
                        continue
                    threshold = config['maximum_population_af_ar'] if mode in {'AR', 'unknown', 'XL', 'XLR'} else config['maximum_population_af_ad']
                    pop = v['population']
                    af = pop['af'] if pop['known'] else pop.get('legacy_af_for_common_variant_exclusion_only')
                    if af is not None and af > threshold:
                        counts['observed_common_alleles_excluded'] += 1
                        continue
                    candidates.append((v, ann, ev))
                combinations = []
                if mode == 'AR':
                    combinations += [([x], 'homozygous', 'Observed diploid homozygous genotype') for x in candidates if x[0]['genotype'] == 'hom_alt']
                    heterozygotes = [x for x in candidates if x[0]['genotype'] == 'het']
                    cap = config.get('maximum_pair_alleles_per_gene_disease', 80)
                    if len(heterozygotes) > cap:
                        # Preserve every candidate as a carrier/second-allele review
                        # item. Do not silently discard low-coverage alleles or
                        # fabricate a pair after an arbitrary score truncation.
                        counts['pair_enumeration_deferred_excess_alleles'] += 1
                        combinations += [([x], 'second_allele_unresolved', 'Pair enumeration deferred; all single alleles retained for local review') for x in heterozygotes]
                    elif len(heterozygotes) == 1:
                        combinations.append(([heterozygotes[0]], 'second_allele_unresolved', 'Single recessive allele; carrier or unresolved second-allele finding'))
                    else:
                        for a, b in itertools.combinations(heterozygotes, 2):
                            phase, reason = phase_call(a[0], b[0])
                            if phase == 'cis':
                                counts['cis_pairs_excluded'] += 1
                                continue
                            combinations.append(([a, b], phase, reason))
                else:
                    for x in candidates:
                        if mode in {'XL', 'XLR', 'XLD'} and canonical_chrom(x[0]['chrom']) != 'X':
                            counts['x_linked_non_X_allele_excluded'] += 1
                            continue
                        if mode == 'AD' and x[0]['genotype'] == 'hemizygous':
                            counts['AD_haploid_context_deferred'] += 1
                            continue
                        combinations.append(([x], 'not_applicable', 'Observed single allele; de novo status, sex and penetrance are not inferred'))
                for members, phase, reason in combinations:
                    weak = min(x[2]['score'] for x in members)
                    strong = max(x[2]['score'] for x in members)
                    genetic = .85 * weak + .15 * strong
                    mechanism = float(all(x[2]['mechanism_status'] == 'supported' for x in members))
                    score = ((1 - pw - mw) * genetic + pw * disease['score'] + mw * mechanism) * (pf if phase == 'unknown' else 1)
                    flags = {flag for x in members for flag in x[2]['flags']}
                    requirements = {'local_manual_review', 'confidence_must_be_supplied_separately'}
                    all_rcv = all(x[2]['clinical_support'] > 0 for x in members)
                    if phase == 'unknown':
                        flags.add('phase_unresolved_not_probability')
                        requirements.add('review_phase_and_recessive_two_copy_hypothesis')
                    if mode in {'XL', 'XLR', 'XLD'}:
                        flags.add('X_linked_sex_ploidy_penetrance_and_inheritance_subtype_require_review')
                        requirements.add('review_X_linked_genotype_and_sex_dependent_expression')
                    if mode == 'AD':
                        flags.add('dominant_segregation_de_novo_status_and_penetrance_unresolved')
                        requirements.add('review_dominant_inheritance_and_allele_specific_mechanism')
                    if mode in {'AD', 'XL', 'XLR', 'XLD'} and not all_rcv:
                        flags.add('novel_single_allele_hypothesis_requires_variant_specific_validation')
                    if phase == 'second_allele_unresolved':
                        tier = 'carrier_or_unresolved_second_allele'
                        requirements.add('identify_or_exclude_second_allele_before_primary_recessive_claim')
                    elif mode == 'unknown':
                        tier = 'inheritance_unresolved_review'
                        requirements.add('curate_disease_inheritance')
                    elif all_rcv:
                        tier = 'condition_supported_review'
                    elif mechanism:
                        tier = 'mechanism_supported_review'
                    elif all(x[2]['mechanism_status'] in {'supported', 'plausible'} for x in members):
                        tier = 'mechanism_plausible_review'
                        requirements.add('review_allele_specific_residual_function_or_hypomorph_effect')
                    else:
                        tier = 'variant_mechanism_unresolved_review'
                        requirements.add('review_variant_function_and_disease_mechanism')
                    measured = sum(x[0]['population']['known'] for x in members)
                    if measured != len(members):
                        requirements.add('review_missing_population_callability_and_frequency')
                    rows.append({'gene': gene, 'disease_id': disease['disease_id'], 'inheritance': mode,
                        'members': [x[0]['id'] for x in members], 'phase': phase, 'phase_evidence': reason, 'score': score,
                        'components': {'weakest_allele': weak, 'strongest_allele': strong, 'phenotype': disease['score'],
                                       'literature_mechanism': mechanism, 'phase_factor': pf if phase == 'unknown' else 1},
                        'allele_evidence': [x[2] for x in members], 'flags': sorted(flags),
                        'selected_transcripts': [x[1]['transcript'] for x in members],
                        'coverage': {'population_measured_members': measured, 'member_count': len(members),
                            'population_complete': measured == len(members),
                            'predictor_available_members': sum(x[2]['predictor_evidence']['available'] for x in members),
                            'condition_RCV_supported_members': sum(x[2]['clinical_support'] > 0 for x in members)},
                        'review_tier': tier, 'submission_eligible': False,
                        'eligible_for_manual_selection': phase != 'second_allele_unresolved' and mode != 'unknown',
                        'selection_requirements': sorted(requirements), 'confidence': None})
    # A pathogenic assertion for an unrelated condition must not trump phenotype
    # compatibility merely because its review tier has stronger evidence.
    rows.sort(key=lambda r: (-r['score'], -r['coverage']['condition_RCV_supported_members'],
                            r['gene'], r['members'], r['inheritance'], r['disease_id']))
    # Different diseases/inheritance explanations are retained: same alleles do
    # not establish that their highest HPO-scoring condition is the correct one.
    seen = set()
    unique = []
    for row in rows:
        identity = (row['gene'], tuple(sorted(row['members'])), row['disease_id'], row['inheritance'])
        if identity in seen:
            continue
        seen.add(identity)
        row['rank'] = len(unique) + 1
        unique.append(row)
    counts['ranked_review_hypotheses'] = len(unique)
    counts['submission_eligible_without_review'] = 0
    return unique, dict(counts)


def split_review_hypotheses(rows):
    ranked, review = [], []
    for row in rows:
        target = ranked if row['eligible_for_manual_selection'] else review
        target.append({**row, 'combined_review_rank': row['rank'], 'rank': len(target) + 1})
    return ranked, review


def merge_annotations(*annotations):
    """Merge complementary source channels without erasing prior predictor matches."""
    result = {}
    list_fields = {'population', 'alphamissense_annotations', 'rcv_assertions'}
    for entry in annotations:
        for key, value in entry.items():
            if key in list_fields:
                values = value if isinstance(value, list) else [value]
                result.setdefault(key, []).extend(v for v in values if v)
            elif key == 'other_predictors':
                result.setdefault(key, {}).update(value)
            else:
                result[key] = value
    return result


def attach_condition_aliases(public, aliases):
    """Only use equivalences from a supplied public condition mapping resource."""
    if not isinstance(aliases, dict):
        return public
    aliases = aliases.get('aliases', aliases)
    for assertion in public.get('rcv_assertions', []):
        mapped = set()
        for identifier in condition_ids(assertion.get('disease_ids', [])):
            equivalents = aliases.get(identifier, [])
            mapped |= condition_ids(equivalents)
        if mapped:
            assertion['mapped_disease_ids'] = sorted(mapped)
            assertion['mapping_scope'] = 'public condition equivalence; no variant-level condition transfer'
    return public


def augment_mechanism_profiles(profiles, matcher, query, config):
    """Add separately cited literature profiles without conflating HPO diseases."""
    for rule in config.get('mechanism_rules', []):
        gene = rule['gene']
        if gene not in profiles or not rule.get('source'):
            continue
        if any(normalize_condition_id(p['disease_id']) == normalize_condition_id(rule['disease_id']) for p in profiles[gene]):
            continue
        terms = rule.get('phenotype_hpo_terms', [])
        if not terms:
            continue
        profiles[gene].append({'disease_id': rule['disease_id'], 'modes': [rule['inheritance']],
            'inheritance_source': rule['source'], 'phenotype_source': 'explicit public literature profile',
            **matcher.score_terms(terms, query)})


def join_aggregate_clinvar(variants, path, fasta=None):
    counts = Counter()
    with pysam.VariantFile(str(path)) as reference:
        for rec in reference:
            for alt in rec.alts or []:
                if not set(rec.ref + alt) <= set('ACGT'):
                    continue
                try:
                    # SNVs already have a unique representation; avoid genome-sized FASTA I/O.
                    key = normalize_allele(rec.contig, rec.pos, rec.ref, alt, fasta if len(rec.ref) != len(alt) else None)
                except ValueError:
                    counts['public_ClinVar_reference_mismatch_skipped'] += 1
                    continue
                if key in variants:
                    variants[key]['clinvar_records'].append({'id': rec.id,
                        **{field: list(rec.info[field]) if isinstance(rec.info[field], tuple) else rec.info[field]
                           for field in CV_FIELDS if field in rec.info}})
                    counts['aggregate_ClinVar_matches'] += 1
    return dict(counts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vcf', type=local_path, default=input_path('vcf'))
    parser.add_argument('--annotated-vcf', type=local_path, default=ROOT / 'work/private_analysis.nosync/track1_v3/annotated.bcf')
    parser.add_argument('--fasta', type=local_path, default=ROOT / 'references/track1_v2/hg38.fa')
    parser.add_argument('--gff', type=local_path, required=True)
    parser.add_argument('--phenotype', type=local_path, default=input_path('phenotype'))
    parser.add_argument('--references', type=local_path, default=ROOT / 'references/track1')
    parser.add_argument('--functional', type=local_path, default=ROOT / 'references/track1_v2_functional')
    parser.add_argument('--reference-index', type=local_path)
    parser.add_argument('--exomiser-mvstore', type=local_path)
    parser.add_argument('--mechanism-config', type=local_path)
    parser.add_argument('--config', type=local_path, default=ROOT / 'config/track1_v3.json')
    parser.add_argument('--output', type=local_path, default=ROOT / 'work/private_analysis.nosync/track1_v3')
    parser.add_argument('--reuse-annotation', action='store_true')
    parser.add_argument('--minimal-normalization-only', action='store_true')
    args = parser.parse_args()
    install_network_guard()
    started = time.perf_counter()
    args.output.mkdir(parents=True, exist_ok=True)
    run_started = datetime.now(timezone.utc).isoformat()
    # A failed or interrupted rerun must not leave an earlier successful
    # manifest authorizing a mixture of old and newly written evidence files.
    write_json(args.output / 'run_manifest.json', {
        'provenance_schema_version': 1, 'completed': False,
        'started_utc': run_started,
    })
    args.annotated_vcf.parent.mkdir(parents=True, exist_ok=True)
    config = json.loads(args.config.read_text())
    if args.mechanism_config:
        supplied = json.loads(args.mechanism_config.read_text())
        config['mechanism_rules'] = supplied.get('mechanism_rules', supplied) if isinstance(supplied, dict) else supplied
    if not isinstance(config.get('mechanism_rules', []), list) or not all(isinstance(r, dict) for r in config.get('mechanism_rules', [])):
        raise ValueError('Mechanism rules must be a list of sourced rule objects')
    catalogue_path = args.functional / 'functional_catalogue.json'
    catalogue = json.loads(catalogue_path.read_text()) if catalogue_path.exists() else {}
    annotation_inputs = {'vcf': digest(args.vcf), 'fasta': digest(args.fasta), 'gff': digest(args.gff), 'pysam': pysam.__version__}
    annotation_command = ['-f', str(args.fasta), '-g', str(args.gff), '-l', '-p', 'a', '-B', '30']
    receipt = args.annotated_vcf.with_suffix('.inputs.json')
    with pysam.VariantFile(str(args.vcf)) as source:
        primary = {str(i) for i in range(1, 23)} | {'X', 'Y', 'M'}
        targets = [c for c in source.header.contigs if canonical_chrom(c) in primary]
        excluded = [c for c in source.header.contigs if c not in targets]
    if not targets:
        raise ValueError('No primary human contigs in input')
    if args.reuse_annotation and args.annotated_vcf.exists():
        if not receipt.exists() or json.loads(receipt.read_text()) != annotation_inputs:
            raise ValueError('Annotation cache input provenance mismatch')
    else:
        pysam.bcftools.csq(*annotation_command, '-t', ','.join(targets), '-Ob', '-o', str(args.annotated_vcf), str(args.vcf), catch_stdout=False)
        write_json(receipt, annotation_inputs)
    write_json(args.output / 'contig_scope.json', {'included': targets, 'excluded': excluded})
    variants, annotation_counts = extract_variants(args.annotated_vcf, config)
    fasta = None if args.minimal_normalization_only else pysam.FastaFile(str(args.fasta))
    variants, normalization_counts = normalized_variants(variants, fasta)
    aggregate_counts = join_aggregate_clinvar(variants, args.references / 'clinvar.vcf.gz', fasta)
    if fasta is not None:
        fasta.close()
    index = PublicIndex(args.reference_index)
    mvstore_annotations = {}
    if args.exomiser_mvstore:
        from prepare_track1_v3_references import annotate_keys
        mvstore_annotations = annotate_keys(list(variants), args.exomiser_mvstore, args.output / 'mvstore_annotation')
    for key, variant in variants.items():
        source_key = normalize_allele(variant['source_chrom'], variant['source_pos'], variant['source_ref'], variant['source_alt'])
        legacy = catalogue.get(':'.join(map(str, key)), catalogue.get(':'.join(map(str, source_key)), {}))
        indexed = index.lookup(key)
        mvstore = mvstore_annotations.get(':'.join(map(str, key)), {})
        variant['public_annotation'] = attach_condition_aliases(merge_annotations(legacy, indexed, mvstore), index.metadata.get('condition_aliases', {}))
        variant['public_annotation_sources'] = {'legacy_catalogue_match': bool(legacy), 'reference_index_match': bool(indexed),
                                               'exomiser_mvstore_match': bool(mvstore.get('database_match'))}
        variant['population'] = population_evidence(variant)
    index_metadata = index.metadata
    index.close()
    phenotypes = extract_phenotypes(args.phenotype)
    query = positive_proband_terms(phenotypes)
    matcher = DiseaseMatcher(args.references / 'hp.obo', args.references / 'genes_to_phenotype.txt', config['aliases'])
    genes = {a['gene'] for v in variants.values() for a in v['annotations']}
    # Disease inheritance comes from disease annotations or exact sourced rules;
    # do not invent a gene-wide recessive profile from the broad MVA review panel.
    profiles = {g: matcher.gene(g, query, []) for g in sorted(genes)}
    augment_mechanism_profiles(profiles, matcher, query, config)
    all_rows, review = model(list(variants.values()), profiles, config)
    ranked, review_only = split_review_hypotheses(all_rows)
    scenarios = {}
    for name, options in [('without_phenotype', {'pheno_weight': 0}), ('phenotype_0.2', {'pheno_weight': .2}),
        ('phenotype_0.5', {'pheno_weight': .5}), ('without_condition_RCV', {'no_clinvar': True}),
        ('AlphaMissense_weight_0.5', {'predictor_weight': .5}), ('without_AlphaMissense', {'no_am': True}),
        ('unknown_phase_utility_factor_0.8', {'unknown_phase_factor': .8})]:
        result, counts = model(list(variants.values()), profiles, config, **options)
        scenarios[name] = split_review_hypotheses(result)[0][:50]
    output_hashes = {}
    for name, obj in [('candidate_evidence', list(variants.values())), ('ranked_hypotheses', ranked),
        ('review_hypotheses', review_only), ('disease_phenotype', profiles), ('sensitivity', scenarios), ('review_counts', review),
        ('annotation_counts', {**annotation_counts, **normalization_counts, **aggregate_counts})]:
        output_path = args.output / (name + '.json')
        write_json(output_path, obj)
        output_hashes[output_path.name] = digest(output_path)
    manifest = {'provenance_schema_version': 1, 'completed': True,
        'started_utc': run_started, 'output_hashes': output_hashes,
        'model': config['model'], 'run_utc': datetime.now(timezone.utc).isoformat(), 'config': config,
        'input_hashes': annotation_inputs, 'phenotype_sha256': digest(args.phenotype), 'config_sha256': digest(args.config),
        'code_sha256': digest(Path(__file__)), 'shared_baseline_sha256': digest(ROOT / 'scripts/track1_baseline.py'),
        'shared_v2_sha256': digest(ROOT / 'scripts/track1_v2.py'), 'samtools_bcftools_version': pysam.__samtools_version__,
        'annotation_command': annotation_command + ['-t', ','.join(targets)], 'annotated_vcf_sha256': digest(args.annotated_vcf),
        'annotation_reused': bool(args.reuse_annotation), 'reference_index_sha256': digest(args.reference_index) if args.reference_index else None,
        'reference_index_metadata': index_metadata, 'mechanism_config_sha256': digest(args.mechanism_config) if args.mechanism_config else None,
        'exomiser_mvstore_sha256': digest(args.exomiser_mvstore) if args.exomiser_mvstore else None,
        'exomiser_adapter_code_sha256': digest(ROOT / 'scripts/prepare_track1_v3_references.py') if args.exomiser_mvstore else None,
        'functional_catalogue_sha256': digest(catalogue_path) if catalogue_path.exists() else None,
        'reference_hashes': {name: digest(args.references / name) for name in ['clinvar.vcf.gz', 'hp.obo', 'genes_to_phenotype.txt']},
        'elapsed_seconds': time.perf_counter() - started, 'training': 'none', 'network_during_analysis': False,
        'genome_wide_csq': True, 'annotation_model': 'MANE transcript set; local per-record effects',
        'normalization': 'minimal_trim_only' if args.minimal_normalization_only else 'reference_checked_left_aligned_sequence_indels',
        'score_interpretation': 'uncalibrated review utility; no EPCR or confidence generated',
        'limitations': ['No new clinical classification or diagnostic validation', 'Phase, segregation, sex and de novo status are not inferred',
            'Population annotation coverage and callability may be incomplete; missingness provides no rarity credit',
            'Only exact disease identifiers in actual RCV or SCV assertions provide condition-specific clinical support; aggregate labels are not transferred',
            'Evidence review tiers are labels; combined utility determines rank and no tier guarantees phenotype causality',
            'AlphaMissense is not part of the default common-evidence score; predictor sensitivity is reported separately',
            'MANE and per-record csq do not cover every transcript or multi-variant codon effect',
            'Structural, repeat and noncoding causes remain outside the main model',
            'Large recessive pair searches are explicitly deferred with single-allele review records retained',
            'Every selection and confidence requires separate local review; no submission is authorized by ranking',
            'Public participant exposure preceded development; this is not blinded discovery']}
    write_json(args.output / 'run_manifest.json', manifest)
    print('Offline v3 ranking completed; private outputs saved locally.')


if __name__ == '__main__':
    main()
