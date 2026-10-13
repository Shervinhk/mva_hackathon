"""Local, conservative phenotype extraction and bounded disease similarity for V4.

No network calls occur here. Document text and derived observations are private
and must be written only to the local analysis directory by the caller.
"""
from __future__ import annotations

from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import re
from zipfile import ZipFile
import xml.etree.ElementTree as ET

from scripts.track1_baseline import PhenotypeMatcher

NS = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
W = '{' + NS['w'] + '}'
HPO = re.compile(r'HP:\d{7}')
FAMILY = re.compile(r'\b(?:family|familial|paternal|maternal|father|mother|parent[s]?|sibling[s]?|sister[s]?|brother[s]?|grandmother|grandfather|aunt|uncle|cousin[s]?|relative[s]?)\b', re.I)
PROBAND = re.compile(r'\b(?:proband|patient|index case|index child|child)\b', re.I)
UNCERTAIN = re.compile(r'\b(?:unknown|uncertain|unclear|possible|possibly|suspected|unconfirmed|not assessed|not evaluated|not known|cannot exclude|cannot be excluded|cannot rule out|cannot be ruled out|not excluded|may have|might have|unavailable|not available|not reported|not documented)\b|\?', re.I)
NEGATIVE = re.compile(r'\b(?:absent|not present|no evidence of|no history of|negative for|negative|denies|without|absence of|not observed|not detected)\b|\bno(?:\s+(?!\d)|$)', re.I)
POSITIVE = re.compile(r'\b(?:present|observed|detected|positive|has|affected)\b', re.I)
MODES = {'HP:0000007': 'AR', 'HP:0000006': 'AD', 'HP:0001417': 'XL',
         'HP:0001419': 'XLR', 'HP:0001423': 'XLD'}


def _path(path):
    if '://' in str(path):
        raise ValueError('Only local files are accepted')
    return Path(path)


def _sha(path):
    with _path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def _text(element):
    """Join text runs without corrupting IDs split by Word formatting."""
    pieces = []
    for node in element.iter():
        if node.tag == W + 't':
            pieces.append(node.text or '')
        elif node.tag in {W + 'tab', W + 'br', W + 'cr'}:
            pieces.append(' ')
        elif node.tag == W + 'p' and pieces:
            pieces.append(' ')
    return ' '.join(''.join(pieces).split())


def _subject(text, default='proband'):
    family, proband = bool(FAMILY.search(text)), bool(PROBAND.search(text))
    if family and proband:
        return 'unknown', ['mixed_family_and_proband_scope']
    if family:
        return 'family', []
    if proband:
        return 'proband', []
    return default, ([] if default != 'unknown' else ['subject_unresolved'])


def _present(text):
    if UNCERTAIN.search(text):
        return None, ['presence_uncertain']
    negatives = list(NEGATIVE.finditer(text))
    # Remove the entire negative phrase before looking for a contradictory
    # positive. For example, 'not present' must not also count as 'present'.
    positive_text = NEGATIVE.sub(' ', text)
    if negatives and POSITIVE.search(positive_text):
        return None, ['mixed_presence_assertions']
    if negatives:
        return False, []
    return True, []


def _context(text, previous):
    """Only short standalone headings change scope for following rows."""
    lowered = text.strip(' :.-').lower()
    if lowered in {'family history', 'family findings', 'family phenotypes',
                   'paternal findings', 'maternal findings', 'sibling findings',
                   'paternal history', 'maternal history', 'siblings',
                   'father', 'mother', 'parents', 'brother', 'sister'}:
        return 'family'
    if lowered in {'proband', 'patient', 'proband findings', 'patient findings',
                   'clinical findings', 'clinical phenotype', 'phenotypes',
                   'proband phenotypes', 'patient phenotypes'}:
        return 'proband'
    return previous


def _row_observations(text, cells, row, table, default):
    terms = list(dict.fromkeys(HPO.findall(text)))
    if not terms:
        return []
    subject, scope_flags = _subject(text, default)
    present, presence_flags = _present(text)
    flags = scope_flags + presence_flags
    if re.search(r'\b(?:example|illustration|template)\b|\be\.g\.', text, re.I):
        present = None
        flags.append('possible_instruction_or_example')
    # Scope and polarity cannot reliably be assigned separately to each term
    # in a mixed multi-term row. Preserve it for local clarification.
    if len(terms) > 1 and (NEGATIVE.search(text) or UNCERTAIN.search(text)):
        present = None
        flags.append('multiple_terms_with_qualifier')
    if subject == 'unknown':
        present = None
    label = next((c for c in cells if c and not HPO.search(c)), '')
    onset = None
    onset_match = re.search(r'\b(?:age at onset|onset)\s*[:=]\s*([^;|]+)', text, re.I)
    if onset_match:
        onset = onset_match.group(1).strip()
    return [{'hpo_id': term, 'label': label, 'subject': subject,
             'present': present, 'onset': onset, 'source_row': row,
             'source_table': table, 'source_text': text,
             'review_required': subject == 'unknown' or present is None,
             'flags': sorted(set(flags))} for term in terms]


def _resolve_contradictions(observations):
    states = defaultdict(set)
    for observation in observations:
        if observation['subject'] == 'proband' and observation['present'] is not None:
            states[observation['hpo_id']].add(observation['present'])
    contradictory = {term for term, values in states.items() if len(values) > 1}
    for observation in observations:
        if observation['subject'] == 'proband' and observation['hpo_id'] in contradictory:
            observation['present'] = None
            observation['review_required'] = True
            observation['flags'] = sorted(set(observation.get('flags', [])) | {'contradictory_proband_observations'})
    return sorted(contradictory)


def extract_observations(path):
    """Extract explicit HPO IDs, retaining negation, family scope and ambiguity.

    An ordinary unqualified clinical table entry is provisionally a positive
    proband observation. This assumption is disclosed in the returned payload;
    extraction never claims a human reviewed the source.
    """
    path = _path(path)
    with ZipFile(path) as archive:
        root = ET.fromstring(archive.read('word/document.xml'))
    body = root.find('w:body', NS)
    if body is None:
        raise ValueError('Clinical document has no Word body')
    observations = []
    context = 'proband'
    row_number = table_number = 0
    for element in body:
        if element.tag == W + 'p':
            text = _text(element)
            context = _context(text, context)
            if HPO.search(text):
                row_number += 1
                # Unscoped prose could be a differential, instruction or family
                # narrative. It requires explicit subject language; clinical
                # table entries retain the disclosed default assumption.
                observations.extend(_row_observations(text, [text], row_number, None, 'unknown'))
        elif element.tag == W + 'tbl':
            table_number += 1
            table_context = context
            for row in element.findall('./w:tr', NS):
                row_number += 1
                cells = [_text(cell) for cell in row.findall('./w:tc', NS)]
                text = ' | '.join(cells)
                if not HPO.search(text):
                    # A table-wide family heading must not turn later rows
                    # into the child's symptoms. Mixed subject headers leave
                    # unqualified rows unresolved rather than guessing.
                    if FAMILY.search(text) and PROBAND.search(text):
                        table_context = 'unknown'
                    elif FAMILY.search(text):
                        table_context = 'family'
                    elif PROBAND.search(text):
                        table_context = 'proband'
                    elif len([cell for cell in cells if cell]) == 1:
                        table_context = _context(text, table_context)
                    continue
                observations.extend(_row_observations(text, cells, row_number, table_number, table_context))
    if not observations:
        raise ValueError('No explicit HPO observations found; provide local structured curation')
    contradictions = _resolve_contradictions(observations)
    return {'schema_version': 1, 'source_sha256': _sha(path),
            'reviewed_locally': False, 'extraction': 'conservative_docx_rules_v4',
            'default_unqualified_clinical_subject': 'proband',
            'default_unqualified_clinical_presence': True,
            'contradictory_terms': contradictions, 'observations': observations}


def _validate_override(value, source_sha256):
    if not isinstance(value, dict) or type(value.get('schema_version')) is not int or value['schema_version'] != 1:
        raise ValueError('Observation override requires schema_version 1')
    if value.get('source_sha256') != source_sha256:
        raise ValueError('Observation override is not bound to this document')
    if type(value.get('reviewed_locally')) is not bool:
        raise ValueError('Observation override requires explicit reviewed_locally boolean')
    observations = value.get('observations')
    if not isinstance(observations, list) or not observations:
        raise ValueError('Observation override requires a nonempty observation list')
    clean = []
    for row in observations:
        if not isinstance(row, dict) or not isinstance(row.get('hpo_id'), str) or not HPO.fullmatch(row['hpo_id']):
            raise ValueError('Invalid HPO observation identifier')
        if row.get('subject') not in {'proband', 'family', 'unknown'}:
            raise ValueError('Invalid observation subject')
        if 'present' not in row or (row['present'] is not None and type(row['present']) is not bool):
            raise ValueError('Observation presence must be true, false, or null')
        if type(row.get('source_row')) is not int or row['source_row'] < 1:
            raise ValueError('Observation requires a positive integer source_row')
        if row.get('onset') is not None and not isinstance(row['onset'], str):
            raise ValueError('Observation onset must be text or null')
        if any(key in row and not isinstance(row[key], str) for key in ('label', 'source_text')):
            raise ValueError('Observation text fields must be strings')
        item = {key: row[key] for key in ('hpo_id', 'subject', 'present', 'source_row')}
        item.update(label=row.get('label', ''), onset=row.get('onset'), source_text=row.get('source_text', ''),
                    review_required=row['subject'] == 'unknown' or row['present'] is None, flags=[])
        if row.get('source_table') is not None:
            if type(row['source_table']) is not int or row['source_table'] < 1:
                raise ValueError('Invalid source_table')
            item['source_table'] = row['source_table']
        clean.append(item)
    result = dict(value)
    result.update(observations=clean, extraction='explicit_local_json_override',
                  contradictory_terms=_resolve_contradictions(clean))
    return result


def load_observations(docx, override=None):
    """Use a hash-bound local override when supplied; otherwise parse the DOCX."""
    if override is None:
        return extract_observations(docx)
    path = _path(override)
    payload = _validate_override(json.loads(path.read_text()), _sha(docx))
    payload['override_sha256'] = _sha(path)
    return payload


def positive_terms(payload, aliases=None):
    """Return only consistent, explicit positive observations of the proband."""
    observations = payload['observations'] if isinstance(payload, dict) else payload
    aliases = aliases or {}
    def canonical(term):
        seen = set()
        while term in aliases:
            if term in seen:
                raise ValueError('Cyclic ontology alternate identifiers')
            seen.add(term)
            term = aliases[term]
        return term
    positive, excluded = set(), set()
    for row in observations:
        if row.get('subject') != 'proband':
            continue
        if row.get('present') is True:
            positive.add(canonical(row['hpo_id']))
        else:
            excluded.add(canonical(row['hpo_id']))
    return sorted(positive - excluded)


class DiseaseMatcherV4:
    """Reuse V3's public ontology/IC resources with corrected scoring semantics."""

    def __init__(self, obo, associations, aliases):
        self.base = PhenotypeMatcher(obo, associations, aliases)
        self.aliases = aliases
        self.profiles = defaultdict(lambda: defaultdict(set))
        header = None
        with _path(associations).open() as stream:
            for line in stream:
                if not line.strip():
                    continue
                values = line.lstrip('#').rstrip('\n').split('\t')
                if {'gene_symbol', 'hpo_id', 'disease_id'} <= set(values):
                    header = values
                    continue
                if line.startswith('#'):
                    continue
                if header is None:
                    raise ValueError('Disease association header not recognized')
                row = dict(zip(header, values))
                gene = aliases.get(row['gene_symbol'], row['gene_symbol'])
                if gene and gene != '-':
                    self.profiles[gene][row['disease_id']].add(self.canonical(row['hpo_id']))

    def canonical(self, term):
        seen = set()
        while term in self.base.alts:
            if term in seen:
                raise ValueError('Cyclic ontology alternate identifiers')
            seen.add(term)
            term = self.base.alts[term]
        return term

    def score_terms(self, terms, query):
        query = {self.canonical(term) for term in query}
        terms = {self.canonical(term) for term in terms}
        zero_ic = {term for term in query if self.base.ic.get(term, 0) <= 0}
        non_phenotype = query & (set(MODES) | {'HP:0000005'})
        unsupported = zero_ic | non_phenotype
        supported = query - unsupported
        # An unsupported descendant must not suppress a supported parent.
        redundant = {term for term in supported if any(term != other and term in self.base.ancestors(other) for other in supported)}
        closure = set().union(*(self.base.ancestors(term) for term in terms)) if terms else set()
        numerator = denominator = 0.0
        details = []
        for term in sorted(query - redundant):
            weight = self.base.ic.get(term, 0)
            common = self.base.ancestors(term) & closure
            match = max(common, key=lambda ancestor: (self.base.ic.get(ancestor, 0), ancestor)) if common else None
            if term in unsupported:
                support = 0.0
                reason = 'non_phenotype_inheritance_term' if term in non_phenotype else 'zero_or_missing_information_content'
            else:
                support = self.base.ic.get(match, 0)
                reason = None
                # Ontology IC must be monotonic; fail explicitly rather than
                # hiding a formula/resource problem by clipping the score.
                if not math.isfinite(weight + support) or support > weight + 1e-12:
                    raise ValueError('Ontology information-content invariant violated')
                numerator += support
                denominator += weight
            details.append({'term': term, 'matched_ancestor': match,
                            'information': weight, 'matched_information': support,
                            'scored': term not in unsupported, 'unscored_reason': reason})
        return {'score': numerator / denominator if denominator else 0.0,
                'terms': details, 'unscored_terms': sorted(unsupported),
                'removed_redundant_terms': sorted(redundant),
                'scored_terms': sorted(supported - redundant),
                'score_interpretation': 'information-weighted phenotype similarity; not a causal probability'}

    def gene(self, gene, query, core):
        gene = self.aliases.get(gene, gene)
        rows = []
        for disease, terms in self.profiles.get(gene, {}).items():
            modes = [mode for term, mode in MODES.items() if term in terms]
            # Generic XL adds no distinct hypothesis when a subtype is given.
            if ('XLR' in modes or 'XLD' in modes) and 'XL' in modes:
                modes.remove('XL')
            rows.append({'disease_id': disease, 'modes': modes or ['unknown'],
                         **self.score_terms(terms, query)})
        if gene in core and not any('AR' in row['modes'] for row in rows):
            rows.append({'disease_id': 'literature_mva_recessive', 'modes': ['AR'],
                         'inheritance_source': 'literature panel; not inferred from genotype',
                         **self.score_terms(self.base.direct.get(gene, set()), query)})
        if not rows:
            rows = [{'disease_id': 'unannotated', 'modes': ['unknown'],
                     **self.score_terms(set(), query)}]
        return sorted(rows, key=lambda row: (-row['score'], row['disease_id']))
