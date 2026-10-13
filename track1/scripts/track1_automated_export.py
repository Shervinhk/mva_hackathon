"""Deterministic automated shortlist; ordinal values are not probabilities.

This module performs no file reads or network access. The caller must verify
analysis provenance, original genotypes, and reference alleles before export.
"""
from __future__ import annotations

import math

POLICY = (
    'Automated top-ten distinct complete allele hypotheses in the existing v3 '
    'ranking. The CSV epcr field uses 1/(submission_rank+1) solely to encode '
    'strict relative order for the challenge scorer. These ordinal weights are '
    'not estimated or calibrated probabilities of causality. No case-level human '
    'review, causal confirmation, or phase confirmation is claimed.'
)


def automated_choices(rows, analysis_hash, identifier, limit=10):
    if type(limit) is not int or not 1 <= limit <= 10:
        raise ValueError('Automated shortlist limit must be between one and ten.')
    selected, seen, previous_rank = [], set(), 0
    previous_score = math.inf
    for row in rows:
        rank, score = row.get('rank'), row.get('score')
        if (type(rank) is not int or rank <= previous_rank
                or isinstance(score, bool) or not isinstance(score, (int, float))
                or not math.isfinite(score) or score > previous_score):
            raise ValueError('Automated export requires the unchanged finite model ranking.')
        previous_rank, previous_score = rank, score
        if (row.get('eligible_for_manual_selection') is False
                or row.get('inheritance') in {None, '', 'unknown'}
                or row.get('phase') == 'second_allele_unresolved'):
            continue
        members = row.get('members')
        if (not isinstance(members, list) or not 1 <= len(members) <= 2
                or not all(isinstance(member, str) and member for member in members)
                or len(set(members)) != len(members)):
            raise ValueError('Automated export requires complete one-or-two-allele hypotheses.')
        key = tuple(sorted(members))
        if key in seen:
            continue
        seen.add(key)
        if len(selected) >= limit:
            continue
        number = len(selected) + 1
        choice = {
            'hypothesis_id': identifier(row), 'selected': True,
            'epcr': 1.0 / (number + 1), 'finding_type': 'primary',
            'reason': ('Automated v3 model rank; EPCR is an ordinal ordering weight, '
                       'not a causal probability. No case-level human review.'),
            'evidence_reviewed': False, 'submission_rank': number,
        }
        selected.append((row, choice))
    if not selected:
        raise ValueError('No complete hypotheses are available for automated export.')
    manifest = {
        'schema_version': 1, 'analysis_sha256': analysis_hash,
        'selection_mode': 'automated_ordinal', 'reviewed_locally': False,
        'confidence_policy': POLICY,
        'identity': {'team': '', 'github_url': '', 'ai_disclosure': ''},
        'choices': [choice for _, choice in selected],
    }
    return selected, manifest
