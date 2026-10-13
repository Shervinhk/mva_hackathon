#!/usr/bin/env python3
"""Validate a Track 1 package locally, without scoring or making a submission.

The official scorer is permissive in some places. This validator additionally
rejects incomplete pairs, duplicate hypotheses, non-sequence alleles and malformed
CSV records before upload. These are local quality checks, not new official rules.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import re
from urllib.parse import urlparse


FIELDS = (
    'proband_id', 'chrom_1', 'pos_1', 'ref_1', 'alt_1',
    'chrom_2', 'pos_2', 'ref_2', 'alt_2', 'epcr', 'finding_type', 'notes',
)
CHROMS = {f'chr{i}' for i in range(1, 23)} | {'chrX', 'chrY', 'chrM'}
DNA = re.compile(r'^[ACGTN]+$')
PLACEHOLDER = re.compile(
    r'\b(?:TODO|TBD|REPLACE_ME|FILL_IN)\b|\[(?:TO CONFIRM|REQUIRED)[^\]]*\]',
    re.IGNORECASE,
)


def validate_csv(path: Path | str) -> dict:
    """Return errors/warnings; do not expose candidate alleles in the report."""
    path = Path(path)
    errors, warnings, rows = [], [], []
    try:
        with path.open(newline='', encoding='utf-8-sig') as handle:
            reader = csv.DictReader(handle, strict=True)
            header = reader.fieldnames
            if header is None:
                return {'valid': False, 'rows': 0, 'errors': ['CSV is empty.'], 'warnings': []}
            if len(header) != len(set(header)):
                errors.append('Duplicate column names are not allowed.')
            missing = sorted(set(FIELDS) - set(header))
            extra = sorted(set(header) - set(FIELDS))
            if missing:
                errors.append('Missing template columns: ' + ', '.join(missing))
            if extra:
                errors.append('Unexpected columns: ' + ', '.join(extra))
            if not missing and not extra and tuple(header) != FIELDS:
                warnings.append('Column order differs from the official template; names are valid.')
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as exc:
        return {'valid': False, 'rows': 0, 'errors': [f'Cannot read CSV: {exc}'], 'warnings': []}

    if not 1 <= len(rows) <= 10:
        errors.append(f'Expected 1–10 candidate rows; found {len(rows)}.')
    hypotheses, epcrs = set(), []
    has_secondary = False
    for number, raw in enumerate(rows, 1):
        prefix = f'Row {number}: '
        if None in raw or any(raw.get(field) is None for field in FIELDS):
            errors.append(prefix + 'Field count does not match the template header.')
            continue
        row = {key: value.strip() for key, value in raw.items()}
        if row['proband_id'] != 'PROBAND01':
            errors.append(prefix + 'proband_id must be PROBAND01.')
        variants = []
        for suffix in ('1', '2'):
            keys = [f'{name}_{suffix}' for name in ('chrom', 'pos', 'ref', 'alt')]
            values = [row[key] for key in keys]
            if suffix == '2' and not any(values):
                continue
            if not all(values):
                errors.append(prefix + f'Variant {suffix} requires all four coordinate/allele fields.')
                continue
            chrom, pos, ref, alt = values
            valid = True
            if chrom not in CHROMS:
                errors.append(prefix + f'Variant {suffix} chromosome must be chr1–chr22, chrX, chrY or chrM.')
                valid = False
            if not re.fullmatch(r'[1-9][0-9]*', pos):
                errors.append(prefix + f'Variant {suffix} position must be a positive integer.')
                valid = False
            if not DNA.fullmatch(ref.upper()) or not DNA.fullmatch(alt.upper()):
                errors.append(prefix + f'Variant {suffix} requires explicit DNA alleles (ACGTN), not symbolic or multiple ALT values.')
                valid = False
            elif ref.upper() == alt.upper():
                errors.append(prefix + f'Variant {suffix} REF and ALT are identical.')
                valid = False
            elif 'N' in ref.upper() + alt.upper():
                warnings.append(prefix + f'Variant {suffix} contains N; verify the exact submission allele.')
            if valid:
                variants.append((chrom, int(pos), ref.upper(), alt.upper()))
        if len(variants) == 2 and variants[0] == variants[1]:
            errors.append(prefix + 'The two pair members are identical.')
        hypothesis = frozenset(variants)
        if hypothesis:
            if hypothesis in hypotheses:
                errors.append(prefix + 'Duplicate hypothesis, including reversed pair order.')
            hypotheses.add(hypothesis)
        try:
            epcr = float(row['epcr'])
            if not math.isfinite(epcr) or not 0 < epcr <= 1:
                raise ValueError
            epcrs.append(epcr)
        except ValueError:
            errors.append(prefix + 'epcr must be a finite number in (0, 1].')
        if row['finding_type'] not in {'primary', 'secondary'}:
            errors.append(prefix + 'finding_type must be primary or secondary.')
        has_secondary |= row['finding_type'] == 'secondary'
        if PLACEHOLDER.search(row['notes']):
            errors.append(prefix + 'Notes contain an unresolved placeholder.')
    if epcrs != sorted(epcrs, reverse=True):
        warnings.append('Rows are not in descending EPCR order; the official scorer will sort them.')
    if len(set(epcrs)) != len(epcrs):
        warnings.append('EPCR ties preserve row order for rank and include all tied rows at an F-max threshold.')
    if has_secondary:
        warnings.append('The current public scorer includes secondary rows despite the prose claiming exclusion; review their evidence and confidence.')
    return {'valid': not errors, 'rows': len(rows), 'errors': errors, 'warnings': warnings}


def validate_package(csv_path: Path | str, methods: Path | str | None = None,
                     github: str | None = None, require_package: bool = False) -> dict:
    """File checks cannot establish authentication, clinical accuracy or consent."""
    result = validate_csv(csv_path)
    errors, warnings = list(result['errors']), list(result['warnings'])
    methods_present, github_supplied = False, False
    if methods is not None:
        methods = Path(methods)
        if not methods.is_file() or methods.stat().st_size == 0:
            errors.append('Methods report is missing or empty.')
        elif methods.suffix.lower() not in {'.md', '.pdf'}:
            errors.append('Methods report must use .md or .pdf.')
        elif methods.suffix.lower() == '.md':
            try:
                report = methods.read_text(encoding='utf-8')
                if not report.strip():
                    errors.append('Methods report has no text.')
                elif PLACEHOLDER.search(report):
                    errors.append('Methods report contains unresolved TODO/TBD/confirmation placeholders.')
                else:
                    methods_present = True
            except (OSError, UnicodeError) as exc:
                errors.append(f'Cannot read methods report: {exc}')
        else:
            try:
                with methods.open('rb') as handle:
                    if handle.read(5) != b'%PDF-':
                        errors.append('Methods file is not a PDF despite its extension.')
                    else:
                        methods_present = True
                        warnings.append('PDF content, disclosure and placeholders need human review.')
            except OSError as exc:
                errors.append(f'Cannot read methods report: {exc}')
    elif require_package:
        errors.append('--methods is required for a complete local package.')
    if github:
        url = urlparse(github)
        parts = [part for part in url.path.split('/') if part]
        if (url.scheme != 'https' or url.netloc != 'github.com' or len(parts) != 2
                or url.query or url.fragment or any(PLACEHOLDER.search(p) for p in parts)
                or parts[0].lower() in {'your-org', 'your-username', 'username'}
                or parts[1].lower() in {'your-repo', 'repository'}):
            errors.append('Supply a real https://github.com/OWNER/REPOSITORY URL without query or fragment.')
        else:
            github_supplied = True
    elif require_package:
        errors.append('--github is required for a complete local package.')
    return {
        'local_validation_passed': not errors,
        'csv_valid': result['valid'], 'rows': result['rows'],
        'package_files_present': methods_present and github_supplied,
        'remote_submission_status': 'not_submitted_by_this_tool',
        'errors': errors, 'warnings': warnings,
        'not_verified': [
            'GRCh38 reference allele correctness, clinical validity and phase',
            'EPCR calibration and scientific evidence',
            'AI provider, model, actual plan/tier and data-handling disclosure accuracy',
            'GitHub repository existence, accessibility and reproducibility',
            'HF authentication, remaining quota and rule acceptance',
            'Service acceptance or challenge score',
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--csv', type=Path, required=True)
    parser.add_argument('--methods', type=Path)
    parser.add_argument('--github')
    parser.add_argument('--require-package', action='store_true',
                        help='Require report and GitHub URL in addition to the CSV.')
    args = parser.parse_args()
    result = validate_package(args.csv, args.methods, args.github, args.require_package)
    print(json.dumps(result, indent=2))
    return 0 if result['local_validation_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
