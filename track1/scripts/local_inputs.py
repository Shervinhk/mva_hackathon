"""Resolve private input locations without embedding case filenames in source."""
from __future__ import annotations

import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULTS = {'vcf': 'data/input.vcf.gz', 'phenotype': 'data/phenotype.docx'}


def input_path(kind):
    if kind not in DEFAULTS:
        raise ValueError('Unknown local input kind.')
    value = os.environ.get('MVA_' + kind.upper() + '_PATH')
    if not value:
        config = ROOT / 'work/private_analysis.nosync/local_inputs.json'
        if config.is_file():
            try:
                settings = json.loads(config.read_text())
                value = settings.get(kind)
                if value is not None and (not isinstance(value, str) or not value):
                    raise ValueError()
            except (ValueError, TypeError, AttributeError):
                raise ValueError('Invalid private input-path configuration.') from None
    path = Path(value or DEFAULTS[kind]).expanduser()
    return path if path.is_absolute() else ROOT / path
