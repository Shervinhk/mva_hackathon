"""Participant-supplied entry metadata; never infer account or plan details."""
from __future__ import annotations

import re
from urllib.parse import urlparse

KEYS = ('team', 'github_url', 'ai_disclosure')
PLACEHOLDER = re.compile(r'\b(?:TODO|TBD|REPLACE_ME|FILL_IN)\b|\[(?:TO CONFIRM|REQUIRED)', re.I)


def validate_identity(value):
    if not isinstance(value, dict):
        raise ValueError('Entry identity must be a JSON object.')
    result = {}
    for key in KEYS:
        item = value.get(key, '')
        if not isinstance(item, str):
            raise ValueError('Entry identity values must be text.')
        result[key] = item.strip()
    errors = []
    for key, item in result.items():
        if not item or PLACEHOLDER.search(item):
            errors.append(key)
    if result['github_url']:
        parsed = urlparse(result['github_url'])
        parts = [p for p in parsed.path.split('/') if p]
        if (parsed.scheme != 'https' or parsed.netloc != 'github.com'
                or len(parts) != 2 or parsed.query or parsed.fragment
                or any(not re.fullmatch(r'[A-Za-z0-9_.-]+', p) for p in parts)
                or parts[0].lower() in {'your-org', 'your-username', 'username'}
                or parts[1].lower() in {'your-repo', 'repository'}):
            errors.append('github_url')
    return result, sorted(set(errors))


def participant_lines(value):
    identity, missing = validate_identity(value)
    def text(key):
        item = identity[key]
        if not item or key in missing:
            return '[TO CONFIRM]'
        return item.replace('\r', '').replace('\n', ' ').replace('<', '&lt;').replace('>', '&gt;')
    return ['- Team: ' + text('team'), '- GitHub repository: ' + text('github_url'),
            '- Participant-supplied AI-use disclosure: ' + text('ai_disclosure')]
