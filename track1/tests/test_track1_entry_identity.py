import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from track1_entry_identity import validate_identity, participant_lines


class EntryIdentityTests(unittest.TestCase):
    def test_actual_participant_fields_are_retained(self):
        identity = {'team': 'Example Team', 'github_url': 'https://github.com/example/mva-code',
                    'ai_disclosure': 'Provider and tier recorded by this synthetic participant; setting not independently verified.'}
        result, errors = validate_identity(identity)
        self.assertEqual(errors, [])
        self.assertEqual(result, identity)
        self.assertIn('Example Team', '\n'.join(participant_lines(identity)))

    def test_unknown_fields_remain_explicit(self):
        self.assertEqual(validate_identity({})[1], ['ai_disclosure', 'github_url', 'team'])
        self.assertIn('[TO CONFIRM]', '\n'.join(participant_lines({})))
        for value in ('https://github.com/your-org/your-repo', 'https://github.com/example/repo/issues',
                      'https://github.com.evil.example/a/b', 'https://github.com/a/b?x=1'):
            self.assertIn('github_url', validate_identity({'github_url': value})[1])

    def test_no_false_attestations_or_html_insertion(self):
        result = '\n'.join(participant_lines({'team': '<script>\nexample',
                            'ai_disclosure': 'Plan unknown; no independent verification.'}))
        self.assertNotIn('<script>', result)
        self.assertIn('Plan unknown; no independent verification.', result)
        self.assertIn('ai_disclosure', validate_identity({'ai_disclosure': '[TO CONFIRM]'} )[1])


if __name__ == '__main__':
    unittest.main()
