"""Synthetic orchestration checks; no challenge files or results are opened."""
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import run_track1_v3_local as runner


class PrivateOutputRoutingTests(unittest.TestCase):
    def test_all_pipeline_steps_use_the_same_excluded_output_root(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            index, store = root / 'public.sqlite', root / 'public.mv.db'
            index.write_text('synthetic'); store.write_text('synthetic')
            status_path = root / 'work/private_analysis.nosync/track1_v3/execution_status.json'
            status_path.parent.mkdir(parents=True)
            status_path.write_text(json.dumps({'completed': True, 'steps': [{'step': 'old_run'}]}))
            calls = []
            def fake_run(command, **kwargs):
                if not calls:
                    status = json.loads(status_path.read_text())
                    self.assertFalse(status['completed'])
                    self.assertEqual(status['steps'], [])
                calls.append(command)
                return SimpleNamespace(returncode=0)
            argv = ['run', '--reference-index', str(index), '--exomiser-mvstore', str(store)]
            with mock.patch.object(runner, 'ROOT', root), mock.patch.object(sys, 'argv', argv), \
                    mock.patch.object(runner.subprocess, 'run', side_effect=fake_run), \
                    contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(runner.main(), 0)
            private = root.resolve() / 'work/private_analysis.nosync'
            self.assertEqual(len(calls), 3)
            for command, name in zip(calls, ['track1_v3', 'track1_exomiser', 'submission_track1_v3']):
                self.assertEqual(command[command.index('--output') + 1], str(private / name))
            self.assertEqual(calls[2][calls[2].index('--analysis') + 1], str(private / 'track1_v3'))
            self.assertEqual(calls[2][calls[2].index('--comparator-output') + 1], str(private / 'track1_exomiser'))
            self.assertIn('--automated', calls[2])
            self.assertTrue(json.loads((private / 'track1_v3/execution_status.json').read_text())['completed'])
            self.assertFalse((root / 'outputs').exists())


if __name__ == '__main__':
    unittest.main()
