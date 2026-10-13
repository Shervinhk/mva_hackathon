import ast
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import local_inputs


class LocalInputsTests(unittest.TestCase):
    def test_generic_defaults_and_environment_override(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(local_inputs, 'ROOT', Path(directory)):
            with patch.dict(os.environ, {}, clear=True):
                self.assertEqual(local_inputs.input_path('vcf'), Path(directory) / 'data/input.vcf.gz')
            with patch.dict(os.environ, {'MVA_VCF_PATH': 'fixture.vcf.gz'}, clear=True):
                self.assertEqual(local_inputs.input_path('vcf'), Path(directory) / 'fixture.vcf.gz')

    def test_private_configuration_and_environment_precedence(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(local_inputs, 'ROOT', Path(directory)):
            config = Path(directory) / 'work/private_analysis.nosync/local_inputs.json'
            config.parent.mkdir(parents=True)
            config.write_text(json.dumps({'vcf': 'fabricated.vcf.gz'}))
            with patch.dict(os.environ, {}, clear=True):
                self.assertEqual(local_inputs.input_path('vcf'), Path(directory) / 'fabricated.vcf.gz')
            with patch.dict(os.environ, {'MVA_VCF_PATH': '/tmp/override.vcf.gz'}, clear=True):
                self.assertEqual(local_inputs.input_path('vcf'), Path('/tmp/override.vcf.gz'))
            config.write_text(json.dumps({'vcf': ['invalid']}))
            with patch.dict(os.environ, {}, clear=True):
                with self.assertRaisesRegex(ValueError, '^Invalid private input-path configuration.$'):
                    local_inputs.input_path('vcf')

    def test_unknown_input_rejected(self):
        with self.assertRaises(ValueError):
            local_inputs.input_path('unknown')

    def test_source_has_no_case_specific_input_path_defaults(self):
        scripts = Path(__file__).resolve().parents[1] / 'scripts'
        affected = []
        for path in scripts.glob('*.py'):
            if path.name == 'local_inputs.py':
                continue
            for node in ast.walk(ast.parse(path.read_text())):
                if (isinstance(node, ast.Constant) and isinstance(node.value, str)
                        and node.value.startswith('data/')
                        and node.value.endswith(('.vcf.gz', '.docx'))):
                    affected.append(path.name)
        self.assertEqual(affected, [])


if __name__ == '__main__':
    unittest.main()
