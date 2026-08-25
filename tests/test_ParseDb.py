"""
Unit tests for ParseDb
"""

# Info
__author__ = 'Jason Anthony Vander Heiden'
from changeo import __version__, __date__

# Imports
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')
bin_path = os.path.join(test_path, os.pardir, 'bin')

# Import script
sys.path.append(bin_path)
import ParseDb


class Test_ParseDb(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)
        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print("<- %s() %.3f" % (self._testMethodName, t))


class Test_ParseDbGzip(unittest.TestCase):
    """Tests for compressed input and output."""
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        source_file = os.path.join(data_path, 'imgt_ig_db-pass.tsv')
        self.plain_file = os.path.join(self.temp_dir, 'input.tsv')
        self.gzip_file = '%s.gz' % self.plain_file
        shutil.copyfile(source_file, self.plain_file)
        with open(self.plain_file, 'rb') as source, gzip.open(self.gzip_file, 'wb') as target:
            shutil.copyfileobj(source, target)

        with open(self.plain_file, 'r') as handle:
            self.expected_lines = handle.read().splitlines()[:2]
        self.sequence_id = self.expected_lines[1].split('\t')[0]

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def _select(self, input_file, gzip_output=False):
        command = [sys.executable, os.path.join(bin_path, 'ParseDb.py'), 'select',
                   '-d', input_file, '-f', 'SEQUENCE_ID', '-u', self.sequence_id,
                   '--outdir', self.temp_dir]
        if gzip_output:
            command.append('--gzip-output')
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        self.assertEqual(result.returncode, 0, result.stdout.decode('utf-8'))

        input_names = {os.path.basename(self.plain_file), os.path.basename(self.gzip_file)}
        output_names = [name for name in os.listdir(self.temp_dir) if name not in input_names]
        self.assertEqual(len(output_names), 1)
        return os.path.join(self.temp_dir, output_names[0])

    def _assert_selected_gzip_output(self, output_file):
        self.assertTrue(output_file.endswith('.gz'))
        with gzip.open(output_file, 'rt') as handle:
            self.assertEqual(handle.read().splitlines(), self.expected_lines)

    def test_compressed_input_preserves_compressed_output(self):
        output_file = self._select(self.gzip_file)
        self._assert_selected_gzip_output(output_file)

    def test_gzip_output_compresses_plain_input_result(self):
        output_file = self._select(self.plain_file, gzip_output=True)
        self._assert_selected_gzip_output(output_file)


if __name__ == '__main__':
    unittest.main()
