"""
Unit tests for ParseDb
"""

# Info
__author__ = 'Jason Anthony Vander Heiden'
from changeo import __version__, __date__

# Imports
import csv
import os
import shutil
import sys
import tempfile
import time
import unittest

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')

# Import script
sys.path.append(os.path.join(test_path, os.pardir, 'bin'))
import ParseDb


class Test_ParseDb(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)
        self.start = time.time()
        self.tmp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)
        t = time.time() - self.start
        print("<- %s() %.3f" % (self._testMethodName, t))

    def _write_tsv(self, name, rows):
        path = os.path.join(self.tmp_dir, name)
        with open(path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        return path

    def _read_tsv(self, path):
        with open(path, 'r', newline='') as f:
            return list(csv.DictReader(f, delimiter='\t'))

    # @unittest.skip("-> mergeDbFiles() skipped\n")
    def test_mergeDbFiles(self):
        file1 = self._write_tsv('sample1.tsv', [{'SEQUENCE_ID': 'A1', 'V_CALL': 'IGHV1'},
                                                 {'SEQUENCE_ID': 'A2', 'V_CALL': 'IGHV2'}])
        file2 = self._write_tsv('sample2.tsv', [{'SEQUENCE_ID': 'B1', 'V_CALL': 'IGHV3'}])
        out_file = os.path.join(self.tmp_dir, 'merged.tsv')

        # Merge without sample annotation
        ParseDb.mergeDbFiles([file1, file2], out_file=out_file)
        records = self._read_tsv(out_file)
        self.assertEqual(len(records), 3)
        self.assertNotIn('SAMPLE', records[0])

        # Merge with sample annotation
        ParseDb.mergeDbFiles([file1, file2], field='SAMPLE', values=['sample1', 'sample2'],
                              out_file=out_file)
        records = self._read_tsv(out_file)
        results = {r['SEQUENCE_ID']: r['SAMPLE'] for r in records}
        self.assertDictEqual(results, {'A1': 'sample1', 'A2': 'sample1', 'B1': 'sample2'})


if __name__ == '__main__':
    unittest.main()