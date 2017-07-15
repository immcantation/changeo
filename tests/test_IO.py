"""
Unit tests for MakeDb
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from changeo import __version__, __date__

# Imports
import os
import sys
import time
import unittest
from Bio import SeqIO

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')

# Import script
sys.path.append(os.path.join(test_path, os.pardir, 'bin'))
from changeo.IO import extractIMGT, readRepo, countDbFile, getDbFields


class Test_MakeDb(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)

        # Germline files
        self.ig_repo_dir = '/usr/local/share/germlines/imgt/human/vdj'
        # Read files
        self.ig_imgt_file = os.path.join(data_path, 'imgt_ig.txz')
        # Change-O files
        self.ig_db_file = os.path.join(data_path, 'imgt_ig_db-pass.tsv')

        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print("<- %s() %.3f" % (self._testMethodName, t))

    #@unittest.skip("-> ChangeoReader() skipped\n")
    def test_getDbFields(self):
        # Get fields
        x = getDbFields(self.ig_db_file)
        print(x)

        # Add fields
        x = getDbFields(self.ig_db_file, add=['A', 'B', 'C'])
        print(x)

        # Exclude fields
        x = getDbFields(self.ig_db_file, exclude=['V_SCORE', 'V_IDENTITY', 'J_SCORE', 'J_IDENTITY'])
        print(x)

        # Add and exclude fields
        x = getDbFields(self.ig_db_file, add=['A', 'B', 'C'],
                        exclude=['V_SCORE', 'V_IDENTITY', 'J_SCORE', 'J_IDENTITY'])
        print(x)

        self.fail('TODO')


if __name__ == '__main__':
    unittest.main()