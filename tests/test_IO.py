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

# Presto and changeo imports
from changeo.IO import getDbFields, extractIMGT, readGermlines, ChangeoReader, \
                       IgBLASTReader, IgBLASTReaderAA, IHMMuneReader, IMGTReader

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')

# Import script
sys.path.append(os.path.join(test_path, os.pardir, 'bin'))
import MakeDb


class Test_MakeDb(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)

        # Germline files
        self.repo_ig = '/usr/local/share/germlines/imgt/human/vdj'
        self.repo_ig_aa = '/usr/local/share/germlines/imgt/human/vdj_aa'
        # Read files
        self.reads_ig = os.path.join(data_path, 'reads_ig.fasta')
        self.reads_ig_aa = os.path.join(data_path, 'reads_ig_aa.fasta')
        # IMGT output
        self.imgt_ig = os.path.join(data_path, 'imgt_ig.txz')
        # IgBLAST output
        self.igblast_ig = os.path.join(data_path, 'igblast1.7_ig.fmt7')
        self.igblast_ig_aa = os.path.join(data_path, 'igblast1.14_ig_aa.fmt7')
        # iHMMune-Align output
        self.ihmmune_ig = os.path.join(data_path, 'ihmmune_ig.csv')
        # Change-O files
        self.db_ig = os.path.join(data_path, 'imgt_ig_db-pass.tsv')

        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print("<- %s() %.3f" % (self._testMethodName, t))

    @unittest.skip("-> ChangeoReader() skipped\n")
    def test_getDbFields(self):
        # Get fields
        x = getDbFields(self.db_ig)
        print(x)

        # Add fields
        x = getDbFields(self.db_ig, add=['A', 'B', 'C'])
        print(x)

        # Exclude fields
        x = getDbFields(self.db_ig, exclude=['V_SCORE', 'V_IDENTITY', 'J_SCORE', 'J_IDENTITY'])
        print(x)

        # Add and exclude fields
        x = getDbFields(self.db_ig, add=['A', 'B', 'C'],
                        exclude=['V_SCORE', 'V_IDENTITY', 'J_SCORE', 'J_IDENTITY'])
        print(x)

        self.fail('TODO')

    @unittest.skip("-> ChangeoReader() skipped\n")
    def test_ChangeoReader(self):
        # Parse
        with open(self.db_ig, 'r') as f:
            result = ChangeoReader(f)
            for x in result: print(x.toDict())

        self.fail('TODO')

    @unittest.skip("-> IMGTReader() skipped\n")
    def test_IMGTReader(self):
        # Extract IMGT files
        temp_dir, files = extractIMGT(self.imgt_ig)

        # Parse
        with open(files['summary'], 'r') as summary, \
                open(files['gapped'], 'r') as gapped, \
                open(files['ntseq'], 'r') as ntseq, \
                open(files['junction'], 'r') as junction:
            result = IMGTReader(summary, gapped, ntseq, junction, receptor=False)
            for x in result: print(x)

        # Remove IMGT temporary directory
        temp_dir.cleanup()

        self.fail('TODO')

    @unittest.skip("-> IgBLASTReader() skipped\n")
    def test_IgBLASTReader(self):
        # Load germlines and sequences
        seq_dict = MakeDb.getSeqDict(self.reads_ig_aa)
        repo_dict = readGermlines([self.repo_ig])

        # Parse
        with open(self.igblast_ig, 'r') as f:
            result = IgBLASTReader(f, seq_dict, repo_dict, receptor=False)
            for x in result: print(x)

        self.fail('TODO')

    @unittest.skip("-> IgBLASTAAReader() skipped\n")
    def test_IgBLASTAAReader(self):
        # Load germlines and sequences
        seq_dict = MakeDb.getSeqDict(self.reads_ig_aa)
        repo_dict = readGermlines([self.repo_ig_aa])

        # Parse
        with open(self.igblast_ig_aa, 'r') as f:
            result = IgBLASTReaderAA(f, seq_dict, repo_dict, receptor=False)
            for x in result: print(x)

        # self.fail('TODO')

    @unittest.skip("-> IHMMReader() skipped\n")
    def test_IHMMReader(self):
        # Load germlines and sequences
        seq_dict = MakeDb.getSeqDict(self.reads_ig)
        repo_dict = readGermlines([self.repo_ig])

        # Parse
        with open(self.ihmmune_ig, 'r') as f:
            result = IHMMuneReader(f, seq_dict, repo_dict, receptor=False)
            for x in result: print(x)

        self.fail('TODO')

if __name__ == '__main__':
    unittest.main()