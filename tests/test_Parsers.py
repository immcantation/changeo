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

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')

# Import script
sys.path.append(os.path.join(test_path, os.pardir, 'bin'))
import MakeDb
from changeo.IO import extractIMGT, getRepo
from changeo.Parsers import IgBLASTReader, IMGTReader


class Test_MakeDb(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)

        # Define data files
        # Created by: ./igblastn -germline_db_V database/IMGT_Human_IGHV -germline_db_D database/IMGT_Human_IGHD
        # -germline_db_J database/IMGT_Human_IGHJ -domain_system imgt -ig_seqtype Ig -outfmt '7 std qseq sseq btop'
        # -auxiliary_data optional_file/human_gl.aux -query ~/workspace/changeo/tests/data/igblast_test_ig.fasta
        # -out ~/workspace/changeo/tests/data/igblast_test_ig.fmt7
        self.igblast_ig_fmt7_file = os.path.join(data_path, 'igblast1.4_ig.fmt7')
        self.igblast_ig_seq_dict = MakeDb.getInputSeq(os.path.join(data_path, 'reads_ig.fasta'))

        # Created by: ./igblastn -germline_db_V database/imgt_human_tr_v -germline_db_D database/imgt_human_tr_d
        # -germline_db_J database/imgt_human_tr_j -domain_system imgt -ig_seqtype TCR -outfmt '7 std qseq sseq btop'
        # -auxiliary_data optional_file/human_gl.aux -query ~/workspace/changeo/tests/data/igblast_test_tr.fasta
        # -out ~/workspace/changeo/tests/data/igblast_test_tr.fmt
        self.igblast_tr_fmt7_file = os.path.join(data_path, 'igblast1.4_tr.fmt7')
        self.igblast_tr_seq_dict = MakeDb.getInputSeq(os.path.join(data_path, 'reads_tr.fasta'))

        self.repo_dict = getRepo(['/home/jason/share/germlines/imgt/human/vdj'])

        self.imgt_file = os.path.join(data_path, 'imgt_ig.txz')

        self.start = time.time()

    def tearDown(self):


        t = time.time() - self.start
        print("<- %s() %.3f" % (self._testMethodName, t))

    @unittest.skip("-> IMGTReader() skipped\n")
    def test_IMGTReader(self):
        print('Testing IG\n')
        # Extract IMGT files
        temp_dir, files = extractIMGT(self.imgt_file)
        # Parse IMGT files
        with open(files['summary'], 'r') as summary, \
                open(files['gapped'], 'r') as gapped, \
                open(files['ntseq'], 'r') as ntseq, \
                open(files['junction'], 'r') as junction:
            result = IMGTReader(summary, gapped, ntseq, junction, ig=False)
            for x in result: print(x)
        # Remove IMGT temporary directory
        temp_dir.cleanup()

        self.fail('TODO')

    @unittest.skip("-> IgBlastReader() skipped\n")
    def test_IgBlastReader(self):
        print('Testing IG\n')
        with open(self.igblast_ig_fmt7_file, 'r') as f:
            result = IgBLASTReader(f, self.igblast_ig_seq_dict, self.repo_dict)
            #for x in result: print(list(x))
            for x in result: print(x.toDict())

        with open(self.igblast_ig_fmt7_file, 'r') as f:
            result = IgBLASTReader(f, self.igblast_ig_seq_dict, self.repo_dict, ig=False)
            for x in result: print(list(x))

        self.fail('TODO')

if __name__ == '__main__':
    unittest.main()