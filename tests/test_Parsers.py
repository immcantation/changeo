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
from changeo.IO import getRepo
from changeo.Parsers import IgBLASTReader, IMGTReader


class Test_MakeDb(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)

        # Define data files
        # Created by: ./igblastn -germline_db_V database/IMGT_Human_IGHV -germline_db_D database/IMGT_Human_IGHD
        # -germline_db_J database/IMGT_Human_IGHJ -domain_system imgt -ig_seqtype Ig -outfmt '7 std qseq sseq btop'
        # -auxiliary_data optional_file/human_gl.aux -query ~/workspace/changeo/tests/data/igblast_test_ig.fasta
        # -out ~/workspace/changeo/tests/data/igblast_test_ig.fmt7
        self.igblast_ig_fmt7_file = os.path.join(data_path, 'igblast_test_ig.fmt7')
        self.igblast_ig_seq_dict = MakeDb.getInputSeq(os.path.join(data_path, 'igblast_test_ig.fasta'))

        # Created by: ./igblastn -germline_db_V database/imgt_human_tr_v -germline_db_D database/imgt_human_tr_d
        # -germline_db_J database/imgt_human_tr_j -domain_system imgt -ig_seqtype TCR -outfmt '7 std qseq sseq btop'
        # -auxiliary_data optional_file/human_gl.aux -query ~/workspace/changeo/tests/data/igblast_test_tr.fasta
        # -out ~/workspace/changeo/tests/data/igblast_test_tr.fmt
        self.igblast_tr_fmt7_file = os.path.join(data_path, 'igblast_test_tr.fmt7')
        self.igblast_tr_seq_dict = MakeDb.getInputSeq(os.path.join(data_path, 'igblast_test_tr.fasta'))

        self.repo_dict = getRepo(['/home/jason/share/germlines/imgt/human/vdj'])

        self.imgt_files = [os.path.join(data_path, 'run_21_TR_all', '1_Summary.txt'),
                           os.path.join(data_path, 'run_21_TR_all', '2_IMGT-gapped-nt-sequences.txt'),
                           os.path.join(data_path, 'run_21_TR_all', '3_Nt-sequences.txt'),
                           os.path.join(data_path, 'run_21_TR_all', '6_Junction.txt')]

        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print("<- %s() %.3f" % (self._testMethodName, t))

    @unittest.skip("-> IMGTReader() skipped\n")
    def test_IMGTReader(self):
        print('Testing IG\n')
        handles = [open(f, 'r') for f in self.imgt_files]
        result = IMGTReader(*handles, ig=False)
        for x in result: print(x)

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