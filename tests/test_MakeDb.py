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


class Test_MakeDb(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)

        # Define data files
        # Created by: ./igblastn -germline_db_V database/IMGT_Human_IGHV -germline_db_D database/IMGT_Human_IGHD
        # -germline_db_J database/IMGT_Human_IGHJ -domain_system imgt -ig_seqtype Ig -outfmt '7 std qseq'
        # -auxiliary_data optional_file/human_gl.aux -query ~/workspace/changeo/tests/data/igblast_test_ig.fasta
        # -out ~/workspace/changeo/tests/data/igblast_test_ig.fmt7
        self.igblast_ig_fmt7_file = os.path.join(data_path, 'igblast_test_ig.fmt7')
        self.igblast_ig_seq_dict = MakeDb.getSeqforIgBlast(os.path.join(data_path, 'igblast_test_ig.fasta'))

        # Created by: ./igblastn -germline_db_V database/imgt_human_tr_v -germline_db_D database/imgt_human_tr_d
        # -germline_db_J database/imgt_human_tr_j -domain_system imgt -ig_seqtype TCR -outfmt '7 std qseq'
        # -auxiliary_data optional_file/human_gl.aux -query ~/workspace/changeo/tests/data/igblast_test_tr.fasta
        # -out ~/workspace/changeo/tests/data/igblast_test_tr.fmt
        self.igblast_tr_fmt7_file = os.path.join(data_path, 'igblast_test_tr.fmt7')
        self.igblast_tr_seq_dict = MakeDb.getSeqforIgBlast(os.path.join(data_path, 'igblast_test_tr.fasta'))

        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print("<- %s() %.3f" % (self._testMethodName, t))

    #@unittest.skip("-> readIgBlast() skipped\n")
    def test_readIgBlast(self):

        print('Testing IG\n')
        result = MakeDb.readIgBlast(self.igblast_ig_fmt7_file, self.igblast_ig_seq_dict)
        for x in result:
            print('   ID> %s' % x.id)
            print('VCALL> %s' % x.v_call)
            print('INPUT> %s' % x.seq_input)
            print('  VDJ> %s' % x.seq_vdj)
            print(' JUNC> %s\n' % x.junction)

        print('Testing TCR\n')
        result = MakeDb.readIgBlast(self.igblast_tr_fmt7_file, self.igblast_tr_seq_dict)
        for x in result:
            print('   ID> %s' % x.id)
            print('VCALL> %s' % x.v_call)
            print('INPUT> %s' % x.seq_input)
            print('  VDJ> %s' % x.seq_vdj)
            print(' JUNC> %s\n' % x.junction)

        self.fail('TODO')


if __name__ == '__main__':
    unittest.main()