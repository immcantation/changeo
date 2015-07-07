"""
Unit tests for MakeDb
"""

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.2.1'
__date__      = '2015.04.03'

# Imports
import os, time, unittest
from Bio import SeqIO
from DbCore import IgRecord
import MakeDb as mod

# Globals
data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')

class Test_MakeDb(unittest.TestCase):
    def setUp(self):
        print '-> %s()' % self._testMethodName

        # Define data files
        # Created by: ./igblastn -germline_db_V database/IMGT_Human_IGHV -germline_db_D database/IMGT_Human_IGHD
        # -germline_db_J database/IMGT_Human_IGHJ -domain_system imgt -ig_seqtype Ig -outfmt '7 std qseq'
        # -auxiliary_data optional_file/human_gl.aux -query ~/workspace/changeo/tests/data/igblast_test_ig.fasta
        # -out ~/workspace/changeo/tests/data/igblast_test_ig.fmt7
        self.igblast_ig_fmt7_file = os.path.join(data_path, 'igblast_test_ig.fmt7')
        self.igblast_ig_seq_dict = mod.getSeqforIgBlast(os.path.join(data_path, 'igblast_test_ig.fasta'))

        # Created by: ./igblastn -germline_db_V database/imgt_human_tr_v -germline_db_D database/imgt_human_tr_d
        # -germline_db_J database/imgt_human_tr_j -domain_system imgt -ig_seqtype TCR -outfmt '7 std qseq'
        # -auxiliary_data optional_file/human_gl.aux -query ~/workspace/changeo/tests/data/igblast_test_tr.fasta
        # -out ~/workspace/changeo/tests/data/igblast_test_tr.fmt
        self.igblast_tr_fmt7_file = os.path.join(data_path, 'igblast_test_tr.fmt7')
        self.igblast_tr_seq_dict = mod.getSeqforIgBlast(os.path.join(data_path, 'igblast_test_tr.fasta'))

        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print "<- %s() %.3f" % (self._testMethodName, t)

    #@unittest.skip("-> readIgBlast() skipped\n")
    def test_readIgBlast(self):

        print 'Testing IG\n'
        result = mod.readIgBlast(self.igblast_ig_fmt7_file, self.igblast_ig_seq_dict)
        for x in result:
            print '   ID> %s' % x.id
            print 'VCALL> %s' % x.v_call
            print 'INPUT> %s' % x.seq_input
            print '  VDJ> %s' % x.seq_vdj
            print ' JUNC> %s\n' % x.junction

        print 'Testing TCR\n'
        result = mod.readIgBlast(self.igblast_tr_fmt7_file, self.igblast_tr_seq_dict)
        for x in result:
            print '   ID> %s' % x.id
            print 'VCALL> %s' % x.v_call
            print 'INPUT> %s' % x.seq_input
            print '  VDJ> %s' % x.seq_vdj
            print ' JUNC> %s\n' % x.junction


        self.fail()


if __name__ == '__main__':
    unittest.main()