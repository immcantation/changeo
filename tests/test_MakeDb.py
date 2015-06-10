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
        print '-> %s()' % self._testMethodName

        # Define data files
        self.igblast_fmt7_file = os.path.join(data_path, 'igblast_test.fmt7')
        self.igblast_seq_dict = MakeDb.getSeqforIgBlast(os.path.join(data_path, 'igblast_test.fasta'))

        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print "<- %s() %.3f" % (self._testMethodName, t)

    #@unittest.skip("-> readIgBlast() skipped\n")
    def test_readIgBlast(self):
        result = MakeDb.readIgBlast(self.igblast_fmt7_file, self.igblast_seq_dict)
        for x in result:
            print '   ID> %s' % x.id
            print 'VCALL> %s' % x.v_call
            print 'INPUT> %s' % x.seq_input
            print '  VDJ> %s' % x.seq_vdj
            print ' JUNC> %s' % x.junction

        self.fail()


if __name__ == '__main__':
    unittest.main()