"""
Unit tests for MakeDb
"""

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.0'
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
        self.igblast_fmt7_file = os.path.join(data_path, 'issue_29.fmt7')
        self.igblast_seq_dict = SeqIO.index(os.path.join(data_path, 'issue_29.fasta'),
                                       format='fasta')

        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print "<- %s() %.3f" % (self._testMethodName, t)

    #@unittest.skip("-> readIgBlast() skipped\n")
    def test_readIgBlast(self):
        #print self.igblast_fmt7_file
        #print self.igblast_seq_dict
        result = mod.readIgBlast(self.igblast_fmt7_file, self.igblast_seq_dict)
        for x in result:
            print x

        self.fail()


if __name__ == '__main__':
    unittest.main()