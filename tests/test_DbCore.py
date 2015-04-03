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
import DbCore as mod

# Globals
data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')

class Test_DbCore(unittest.TestCase):
    def setUp(self):
        print '-> %s()' % self._testMethodName

        row = {'SEQUENCE_ID':    'TEST1',
               'SEQUENCE_INPUT': 'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGNNNNN',
               'V_CALL':         'IGHV6-1*01,IGHV6-1*02',
               'D_CALL':         'IGHD6-6*01',
               'J_CALL':         'IGHJ6*02',}
        self.igrec = mod.IgRecord(row)

        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print "<- %s() %.3f" % (self._testMethodName, t)

    #@unittest.skip("-> IgRecord() skipped\n")
    def test_IgRecord(self):
        print self.igrec.getAlleleCalls(['v','d','j'])
        print self.igrec.getAlleleCalls(['v','j'])
        print self.igrec.getAlleleCalls(['j','v'])
        print self.igrec.getGeneCalls(['v','d','j'])
        print self.igrec.getGeneCalls(['v','j'])
        print self.igrec.getGeneCalls(['j','v'])
        print self.igrec.getFamilyCalls(['v','d','j'])
        print self.igrec.getFamilyCalls(['v','j'])
        print self.igrec.getFamilyCalls(['j','v'])

        self.fail()


if __name__ == '__main__':
    unittest.main()