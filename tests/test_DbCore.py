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
        self.ig_rec = mod.IgRecord(row)

        row = {'SEQUENCE_ID':    'TEST2',
               'SEQUENCE_INPUT': 'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGNNNNN',
               'V_CALL':         'TRAV6-1*01,TRAV6-1*02',
               'D_CALL':         'TRAD6-6*01',
               'J_CALL':         'TRAJ6*02',}
        self.tr_rec = mod.IgRecord(row)

        row = {'SEQUENCE_ID':    'TEST3',
               'SEQUENCE_INPUT': '',
               'V_CALL':         '',
               'D_CALL':         None,
               'J_CALL':         'princess,banana,hammock',}
        self.bad_rec = mod.IgRecord(row)

        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print '<- %s() %.3f' % (self._testMethodName, t)

    #@unittest.skip("-> IgRecord() skipped\n")
    def test_IgRecord(self):
        print 'IG>'
        print self.ig_rec.getAlleleCalls(['v','d','j'], action='first')
        print self.ig_rec.getAlleleCalls(['v','j'], action='first')
        print self.ig_rec.getAlleleCalls(['j','v'], action='first')
        print self.ig_rec.getGeneCalls(['v','d','j'], action='first')
        print self.ig_rec.getGeneCalls(['v','j'], action='first')
        print self.ig_rec.getGeneCalls(['j','v'], action='first')
        print self.ig_rec.getFamilyCalls(['v','d','j'], action='first')
        print self.ig_rec.getFamilyCalls(['v','j'], action='first')
        print self.ig_rec.getFamilyCalls(['j','v'], action='first')

        print self.ig_rec.getAlleleCalls(['v','d','j'], action='set')
        print self.ig_rec.getAlleleCalls(['v','j'], action='set')
        print self.ig_rec.getAlleleCalls(['j','v'], action='set')
        print self.ig_rec.getGeneCalls(['v','d','j'], action='set')
        print self.ig_rec.getGeneCalls(['v','j'], action='set')
        print self.ig_rec.getGeneCalls(['j','v'], action='set')
        print self.ig_rec.getFamilyCalls(['v','d','j'], action='set')
        print self.ig_rec.getFamilyCalls(['v','j'], action='set')
        print self.ig_rec.getFamilyCalls(['j','v'], action='set')

        print self.ig_rec.getAlleleCalls(['v','d','j'], action='list')
        print self.ig_rec.getAlleleCalls(['v','j'], action='list')
        print self.ig_rec.getAlleleCalls(['j','v'], action='list')
        print self.ig_rec.getGeneCalls(['v','d','j'], action='list')
        print self.ig_rec.getGeneCalls(['v','j'], action='list')
        print self.ig_rec.getGeneCalls(['j','v'], action='list')
        print self.ig_rec.getFamilyCalls(['v','d','j'], action='list')
        print self.ig_rec.getFamilyCalls(['v','j'], action='list')
        print self.ig_rec.getFamilyCalls(['j','v'], action='list')

        print 'TR>'
        print self.tr_rec.getAlleleCalls(['v','d','j'], action='first')
        print self.tr_rec.getAlleleCalls(['v','j'], action='first')
        print self.tr_rec.getAlleleCalls(['j','v'], action='first')
        print self.tr_rec.getGeneCalls(['v','d','j'], action='first')
        print self.tr_rec.getGeneCalls(['v','j'], action='first')
        print self.tr_rec.getGeneCalls(['j','v'], action='first')
        print self.tr_rec.getFamilyCalls(['v','d','j'], action='first')
        print self.tr_rec.getFamilyCalls(['v','j'], action='first')
        print self.tr_rec.getFamilyCalls(['j','v'], action='first')

        print self.tr_rec.getAlleleCalls(['v','d','j'], action='set')
        print self.tr_rec.getAlleleCalls(['v','j'], action='set')
        print self.tr_rec.getAlleleCalls(['j','v'], action='set')
        print self.tr_rec.getGeneCalls(['v','d','j'], action='set')
        print self.tr_rec.getGeneCalls(['v','j'], action='set')
        print self.tr_rec.getGeneCalls(['j','v'], action='set')
        print self.tr_rec.getFamilyCalls(['v','d','j'], action='set')
        print self.tr_rec.getFamilyCalls(['v','j'], action='set')
        print self.tr_rec.getFamilyCalls(['j','v'], action='set')

        print self.tr_rec.getAlleleCalls(['v','d','j'], action='list')
        print self.tr_rec.getAlleleCalls(['v','j'], action='list')
        print self.tr_rec.getAlleleCalls(['j','v'], action='list')
        print self.tr_rec.getGeneCalls(['v','d','j'], action='list')
        print self.tr_rec.getGeneCalls(['v','j'], action='list')
        print self.tr_rec.getGeneCalls(['j','v'], action='list')
        print self.tr_rec.getFamilyCalls(['v','d','j'], action='list')
        print self.tr_rec.getFamilyCalls(['v','j'], action='list')
        print self.tr_rec.getFamilyCalls(['j','v'], action='list')

        print 'JUNK>'
        print self.bad_rec.getAlleleCalls(['v','d','j'], action='first')
        print self.bad_rec.getAlleleCalls(['v','j'], action='first')
        print self.bad_rec.getAlleleCalls(['j','v'], action='first')
        print self.bad_rec.getGeneCalls(['v','d','j'], action='first')
        print self.bad_rec.getGeneCalls(['v','j'], action='first')
        print self.bad_rec.getGeneCalls(['j','v'], action='first')
        print self.bad_rec.getFamilyCalls(['v','d','j'], action='first')
        print self.bad_rec.getFamilyCalls(['v','j'], action='first')
        print self.bad_rec.getFamilyCalls(['j','v'], action='first')

        print self.bad_rec.getAlleleCalls(['v','d','j'], action='set')
        print self.bad_rec.getAlleleCalls(['v','j'], action='set')
        print self.bad_rec.getAlleleCalls(['j','v'], action='set')
        print self.bad_rec.getGeneCalls(['v','d','j'], action='set')
        print self.bad_rec.getGeneCalls(['v','j'], action='set')
        print self.bad_rec.getGeneCalls(['j','v'], action='set')
        print self.bad_rec.getFamilyCalls(['v','d','j'], action='set')
        print self.bad_rec.getFamilyCalls(['v','j'], action='set')
        print self.bad_rec.getFamilyCalls(['j','v'], action='set')

        print self.bad_rec.getAlleleCalls(['v','d','j'], action='list')
        print self.bad_rec.getAlleleCalls(['v','j'], action='list')
        print self.bad_rec.getAlleleCalls(['j','v'], action='list')
        print self.bad_rec.getGeneCalls(['v','d','j'], action='list')
        print self.bad_rec.getGeneCalls(['v','j'], action='list')
        print self.bad_rec.getGeneCalls(['j','v'], action='list')
        print self.bad_rec.getFamilyCalls(['v','d','j'], action='list')
        print self.bad_rec.getFamilyCalls(['v','j'], action='list')
        print self.bad_rec.getFamilyCalls(['j','v'], action='list')

        self.fail()


if __name__ == '__main__':
    unittest.main()