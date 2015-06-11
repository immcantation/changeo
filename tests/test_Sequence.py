"""
Unit tests for MakeDb
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from changeo import __version__, __date__

# Imports
import os
import time
import unittest

# Presto and changeo imports
from changeo.Receptor import IgRecord

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')


class Test_DbCore(unittest.TestCase):
    def setUp(self):
        print '-> %s()' % self._testMethodName

        row = {'SEQUENCE_ID':    'TEST1',
               'SEQUENCE_INPUT': 'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGNNNNN',
               'V_CALL':         'IGHV6-1*01,IGHV6-1*02',
               'D_CALL':         'IGHD6-6*01',
               'J_CALL':         'IGHJ6*02',}
        self.ig_rec = IgRecord(row)

        row = {'SEQUENCE_ID':    'TEST2',
               'SEQUENCE_INPUT': 'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGNNNNN',
               'V_CALL':         'TRAV6-1*01,TRAV6-1*02',
               'D_CALL':         'TRAD6-6*01',
               'J_CALL':         'TRAJ6*02',}
        self.tr_rec = IgRecord(row)

        row = {'SEQUENCE_ID':    'TEST3',
               'SEQUENCE_INPUT': '',
               'V_CALL':         '',
               'D_CALL':         None,
               'J_CALL':         'princess,banana,hammock',}
        self.bad_rec = IgRecord(row)

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