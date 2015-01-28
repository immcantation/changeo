"""
Unit tests for DefineClones
"""

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.0'
__date__      = '2015.01.28'

# Imports
import time, unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from IgCore import getScoreDict
from DbCore import IgRecord
import DefineClones as mod


class Test_DefineClones(unittest.TestCase):
    def setUp(self):
        print '-> %s()' % self._testMethodName

        # Define common preclone properties
        clone_dict = {'V_CALL':'IGHV6-1*01',
                      'D_CALL':'IGHD6-6*01',
                      'J_CALL':'IGHJ6*02',
                      'JUNCTION_GAP_LENGTH':48}

        # Define unique sequences
        seq_list = [{'SEQUENCE_ID':'A1',
                     'SEQUENCE':'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGTCTGG',
                     'JUNCTION':'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGTCTGG'},
                    {'SEQUENCE_ID':'A2',
                     'SEQUENCE':'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGTCTGG',
                     'JUNCTION':'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGTCTGG'},
                    {'SEQUENCE_ID':'A3',
                     'SEQUENCE':'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGTCTGG',
                     'JUNCTION':'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGTCTGG'},
                    {'SEQUENCE_ID':'A4',
                     'SEQUENCE':'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGNNNNN',
                     'JUNCTION':'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGNNNNN'},
                    {'SEQUENCE_ID':'B1',
                     'SEQUENCE':'TGTGCAAGATATAGCAGCAGCTACTACTACTACGGTATGGACGTCTGG',
                     'JUNCTION':'TGTGCAAGATATAGCAGCAGCTACTACTACTACGGTATGGACGTCTGG'},
                    {'SEQUENCE_ID':'B2',
                     'SEQUENCE':'TGTGNAAGATNTAGCAGCAGCTACTACTACTACGGTATNGACGTCTGG',
                     'JUNCTION':'TGTGNAAGATNTAGCAGCAGCTACTACTACTACGGTATNGACGTCTGG'}]

        # Build preclone IgRecord list
        for x in seq_list:  x.update(clone_dict)
        self.records = [IgRecord(x) for x in seq_list]

        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print "<- %s() %.3f" % (self._testMethodName, t)

    @unittest.skip("-> indexJunctions() skipped\n")
    def test_indexJunctions(self):
        #mod.indexJunctions(db_iter, fields=None, mode='gene', action='first')
        self.fail()

    def test_distanceClones(self):
        score_dict = getScoreDict(n_score=1, gap_score=0, alphabet='aa')
        results = mod.distanceClones(self.records, model='hs5f', distance=3.0, score_dict=score_dict)

        print 'MODEL> hs5f'
        for k, v in results.iteritems():
            for s in v:
                print '  CLONE-%i> %s' % (k, s.id)
        self.fail()


if __name__ == '__main__':
    unittest.main()