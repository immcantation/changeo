"""
Unit tests for DefineClones
"""

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.2.0'
__date__      = '2015.04.29'

# Imports
import time, unittest, os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from IgCore import getScoreDict
from DbCore import IgRecord, readDbFile, getDistMat
import DefineClones as mod
import pandas as pd

# Globals
data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
model_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, 'models')

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
        self.clones = {1:self.records[0:4], 2:self.records[4:]}

        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print "<- %s() %.3f" % (self._testMethodName, t)

    @unittest.skip("-> indexJunctions() skipped\n")
    def test_indexJunctions(self):
        #mod.indexJunctions(db_iter, fields=None, mode='gene', action='first')
        self.fail()

    def test_distanceClones(self):

        # import cProfile
        # prof = cProfile.Profile()
        # results = prof.runcall(mod.distanceClones, self.records, model='hs5f', distance=1.0, dist_mat=self.dist_mat)
        # prof.dump_stats('hs5f-unit-test-dict.prof')

        # hs5f
        results = mod.distanceClones(self.records, model='hs5f', distance=1.0)
        results = {k: sorted(v, key=lambda x: x.id) for k, v in results.iteritems()}
        print 'MODEL> hs5f'
        for k, v in results.iteritems():
            for s in v:
                print '  CLONE-%i> %s' % (k, s.id)

        self.assertEqual(sorted(self.clones.values()), sorted(results.values()))

        # m3n
        results = mod.distanceClones(self.records, model='m3n', distance=1.0)
        results = {k: sorted(v, key=lambda x: x.id) for k, v in results.iteritems()}
        print 'MODEL> m3n'
        for k, v in results.iteritems():
            for s in v:
                print '  CLONE-%i> %s' % (k, s.id)

        self.assertEqual(sorted(self.clones.values()), sorted(results.values()))

        # m1n
        results = mod.distanceClones(self.records, model='m1n', distance=10)
        results = {k: sorted(v, key=lambda x: x.id) for k, v in results.iteritems()}
        print 'MODEL> m1n'
        for k, v in results.iteritems():
            for s in v:
                print '  CLONE-%i> %s' % (k, s.id)
        self.assertEqual(sorted(self.clones.values()), sorted(results.values()))

        # hs1f
        results = mod.distanceClones(self.records, model='hs1f', distance=10)
        results = {k: sorted(v, key=lambda x: x.id) for k, v in results.iteritems()}
        print 'MODEL> hs1f'
        for k, v in results.iteritems():
            for s in v:
                print '  CLONE-%i> %s' % (k, s.id)
        self.assertEqual(sorted(self.clones.values()), sorted(results.values()))

        # aa
        results = mod.distanceClones(self.records, model='aa', distance=3.0)
        results = {k: sorted(v, key=lambda x: x.id) for k, v in results.iteritems()}
        print 'MODEL> aa'
        for k, v in results.iteritems():
            for s in v:
                print '  CLONE-%i> %s' % (k, s.id)

        self.assertEqual(sorted(self.clones.values()), sorted(results.values()))

         # aa
        results = mod.distanceClones(self.records, model='ham', distance=9.0)
        results = {k: sorted(v, key=lambda x: x.id) for k, v in results.iteritems()}
        print 'MODEL> ham'
        for k, v in results.iteritems():
            for s in v:
                print '  CLONE-%i> %s' % (k, s.id)

        self.assertEqual(sorted(self.clones.values()), sorted(results.values()))


if __name__ == '__main__':
    unittest.main()