"""
Unit tests for DefineClones
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from changeo import __version__, __date__

# Imports
import os
import sys
import time
import unittest
import pandas as pd

# Presto and changeo imports
from changeo.Sequence import IgRecord, getDistMat

# Import script
test_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(test_dir, os.pardir, 'bin'))
import DefineClones

# Paths
data_path = os.path.join(test_dir, 'data')
model_path = os.path.join(test_dir, os.pardir, 'models')


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
        #DefineClones.indexJunctions(db_iter, fields=None, mode='gene', action='first')
        self.fail()

    def test_distanceClones(self):

        # import cProfile
        # prof = cProfile.Profile()
        # results = prof.runcall(DefineClones.distanceClones, self.records, model='hs5f', distance=1.0, dist_mat=self.dist_mat)
        # prof.dump_stats('hs5f-unit-test-dict.prof')

        # Define data files
        model_file = os.path.join(model_path, 'HS5F_Distance.tab')
        dist_dict = pd.read_csv(model_file, sep='\t', index_col=0).to_dict()

        results = DefineClones.distanceClones(self.records, model='hs5f', distance=1.0, dist_mat=dist_dict)
        results = {k: sorted(v, key=lambda x: x.id) for k, v in results.iteritems()}
        print 'MODEL> hs5f'
        for k, v in results.iteritems():
            for s in v:
                print '  CLONE-%i> %s' % (k, s.id)

        self.assertEqual(sorted(self.clones.values()), sorted(results.values()))

        # Define data files
        model_file = os.path.join(model_path, 'M3N_Distance.tab')
        dist_dict = pd.read_csv(model_file, sep='\t', index_col=0).to_dict()

        results = DefineClones.distanceClones(self.records, model='m3n', distance=1.0, dist_mat=dist_dict)
        results = {k: sorted(v, key=lambda x: x.id) for k, v in results.iteritems()}
        print 'MODEL> m3n'
        for k, v in results.iteritems():
            for s in v:
                print '  CLONE-%i> %s' % (k, s.id)

        self.assertEqual(sorted(self.clones.values()), sorted(results.values()))

        smith96 = pd.DataFrame([[0,2.86,1,2.14],[2.86,0,2.14,1],[1,2.14,0,2.86],[2.14,1,2.86,0]],
                                index=['A','C','G','T'], columns=['A','C','G','T'], dtype=float)
        dist_dict = getDistMat(smith96)
        results = DefineClones.distanceClones(self.records, model='m1n', distance=10, dist_mat=dist_dict)
        results = {k: sorted(v, key=lambda x: x.id) for k, v in results.iteritems()}
        print 'MODEL> m1n'
        for k, v in results.iteritems():
            for s in v:
                print '  CLONE-%i> %s' % (k, s.id)
        self.assertEqual(sorted(self.clones.values()), sorted(results.values()))

        results = DefineClones.distanceClones(self.records, model='aa', distance=3.0)
        results = {k: sorted(v, key=lambda x: x.id) for k, v in results.iteritems()}
        print 'MODEL> aa'
        for k, v in results.iteritems():
            for s in v:
                print '  CLONE-%i> %s' % (k, s.id)

        self.assertEqual(sorted(self.clones.values()), sorted(results.values()))


if __name__ == '__main__':
    unittest.main()