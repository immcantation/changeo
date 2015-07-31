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

# Presto and changeo imports
from changeo.Receptor import IgRecord

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')

# Import script
sys.path.append(os.path.join(test_path, os.pardir, 'bin'))
import DefineClones


class Test_DefineClones(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)

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
        print("<- %s() %.3f" % (self._testMethodName, t))

    @unittest.skip("-> indexJunctions() skipped\n")
    def test_indexJunctions(self):
        #DefineClones.indexJunctions(db_iter, fields=None, mode='gene', action='first')
        self.fail()

    def test_distanceClones(self):
        # import cProfile
        # prof = cProfile.Profile()
        # results = prof.runcall(DefineClones.distanceClones, self.records, model='hs5f', distance=1.0, dist_mat=self.dist_mat)
        # prof.dump_stats('hs5f-unit-test-dict.prof')

        # hs5f model
        results = DefineClones.distanceClones(self.records, model='hs5f', distance=0.1)
        results = {k: sorted(v, key=lambda x: x.id) for k, v in results.items()}
        print('MODEL> hs5f')
        for k, v in results.items():
            for s in v:
                print('  CLONE-%i> %s' % (k, s.id))
        self.assertListEqual(sorted(self.clones.values(), key=lambda x: len(x)),
                             sorted(results.values(), key=lambda x: len(x)))

        # m1n model
        results = DefineClones.distanceClones(self.records, model='m1n', distance=10.0, norm='none')
        results = {k: sorted(v, key=lambda x: x.id) for k, v in results.items()}
        print('MODEL> m1n')
        for k, v in results.items():
            for s in v:
                print('  CLONE-%i> %s' % (k, s.id))
        self.assertListEqual(sorted(self.clones.values(), key=lambda x: len(x)),
                             sorted(results.values(), key=lambda x: len(x)))

        # hs1f model
        results = DefineClones.distanceClones(self.records, model='hs1f', distance=0.25)
        results = {k: sorted(v, key=lambda x: x.id) for k, v in results.items()}
        print('MODEL> hs1f')
        for k, v in results.items():
            for s in v:
                print('  CLONE-%i> %s' % (k, s.id))
        self.assertListEqual(sorted(self.clones.values(), key=lambda x: len(x)),
                             sorted(results.values(), key=lambda x: len(x)))

        # aa model
        results = DefineClones.distanceClones(self.records, model='aa', distance=3.0, norm='none')
        results = {k: sorted(v, key=lambda x: x.id) for k, v in results.items()}
        print('MODEL> aa')
        for k, v in results.items():
            for s in v:
                print('  CLONE-%i> %s' % (k, s.id))
        self.assertListEqual(sorted(self.clones.values(), key=lambda x: len(x)),
                             sorted(results.values(), key=lambda x: len(x)))

         # ham model
        results = DefineClones.distanceClones(self.records, model='ham', distance=9.0, norm='none')
        results = {k: sorted(v, key=lambda x: x.id) for k, v in results.items()}
        print('MODEL> ham')
        for k, v in results.items():
            for s in v:
                print('  CLONE-%i> %s' % (k, s.id))
        self.assertListEqual(sorted(self.clones.values(), key=lambda x: len(x)),
                             sorted(results.values(), key=lambda x: len(x)))


if __name__ == '__main__':
    unittest.main()