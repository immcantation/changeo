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
import cProfile
from collections import OrderedDict
from copy import deepcopy

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

        self.mg_ig_file = os.path.join(data_path, 'AB0RF_MG67U_functional-heavy_parse-select_filter-collapse.tab')

        # Define list of preclone properties
        group_list = [{'V_CALL': 'IGHV1-1*01',
                       'D_CALL': 'IGHD6-6*01',
                       'J_CALL': 'IGHJ6*01',
                       'JUNCTION_LENGTH': 48},
                      {'V_CALL': 'IGHV2-1*01',
                       'D_CALL': 'IGHD6-6*01',
                       'J_CALL': 'IGHJ6*01',
                       'JUNCTION_LENGTH': 48},
                      {'V_CALL': 'IGHV3-1*01',
                       'D_CALL': 'IGHD6-6*01',
                       'J_CALL': 'IGHJ6*01',
                       'JUNCTION_LENGTH': 48},
                      {'V_CALL': 'IGHV1-1*01, IGHV2-1*01, IGHV3-1*01',
                       'D_CALL': 'IGHD6-6*01',
                       'J_CALL': 'IGHJ6*01',
                       'JUNCTION_LENGTH': 48}]

        # Define unique sequences
        seq_list = [{'SEQUENCE_ID':'A1',
                     'SEQUENCE_INPUT':'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGTCTGG',
                     'JUNCTION':'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGTCTGG'},
                    {'SEQUENCE_ID':'A2',
                     'SEQUENCE_INPUT':'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGTCTGG',
                     'JUNCTION':'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGTCTGG'},
                    {'SEQUENCE_ID':'A3',
                     'SEQUENCE_INPUT':'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGTCTGG',
                     'JUNCTION':'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGTCTGG'},
                    {'SEQUENCE_ID':'A4',
                     'SEQUENCE_INPUT':'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGNNNNN',
                     'JUNCTION':'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGNNNNN'},
                    {'SEQUENCE_ID':'B1',
                     'SEQUENCE_INPUT':'TGTGCAAGATATAGCAGCAGCTACTACTACTACGGTATGGACGTCTGG',
                     'JUNCTION':'TGTGCAAGATATAGCAGCAGCTACTACTACTACGGTATGGACGTCTGG'},
                    {'SEQUENCE_ID':'B2',
                     'SEQUENCE_INPUT':'TGTGNAAGATNTAGCAGCAGCTACTACTACTACGGTATNGACGTCTGG',
                     'JUNCTION':'TGTGNAAGATNTAGCAGCAGCTACTACTACTACGGTATNGACGTCTGG'}]

        # Build preclone IgRecord list with unambiguous gene calls
        seq_copy = deepcopy(seq_list)
        for x in seq_copy:  x.update(deepcopy(group_list[1]))
        self.unambig_records = [IgRecord(x) for x in seq_copy]
        self.unambig_clones = {1:self.unambig_records[0:4], 2:self.unambig_records[4:]}

        # Build db iterator with ambiguous assignments
        group_copy = deepcopy(group_list)
        seq_copy = deepcopy(seq_list)
        for i, x in enumerate(group_copy):  x.update(seq_copy[i])
        self.ambig_records = [IgRecord(x) for x in group_copy]

        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print("<- %s() %.3f" % (self._testMethodName, t))

    # @unittest.skip("-> indexJunctions() skipped\n")
    def test_indexJunctions(self):
        # db_iter = DefineClones.readDbFile(self.mg_ig_file)
        #
        # # Original version of 'set'
        # prof_old = cProfile.Profile()
        # prof_old.enable()
        # results_old = DefineClones.indexJunctions(db_iter, fields=['CREGION'], mode='gene', action='setold')
        # prof_old.disable()
        # prof_old.print_stats()
        #
        # db_iter = DefineClones.readDbFile(self.mg_ig_file)
        #
        # # New nested dict version of 'set'
        # prof_new = cProfile.Profile()
        # prof_new.enable()
        # results_new = DefineClones.indexJunctions(db_iter, fields=['CREGION'], mode='gene', action='set')
        # prof_new.disable()
        # prof_new.print_stats()

        # results = {}
        # for k,v in results_new.items():
        #     k_sorted = (k[3], k[2], k[0], k[1])
        #     results[k_sorted] = v
        #
        # results_old = OrderedDict(sorted(results_old.items(), key=lambda x:x[0]))
        # results = OrderedDict(sorted(results.items(), key=lambda x:x[0]))
        # print(results)
        # print(results_old)

        # Test ambiguous grouping
        #for x in self.ambig_records: print(x.id)
        results = DefineClones.indexJunctions(self.ambig_records, mode='gene', action='set')
        #print(results)
        #for k, v in results.items(): print('%s: %d' % (k, len(v)))
        for k, w in results.items():
            print(k, ' -> ', end='')
            for j, v in w.items():
                print(j, ' -> ', end='')
                for i, u in v.items():
                    print(i, ' : ', len(u))

        self.fail()

    @unittest.skip("-> distanceClones() skipped\n")
    def test_distanceClones(self):
        # import cProfile
        # prof = cProfile.Profile()
        # results = prof.runcall(DefineClones.distanceClones, self.unambig_records, model='hs5f', distance=1.0, dist_mat=self.dist_mat)
        # prof.dump_stats('hs5f-unit-test-dict.prof')

        # hs5f model
        results = DefineClones.distanceClones(self.unambig_records, model='hs5f', distance=0.1)
        results = {k: sorted(v, key=lambda x: x.id) for k, v in results.items()}
        print('MODEL> hs5f')
        for k, v in results.items():
            for s in v:
                print('  CLONE-%i> %s' % (k, s.id))
        self.assertListEqual(sorted(self.unambig_clones.values(), key=lambda x: len(x)),
                             sorted(results.values(), key=lambda x: len(x)))

        # m1n model
        results = DefineClones.distanceClones(self.unambig_records, model='m1n', distance=10.0, norm='none')
        results = {k: sorted(v, key=lambda x: x.id) for k, v in results.items()}
        print('MODEL> m1n')
        for k, v in results.items():
            for s in v:
                print('  CLONE-%i> %s' % (k, s.id))
        self.assertListEqual(sorted(self.unambig_clones.values(), key=lambda x: len(x)),
                             sorted(results.values(), key=lambda x: len(x)))

        # hs1f model
        results = DefineClones.distanceClones(self.unambig_records, model='hs1f', distance=0.25)
        results = {k: sorted(v, key=lambda x: x.id) for k, v in results.items()}
        print('MODEL> hs1f')
        for k, v in results.items():
            for s in v:
                print('  CLONE-%i> %s' % (k, s.id))
        self.assertListEqual(sorted(self.unambig_clones.values(), key=lambda x: len(x)),
                             sorted(results.values(), key=lambda x: len(x)))

        # aa model
        results = DefineClones.distanceClones(self.unambig_records, model='aa', distance=3.0, norm='none')
        results = {k: sorted(v, key=lambda x: x.id) for k, v in results.items()}
        print('MODEL> aa')
        for k, v in results.items():
            for s in v:
                print('  CLONE-%i> %s' % (k, s.id))
        self.assertListEqual(sorted(self.unambig_clones.values(), key=lambda x: len(x)),
                             sorted(results.values(), key=lambda x: len(x)))

         # ham model
        results = DefineClones.distanceClones(self.unambig_records, model='ham', distance=9.0, norm='none')
        results = {k: sorted(v, key=lambda x: x.id) for k, v in results.items()}
        print('MODEL> ham')
        for k, v in results.items():
            for s in v:
                print('  CLONE-%i> %s' % (k, s.id))
        self.assertListEqual(sorted(self.unambig_clones.values(), key=lambda x: len(x)),
                             sorted(results.values(), key=lambda x: len(x)))


if __name__ == '__main__':
    unittest.main()