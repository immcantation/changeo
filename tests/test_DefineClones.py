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
from itertools import chain

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
                       'JUNCTION_LENGTH': 48},
                      {'V_CALL': 'IGHV4-1*01',
                       'D_CALL': 'IGHD6-6*01',
                       'J_CALL': 'IGHJ6*01',
                       'JUNCTION_LENGTH': 48},
                      {'V_CALL': 'IGHV2-1*01, IGHV4-1*01',
                       'D_CALL': 'IGHD6-6*01',
                       'J_CALL': 'IGHJ6*01',
                       'JUNCTION_LENGTH': 48},
                      {'V_CALL': 'IGHV5-1*01',
                       'D_CALL': 'IGHD6-6*01',
                       'J_CALL': 'IGHJ6*01',
                       'JUNCTION_LENGTH': 48},
                      {'V_CALL': 'IGHV5-1*01, IGHV6-1*01',
                       'D_CALL': 'IGHD6-6*01',
                       'J_CALL': 'IGHJ6*01',
                       'JUNCTION_LENGTH': 48}]

        # Define unique sequences
        seq_list = [{'SEQUENCE_ID': 'A1',
                     'SEQUENCE_INPUT': 'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGTCTGG',
                     'JUNCTION': 'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGTCTGG'},
                    {'SEQUENCE_ID': 'A2',
                     'SEQUENCE_INPUT': 'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGTCTGG',
                     'JUNCTION': 'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGTCTGG'},
                    {'SEQUENCE_ID': 'A3',
                     'SEQUENCE_INPUT': 'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGTCTGG',
                     'JUNCTION': 'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGTCTGG'},
                    {'SEQUENCE_ID': 'A4',
                     'SEQUENCE_INPUT': 'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGNNNNN',
                     'JUNCTION': 'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGNNNNN'},
                    {'SEQUENCE_ID': 'B1',
                     'SEQUENCE_INPUT': 'TGTGCAAGATATAGCAGCAGCTACTACTACTACGGTATGGACGTCTGG',
                     'JUNCTION': 'TGTGCAAGATATAGCAGCAGCTACTACTACTACGGTATGGACGTCTGG'},
                    {'SEQUENCE_ID': 'B2',
                     'SEQUENCE_INPUT': 'TGTGNAAGATNTAGCAGCAGCTACTACTACTACGGTATNGACGTCTGG',
                     'JUNCTION': 'TGTGNAAGATNTAGCAGCAGCTACTACTACTACGGTATNGACGTCTGG'},
                    {'SEQUENCE_ID': 'B3',
                     'SEQUENCE_INPUT': 'TGTGCAAGATATAGCAGCAGCTACTACTACTACGGTATGGACGTCTGG',
                     'JUNCTION': 'TGTGCAAGATATAGCAGCAGCTACTACTACTACGGTATGGACGTCTGG'},
                    {'SEQUENCE_ID': 'B4',
                     'SEQUENCE_INPUT': 'TGTGNAAGATNTAGCAGCAGCTACTACTACTACGGTATNGACGTCTGG',
                     'JUNCTION': 'TGTGNAAGATNTAGCAGCAGCTACTACTACTACGGTATNGACGTCTGG'}]

        # Build preclone IgRecord list with unambiguous gene calls
        seq_copy = deepcopy(seq_list)
        for x in seq_copy:  x.update(deepcopy(group_list[1]))
        self.unambig_records = [IgRecord(x) for x in seq_copy]
        self.unambig_clones = {('A1', 'A2', 'A3', 'A4'),
                               ('B1', 'B2', 'B3', 'B4')}

        # Build db iterator with ambiguous assignments
        group_copy = deepcopy(group_list)
        seq_copy = deepcopy(seq_list)
        for i, x in enumerate(group_copy):  x.update(seq_copy[i])
        self.ambig_records = [IgRecord(x) for x in group_copy]
        self.ambig_groups = {('48', 'IGHJ6', 'IGHV1-1', 'IGHV2-1', 'IGHV3-1', 'IGHV4-1'):
                             ['A1', 'A2', 'A3', 'A4', 'B1', 'B2'],
                             ('48', 'IGHJ6', 'IGHV5-1', 'IGHV6-1'):
                             ['B3', 'B4']}
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
        results = DefineClones.indexJunctions(self.ambig_records, mode='gene', action='set')

        # Extract nested keys and group lengths for comparison
        results_dict = dict()
        for k, w in results.items():
            print('GROUP>\n', ' ->', k)
            for j, v in w.items():
                print('    ->', j)
                for i, u in v.items():
                    print('      ->', i, ':', len(u))
                    nest_key = tuple(sorted(chain([str(k)], j, i)))
                    results_dict[nest_key] = sorted([x.id for x in u])

        self.assertDictEqual(self.ambig_groups, results_dict)

    # @unittest.skip("-> distanceClones() skipped\n")
    def test_distanceClones(self):
        # import cProfile
        # prof = cProfile.Profile()
        # results = prof.runcall(DefineClones.distanceClones, self.unambig_records, model='hs5f', distance=1.0, dist_mat=self.dist_mat)
        # prof.dump_stats('hs5f-unit-test-dict.prof')

        # ham model
        results = DefineClones.distanceClones(self.unambig_records, model='ham', distance=9.0, norm='none')
        results_set = set([tuple(sorted([x.id for x in v])) for v in results.values()])
        print('MODEL> ham')
        for i, x in enumerate(results_set):  print('  CLONE-%i> %s' % (i + 1, x))
        self.assertSetEqual(results_set, self.unambig_clones)

        # m1n model
        results = DefineClones.distanceClones(self.unambig_records, model='m1n', distance=10.0, norm='none')
        results_set = set([tuple(sorted([x.id for x in v])) for v in results.values()])
        print('MODEL> m1n')
        for i, x in enumerate(results_set):  print('  CLONE-%i> %s' % (i + 1, x))
        self.assertSetEqual(results_set, self.unambig_clones)

        # hs1f model
        results = DefineClones.distanceClones(self.unambig_records, model='hs1f', distance=0.25)
        results_set = set([tuple(sorted([x.id for x in v])) for v in results.values()])
        print('MODEL> hs1f')
        for i, x in enumerate(results_set):  print('  CLONE-%i> %s' % (i + 1, x))
        self.assertSetEqual(results_set, self.unambig_clones)

        # hs5f model
        results = DefineClones.distanceClones(self.unambig_records, model='hs5f', distance=0.1)
        results_set = set([tuple(sorted([x.id for x in v])) for v in results.values()])
        print('MODEL> hs5f')
        for i, x in enumerate(results_set):  print('  CLONE-%i> %s' % (i + 1, x))
        self.assertSetEqual(results_set, self.unambig_clones)

        # aa model
        results = DefineClones.distanceClones(self.unambig_records, model='aa', distance=3.0, norm='none')
        results_set = set([tuple(sorted([x.id for x in v])) for v in results.values()])
        print('MODEL> aa')
        for i, x in enumerate(results_set):  print('  CLONE-%i> %s' % (i + 1, x))
        self.assertSetEqual(results_set, self.unambig_clones)



if __name__ == '__main__':
    unittest.main()