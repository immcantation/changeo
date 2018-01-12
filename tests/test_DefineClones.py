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
from copy import deepcopy
from itertools import chain

# Presto and changeo imports
from changeo.Receptor import Receptor
from changeo.Multiprocessing import DbResult

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
                       'JUNCTION_LENGTH': 48,
                       'FIELD1': 'GA',
                       'FIELD2': '1'},
                      {'V_CALL': 'IGHV2-1*01',
                       'D_CALL': 'IGHD6-6*01',
                       'J_CALL': 'IGHJ6*01',
                       'JUNCTION_LENGTH': 48,
                       'FIELD1': 'GA',
                       'FIELD2': '1'},
                      {'V_CALL': 'IGHV3-1*01',
                       'D_CALL': 'IGHD6-6*01',
                       'J_CALL': 'IGHJ6*01',
                       'JUNCTION_LENGTH': 48,
                       'FIELD1': 'GA',
                       'FIELD2': '2'},
                      {'V_CALL': 'IGHV1-1*01, IGHV2-1*01, IGHV3-1*01',
                       'D_CALL': 'IGHD6-6*01',
                       'J_CALL': 'IGHJ6*01',
                       'JUNCTION_LENGTH': 48,
                       'FIELD1': 'GA',
                       'FIELD2': '2'},
                      {'V_CALL': 'IGHV4-1*01',
                       'D_CALL': 'IGHD6-6*01',
                       'J_CALL': 'IGHJ6*01',
                       'JUNCTION_LENGTH': 48,
                       'FIELD1': 'GB',
                       'FIELD2': '1'},
                      {'V_CALL': 'IGHV2-1*01, IGHV4-1*01',
                       'D_CALL': 'IGHD6-6*01',
                       'J_CALL': 'IGHJ6*01',
                       'JUNCTION_LENGTH': 48,
                       'FIELD1': 'GB',
                       'FIELD2': '2'},
                      {'V_CALL': 'IGHV5-1*01',
                       'D_CALL': 'IGHD6-6*01',
                       'J_CALL': 'IGHJ6*01',
                       'JUNCTION_LENGTH': 48,
                       'FIELD1': 'GB',
                       'FIELD2': None},
                      {'V_CALL': 'IGHV5-1*01, IGHV6-1*01',
                       'D_CALL': 'IGHD6-6*01',
                       'J_CALL': 'IGHJ6*01',
                       'JUNCTION_LENGTH': 48,
                       'FIELD1': '',
                       'FIELD2': '1'}]

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

        # Build preclone Receptor list with unambiguous gene calls
        seq_copy = deepcopy(seq_list)
        for x in seq_copy:  x.update(deepcopy(group_list[1]))
        self.unambig_records = [Receptor(x) for x in seq_copy]
        self.unambig_data = DbResult('unambig', self.unambig_records)
        self.unambig_clones = {('A1', 'A2', 'A3', 'A4'),
                               ('B1', 'B2', 'B3', 'B4')}

        # Build db iterator with ambiguous assignments
        group_copy = deepcopy(group_list)
        seq_copy = deepcopy(seq_list)
        for i, x in enumerate(group_copy):  x.update(seq_copy[i])
        self.ambig_records = [Receptor(x) for x in group_copy]
        self.ambig_data = DbResult('ambig', self.ambig_records)
        self.first_nofields = {('IGHV1-1', 'IGHJ6', '48'): ['A1','A4'],
                               ('IGHV2-1', 'IGHJ6', '48'): ['A2','B2'],
                               ('IGHV3-1', 'IGHJ6', '48'): ['A3'],
                               ('IGHV4-1', 'IGHJ6', '48'): ['B1'],
                               ('IGHV5-1', 'IGHJ6', '48'): ['B3','B4']}
        self.set_nofields = {('48', ('IGHJ6',), ('IGHV1-1', 'IGHV2-1', 'IGHV3-1', 'IGHV4-1')):
                             ['A1', 'A2', 'A3', 'A4', 'B1', 'B2'],
                             ('48', ('IGHJ6',), ('IGHV5-1', 'IGHV6-1')):
                             ['B3', 'B4']}
        self.first_fields = {('IGHV1-1', 'IGHJ6', '48', 'GA', '1'): ['A1'],
                             ('IGHV2-1', 'IGHJ6', '48', 'GA', '1'): ['A2'],
                             ('IGHV3-1', 'IGHJ6', '48', 'GA', '2'): ['A3'],
                             ('IGHV1-1', 'IGHJ6', '48', 'GA', '2'): ['A4'],
                             ('IGHV4-1', 'IGHJ6', '48', 'GB', '1'): ['B1'],
                             ('IGHV2-1', 'IGHJ6', '48', 'GB', '2'): ['B2'],
                             None: ['B3','B4']}
        self.set_fields = {('48', 'GA', '1', ('IGHJ6',), ('IGHV1-1',)): ['A1'],
                           ('48', 'GA', '1', ('IGHJ6',), ('IGHV2-1',)): ['A2'],
                           ('48', 'GA', '2', ('IGHJ6',), ('IGHV1-1', 'IGHV2-1', 'IGHV3-1')): ['A3', 'A4'],
                           ('48', 'GB', '1', ('IGHJ6',), ('IGHV4-1',)): ['B1'],
                           ('48', 'GB', '2', ('IGHJ6',), ('IGHV2-1', 'IGHV4-1')): ['B2'],
                           None: ['B3','B4']}
        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print("<- %s() %.3f" % (self._testMethodName, t))

    # @unittest.skip("-> groupByGene() skipped\n")
    def test_groupByGene(self):
        # Test first grouping without fields
        results = DefineClones.groupByGene(self.ambig_records, mode='gene', action='first')
        # Extract nested keys and group lengths for comparison
        results_dict = dict()
        print('FIRST>')
        for k, v in results.items():
            nest_key = tuple([tuple(sorted(chain(x))) if isinstance(x, tuple) else str(x) for x in k])
            results_dict[nest_key] = sorted([x.sequence_id for x in v])
            print('  GROUP>', nest_key, ':', results_dict[nest_key])
        print('')
        self.assertDictEqual(self.first_nofields, results_dict)

        # Test ambiguous grouping without fields
        results = DefineClones.groupByGene(self.ambig_records, mode='gene', action='set')
        # Extract nested keys and group lengths for comparison
        results_dict = dict()
        print('SET>')
        for k, v in results.items():
            nest_key = tuple([tuple(sorted(chain(x))) if isinstance(x, tuple) else str(x) for x in k])
            results_dict[nest_key] = sorted([x.sequence_id for x in v])
            print('  GROUP>', nest_key, ':', results_dict[nest_key])
        print('')
        self.assertDictEqual(self.set_nofields, results_dict)

        # Test first grouping with fields
        results = DefineClones.groupByGene(self.ambig_records, group_fields=['FIELD1', 'FIELD2'],
                                           mode='gene', action='first')
        # Extract nested keys and group lengths for comparison
        results_dict = dict()
        print('FIRST>')
        for k, v in results.items():
            nest_key = tuple([tuple(sorted(chain(x))) if isinstance(x, tuple) else str(x) for x in k]) if k is not None else None
            results_dict[nest_key] = sorted([x.sequence_id for x in v])
            print('  GROUP>', nest_key, ':', results_dict[nest_key])
        print('')
        print('RRRR>', self.first_fields)
        print('QQQQ>', results_dict)

        self.assertDictEqual(self.first_fields, results_dict)

        # Test ambiguous grouping with fields
        results = DefineClones.groupByGene(self.ambig_records, group_fields=['FIELD1', 'FIELD2'],
                                              mode='gene', action='set')
        # Extract nested keys and group lengths for comparison
        results_dict = dict()
        print('SET>')
        for k, v in results.items():
            nest_key = tuple([tuple(sorted(chain(x))) if isinstance(x, tuple) else str(x) for x in k]) if k is not None else None
            results_dict[nest_key] = sorted([x.sequence_id for x in v])
            print('  GROUP>', nest_key, ':', results_dict[nest_key])
        print('')
        self.assertDictEqual(self.set_fields, results_dict)


    # @unittest.skip("-> distanceClones() skipped\n")
    def test_distanceClones(self):
        # import cProfile
        # prof = cProfile.Profile()
        # results = prof.runcall(DefineClones.distanceClones, self.unambig_records, model='hs5f', distance=1.0, dist_mat=self.dist_mat)
        # prof.dump_stats('hs5f-unit-test-dict.prof')

        # ham model
        results = DefineClones.distanceClones(self.unambig_data, model='ham', distance=9.0, norm='none')
        results_set = set([tuple(sorted([x.sequence_id for x in v])) for v in results.results.values()])
        print('MODEL> ham')
        for i, x in enumerate(results_set):  print('  CLONE-%i> %s' % (i + 1, x))
        self.assertSetEqual(results_set, self.unambig_clones)

        # m1n_compat model
        results = DefineClones.distanceClones(self.unambig_data, model='m1n_compat', distance=10.0, norm='none')
        results_set = set([tuple(sorted([x.sequence_id for x in v])) for v in results.results.values()])
        print('MODEL> m1n_compat')
        for i, x in enumerate(results_set):  print('  CLONE-%i> %s' % (i + 1, x))
        self.assertSetEqual(results_set, self.unambig_clones)

        # hs1f_compat model
        results = DefineClones.distanceClones(self.unambig_data, model='hs1f_compat', distance=0.25)
        results_set = set([tuple(sorted([x.sequence_id for x in v])) for v in results.results.values()])
        print('MODEL> hs1f_compat')
        for i, x in enumerate(results_set):  print('  CLONE-%i> %s' % (i + 1, x))
        self.assertSetEqual(results_set, self.unambig_clones)

        # hh_s5f model
        results = DefineClones.distanceClones(self.unambig_data, model='hh_s5f', distance=0.1)
        results_set = set([tuple(sorted([x.sequence_id for x in v])) for v in results.results.values()])
        print('MODEL> hh_s5f')
        for i, x in enumerate(results_set):  print('  CLONE-%i> %s' % (i + 1, x))
        self.assertSetEqual(results_set, self.unambig_clones)

        # aa model
        results = DefineClones.distanceClones(self.unambig_data, model='aa', distance=3.0, norm='none')
        results_set = set([tuple(sorted([x.sequence_id for x in v])) for v in results.results.values()])
        print('MODEL> aa')
        for i, x in enumerate(results_set):  print('  CLONE-%i> %s' % (i + 1, x))
        self.assertSetEqual(results_set, self.unambig_clones)


if __name__ == '__main__':
    unittest.main()