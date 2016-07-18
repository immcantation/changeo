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

    # @unittest.skip("-> indexJunctions() skipped\n")
    def test_indexJunctions(self):
        from collections import OrderedDict
        db_iter = DefineClones.readDbFile(self.mg_ig_file)

        # Original version of 'set'
        prof_old = cProfile.Profile()
        prof_old.enable()
        results_old = DefineClones.indexJunctions(db_iter, fields=['CREGION'], mode='gene', action='setold')
        prof_old.disable()
        prof_old.print_stats()

        db_iter = DefineClones.readDbFile(self.mg_ig_file)

        # New nested dict version of 'set'
        prof_new = cProfile.Profile()
        prof_new.enable()
        results_new = DefineClones.indexJunctions(db_iter, fields=['CREGION'], mode='gene', action='setnew')
        prof_new.disable()
        prof_new.print_stats()

        # results = {}
        # for k,v in results_new.items():
        #     k_sorted = (k[3], k[2], k[0], k[1])
        #     results[k_sorted] = v
        #
        # results_old = OrderedDict(sorted(results_old.items(), key=lambda x:x[0]))
        # results = OrderedDict(sorted(results.items(), key=lambda x:x[0]))
        # print(results)
        # print(results_old)

        self.fail()

    @unittest.skip("-> distanceClones() skipped\n")
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