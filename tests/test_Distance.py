"""
Unit tests for Distance module
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from changeo import __version__, __date__

# Imports
import os
import time
import unittest
from Bio.Seq import Seq

# Presto and changeo imports
import changeo.Distance as Distance
# Paths
test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')


class Test_Distance(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)
        
        self.seq_error = ['TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGTCTGG',
                          'TGTAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGTCTGG']
        self.seq_valid = ['TGTGCAAGGGGGCCA',
                          'TGTGCAAGGGGGCCA',
                          'TGTATTTGGGGGCCA',
                          'ACACGTTCCCCCGGT',
                          'NNNNNNNNNNNNNNN']
        self.aa_valid = [str(Seq(x).translate()) for x in self.seq_valid]

        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print("<- %s() %.3f" % (self._testMethodName, t))
        
    def test_calcDistances(self):
        # aa
        result = Distance.calcDistances(self.aa_valid, n=1, dist_mat=Distance.aa_model,
                                        norm='len', sym='min')
        print('    AA>\n', result)

        # ham
        result = Distance.calcDistances(self.seq_valid, n=1, dist_mat=Distance.ham_model,
                                        norm='len', sym='min')
        print('   HAM>\n', result)

        # hh_s1f
        result = Distance.calcDistances(self.seq_valid, n=1, dist_mat=Distance.hh_s1f_model,
                                        norm='len', sym='min')
        print('HH_S1F>\n', result)

        # hh_s5f
        result = Distance.calcDistances(self.seq_valid, n=5, dist_mat=Distance.hh_s5f_model,
                                        norm='len', sym='min')
        print('HH_S5F>\n', result)

        # Raise error on length mismatch
        with self.assertRaises(IndexError):
            Distance.calcDistances(self.seq_error, n=1, dist_mat=Distance.hh_s1f_model,
                                   norm='len', sym='min')


    @unittest.skip("-> loadModels() skipped\n")
    def test_loadModels(self):
        print('MK_RS1NF> ', Distance.mk_rs1nf_model)
        print('  HH_S1F> ', Distance.hh_s5f_model)

        self.fail('TODO')


if __name__ == '__main__':
    unittest.main()
