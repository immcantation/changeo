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

# Presto and changeo imports
import changeo.Distance as Distance
# Paths
test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')


class Test_Distance(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)
        
        self.seq_list = ['TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGTCTGG',
                         'TGTAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGTCTGG']
                         
        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print("<- %s() %.3f" % (self._testMethodName, t))
        
    def test_calcDistances(self):
	    result_mat = Distance.calcDistances(self.seq_list, n = 1, dist_mat = Distance.hs1f_model,
                                            norm = 'len', sym = 'min')
	    print(result_mat)

    @unittest.skip("-> loadModels() skipped\n")
    def test_loadModels(self):
        print(' M1N> ', Distance.m1n_model)
        print('HS5F> ', Distance.hs5f_model)

        self.fail('TODO')


if __name__ == '__main__':
    unittest.main()
