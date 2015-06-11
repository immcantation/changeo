"""
Unit tests for Distance module
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
from changeo.Distance import getDNADistMatrix, getAADistMatrix, \
                             m1n_distance, m3n_distance, hs5f_distance
# Paths
test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')


class Test_Distance(unittest.TestCase):
    def setUp(self):
        print '-> %s()' % self._testMethodName

        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print "<- %s() %.3f" % (self._testMethodName, t)

    #@unittest.skip("-> loadModels() skipped\n")
    def test_loadModels(self):
        print ' M1N> ', m1n_distance
        print ' M3N> ', m3n_distance
        print 'HS5F> ', hs5f_distance

        self.fail('TODO')


if __name__ == '__main__':
    unittest.main()