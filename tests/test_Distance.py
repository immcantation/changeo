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
from changeo.Distance import m1n_model, hs5f_model
# Paths
test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')


class Test_Distance(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)

        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print("<- %s() %.3f" % (self._testMethodName, t))

    @unittest.skip("-> loadModels() skipped\n")
    def test_loadModels(self):
        print(' M1N> ', m1n_model)
        print('HS5F> ', hs5f_model)

        self.fail('TODO')


if __name__ == '__main__':
    unittest.main()