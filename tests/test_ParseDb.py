"""
Unit tests for ParseDb
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from changeo import __version__, __date__

# Imports
import unittest
import time


class Test_ParseDb(unittest.TestCase):
    def setUp(self):
        print '-> %s()' % self._testMethodName
        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print "<- %s() %.3f" % (self._testMethodName, t)


if __name__ == '__main__':
    unittest.main()