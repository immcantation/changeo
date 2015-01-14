"""
Unit tests for DefineClones
"""

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.0'
__date__      = '2014.11.26'

# Imports
import time, unittest
import DefineClones as mod


class Test_DefineClones(unittest.TestCase):
    def setUp(self):
        print '-> %s()' % self._testMethodName
        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print "<- %s() %.3f" % (self._testMethodName, t)

    def test_indexJunctions(self):
        #mod.indexJunctions(db_iter, fields=None, mode='gene', action='first')
        self.fail()


if __name__ == '__main__':
    unittest.main()