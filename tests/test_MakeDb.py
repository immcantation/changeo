"""
Unit tests for MakeDb
"""

# Info
__author__ = 'Gisela Gabernet'
from changeo import __version__, __date__

# Imports
import os
import sys
import time
import unittest
import airr

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')

# Import script
sys.path.append(os.path.join(test_path, os.pardir, 'bin'))
import MakeDb

class Test_MakeDb(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)
        self.start = time.time()

        self.testdata_indels = os.path.join(data_path, 'makedb_number_test.tsv')

        # Download reference data running wget command
        os.system('wget -c https://github.com/nf-core/test-datasets/raw/airrflow/database-cache/imgtdb_base.zip')
        os.system('unzip -o imgtdb_base.zip -d ' + data_path)
        os.system('rm imgtdb_base.zip')


    def tearDown(self):
        t = time.time() - self.start
        print("<- %s() %.3f" % (self._testMethodName, t))

    def test_numberAIRR(self):
        """Test IMGT numbering insertion"""
        out_file = os.path.join(data_path, 'makedb_number_test_out.tsv')
        MakeDb.numberAIRR(self.testdata_indels,
                          germline_reference=[os.path.join(data_path, 'imgtdb_base/human/vdj')],
                          out_file=out_file, debug_mode=True)
        self.assertTrue(os.path.isfile(out_file))

        out_tab = airr.load_rearrangement(out_file)
        self.assertTrue(out_tab.shape[0] == 6)

        #TODO: potentially add more tests.

        #os.remove(out_file)
        os.system('rm -r ' + os.path.join(data_path, 'imgtdb_base'))


if __name__ == '__main__':
    unittest.main()

