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
import numpy as np
from Bio.Seq import Seq

# Presto and changeo imports
import changeo.Distance as Distance

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')


class Test_Distance(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)

        self.nt_seq = ['TGTGCAAGGGGGCCA',
                       'TGTGCAAGGGGGCCA',
                       'TGTATTTGGGGGCCA',
                       'ACACTTGCCACTGAT',
                       'NNNNNNNNNNNNTGA',
                       'NNNNNNNNNNNNNNN']
        self.nt_seq_short = ['TGTGCAAGG',
                             'TGTGCAAGG',
                             'TGTATTTGG',
                             'NNNNNNNNN']
        self.aa_seq = [str(Seq(x).translate()) for x in self.nt_seq]


        # Amino acid distance from sequence 0 to sequence 0,..,5
        self.aa_len = 5.0
        self.aa_ham = np.array([0.0, 0.0, 2.0, 5.0, 1.0, 0.0])

        # Nucleotide distances sequence 0 to sequence 0,..,5
        self.nt_len = 15.0
        # 0 vs 2
        #   A-C = 0; A-G = 1; A-T = 2;
        #   C-G = 0; C-T = 1; G-T = 0
        # 0 vs 3
        #   A-C = 1; A-G = 2; A-T = 4;
        #   C-G = 6; C-T = 1; G-T = 1
        # 0 vs 4
        #   A-C = 0; A-G = 0; A-T = 0;
        #   C-G = 1; C-T = 1; G-T = 0
        self.nt_ham = np.array([0.0, 0.0, 4.0, 15.0, 2.0, 0.0])
        self.nt_hh_s1f = np.array([0.0,
                                   0.0,
                                   0.0 + 1*0.64 + 2*1.16 + 0.0 + 1*0.64 + 0.0,
                                   1*1.21 + 2*0.64 + 4*1.16 + 6*1.16 + 1*0.64 + 1*1.21,
                                   0.0 + 0.0 + 0.0 + 0.0 + 1*1.16 + 1*0.64 + 0.0,
                                   0.0])
        self.nt_mk_rs1nf = np.array([0.0,
                                     0.0,
                                     0.0 + 1*0.32 + 2*1.17 + 0.0 + 1*0.32 + 0.0,
                                     1*1.51 + 2*0.32 + 4*1.17 + 6*1.17 + 1*0.32 + 1*1.51,
                                     0.0 + 0.0 + 0.0 + 0.0 + 1*1.17 + 1*0.32 + 0.0,
                                     0.0])

        # 5-mer models use shorter sequence
        # For sequence 0 to sequence 0,..,3
        self.nt_short_len = 9.0
        # hh_hs5f
        #   GTGCA-GTATT = [0.97, 0.84]
        #   TGCAA-TATTT = [0.93, 0.83]
        #   GCAAG-ATTTG = [1.08, 1.16]
        #   CAAGG-TTTGG = [0.91, 1.07]
        self.nt_hh_s5f_avg = np.array([0.0,
                                       0.0,
                                       np.mean([0.97, 0.84]) + np.mean([0.93, 0.83]) +
                                       np.mean([1.08, 1.16]) + np.mean([0.91, 1.07]),
                                       0.0])
        self.nt_hh_s5f_min = np.array([0.0,
                                       0.0,
                                       np.min([0.97, 0.84]) + np.min([0.93, 0.83]) +
                                       np.min([1.08, 1.16]) + np.min([0.91, 1.07]),
                                       0.0])
        # mk_rs5nf
        #   GTGCA-GTATT = [0.71, 0.77]
        #   TGCAA-TATTT = [0.71, 0.93]
        #   GCAAG-ATTTG = [1.05, 1.03]
        #   CAAGG-TTTGG = [1.08, 1.13]
        self.nt_mk_rs5nf_avg = np.array([0.0,
                                         0.0,
                                         np.mean([0.71, 0.77]) + np.mean([0.71, 0.93]) +
                                         np.mean([1.05, 1.03]) + np.mean([1.08, 1.13]),
                                         0.0])
        self.nt_mk_rs5nf_min = np.array([0.0,
                                         0.0,
                                         np.min([0.71, 0.77]) + np.min([0.71, 0.93]) +
                                         np.min([1.05, 1.03]) + np.min([1.08, 1.13]),
                                         0.0])

        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print("<- %s() %.3f" % (self._testMethodName, t))
        
    def test_calcDistances(self):
        # aa
        print(' AA>', self.aa_seq)
        result = Distance.calcDistances(self.aa_seq, n=1, dist_mat=Distance.aa_model, norm=None)
        print(' AA>', result[1, :])
        np.testing.assert_almost_equal(result[1, :], self.aa_ham)

        result = Distance.calcDistances(self.aa_seq, n=1, dist_mat=Distance.aa_model, norm='len')
        np.testing.assert_almost_equal(result[1, :], self.aa_ham / self.aa_len)

        # ham
        print(' NT>', self.nt_seq)
        result = Distance.calcDistances(self.nt_seq, n=1, dist_mat=Distance.ham_model, norm=None)
        print('HAM>', result[1, :])
        np.testing.assert_almost_equal(result[1, :], self.nt_ham)

        result = Distance.calcDistances(self.nt_seq, n=1, dist_mat=Distance.ham_model, norm='len')
        np.testing.assert_almost_equal(result[1, :], self.nt_ham / self.nt_len)

        # hh_s1f
        result = Distance.calcDistances(self.nt_seq, n=1, dist_mat=Distance.hh_s1f_model, norm=None)
        print('HH_S1F>\n', result[1, :])
        np.testing.assert_almost_equal(result[1, :], self.nt_hh_s1f)

        result = Distance.calcDistances(self.nt_seq, n=1, dist_mat=Distance.hh_s1f_model, norm='len')
        np.testing.assert_almost_equal(result[1, :], self.nt_hh_s1f / self.nt_len)

        # mk_rs1nf
        result = Distance.calcDistances(self.nt_seq, n=1, dist_mat=Distance.mk_rs1nf_model, norm=None)
        print('MK_RS1NF>\n', result[1, :])
        np.testing.assert_almost_equal(result[1, :], self.nt_mk_rs1nf)

        result = Distance.calcDistances(self.nt_seq, n=1, dist_mat=Distance.mk_rs1nf_model, norm='len')
        np.testing.assert_almost_equal(result[1, :], self.nt_mk_rs1nf / self.nt_len)

        # hh_s5f
        result = Distance.calcDistances(self.nt_seq_short, n=5, dist_mat=Distance.hh_s5f_model, sym='avg', norm=None)
        print('HH_S5F-AVG>\n', result[1, :])
        np.testing.assert_almost_equal(result[1, :], self.nt_hh_s5f_avg)

        result = Distance.calcDistances(self.nt_seq_short, n=5, dist_mat=Distance.hh_s5f_model, sym='avg', norm='len')
        np.testing.assert_almost_equal(result[1, :], self.nt_hh_s5f_avg / self.nt_short_len)

        result = Distance.calcDistances(self.nt_seq_short, n=5, dist_mat=Distance.hh_s5f_model, sym='min', norm=None)
        print('HH_S5F-MIN>\n', result[1, :])
        np.testing.assert_almost_equal(result[1, :], self.nt_hh_s5f_min)

        result = Distance.calcDistances(self.nt_seq_short, n=5, dist_mat=Distance.hh_s5f_model, sym='min', norm='len')
        np.testing.assert_almost_equal(result[1, :], self.nt_hh_s5f_min / self.nt_short_len)

        # mk_rs5nf
        result = Distance.calcDistances(self.nt_seq_short, n=5, dist_mat=Distance.mk_rs5nf_model, sym='avg', norm=None)
        print('MK_RS5NF-AVG>\n', result[1, :])
        np.testing.assert_almost_equal(result[1, :], self.nt_mk_rs5nf_avg)

        result = Distance.calcDistances(self.nt_seq_short, n=5, dist_mat=Distance.mk_rs5nf_model, sym='avg', norm='len')
        np.testing.assert_almost_equal(result[1, :], self.nt_mk_rs5nf_avg / self.nt_short_len)

        result = Distance.calcDistances(self.nt_seq_short, n=5, dist_mat=Distance.mk_rs5nf_model, sym='min', norm=None)
        print('MK_RS5NF-MIN>\n', result[1, :])
        np.testing.assert_almost_equal(result[1, :], self.nt_mk_rs5nf_min)

        result = Distance.calcDistances(self.nt_seq_short, n=5, dist_mat=Distance.mk_rs5nf_model, sym='min', norm='len')
        np.testing.assert_almost_equal(result[1, :], self.nt_mk_rs5nf_min / self.nt_short_len)

        # Raise error on length mismatch
        with self.assertRaises(IndexError):
            Distance.calcDistances(['ATGCAA', 'ATGC'], n=1, dist_mat=Distance.ham_model)


    @unittest.skip("-> loadModels() skipped\n")
    def test_loadModels(self):
        print('MK_RS1NF> ', Distance.mk_rs1nf_model)
        print('  HH_S1F> ', Distance.hh_s5f_model)

        self.fail('TODO')


if __name__ == '__main__':
    unittest.main()
