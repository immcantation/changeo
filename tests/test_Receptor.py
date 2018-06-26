"""
Unit tests for MakeDb
"""

# Imports
import os
import time
import unittest

# Presto and changeo imports
from changeo.Receptor import Receptor
from changeo.Alignment import decodeBTOP, decodeCIGAR, encodeCIGAR, padAlignment, alignmentPositions

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')


class Test_Receptor(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)

        row = {'SEQUENCE_ID':    'TEST1',
               'SEQUENCE_INPUT': 'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGNNNNN',
               'V_CALL':         'IGHV6-1*01,IGHV6-1*02',
               'D_CALL':         'IGHD6-6*01',
               'J_CALL':         'IGHJ6*02',}
        self.ig_rec = Receptor(row)

        row = {'SEQUENCE_ID':    'TEST2',
               'SEQUENCE_INPUT': 'TGTGCAAGGGGGCCATTGGACTACTTCTACTACGGTGTGGACGNNNNN',
               'V_CALL':         'TRAV6-1*01,TRAV6-1*02',
               'D_CALL':         'TRAD6-6*01',
               'J_CALL':         'TRAJ6*02',}
        self.tr_rec = Receptor(row)

        row = {'SEQUENCE_ID':    'TEST3',
               'SEQUENCE_INPUT': '',
               'V_CALL':         '',
               'D_CALL':         None,
               'J_CALL':         'princess,banana,hammock',}
        self.bad_rec = Receptor(row)

        # CIGAR strings
        self.cigar_string = ['30M1I69M3D', '4S16N30M1I69M3D']
        self.cigar_decoded = [[('M', 30), ('I', 1), ('M', 69), ('D', 3)],
                              [('S', 4), ('N', 16), ('M', 30), ('I', 1), ('M', 69), ('D', 3)]]
        self.cigar_len = [{'q_start': 0, 'q_length': 30 + 1 + 69 + 0, 'r_start': 0, 'r_length': 30 + 0 + 69 + 3},
                          {'q_start': 4, 'q_length': 30 + 1 + 69 + 0, 'r_start': 16, 'r_length': 30 + 0 + 69 + 3}]

        self.cigar_pos = [(0, 0), (5, 0), (0, 3), (5, 3)]
        self.cigar_pad_1 = [[('M', 30), ('I', 1), ('M', 69), ('D', 3)],
                            [('S', 5), ('M', 30), ('I', 1), ('M', 69), ('D', 3)],
                            [('N', 3), ('M', 30), ('I', 1), ('M', 69), ('D', 3)],
                            [('S', 5), ('N', 3), ('M', 30), ('I', 1), ('M', 69), ('D', 3)]]
        self.cigar_pad_2 = [[('S', 4), ('N', 16), ('M', 30), ('I', 1), ('M', 69), ('D', 3)],
                            [('S', 9), ('N', 16), ('M', 30), ('I', 1), ('M', 69), ('D', 3)],
                            [('S', 4), ('N', 19), ('M', 30), ('I', 1), ('M', 69), ('D', 3)],
                            [('S', 9), ('N', 19), ('M', 30), ('I', 1), ('M', 69), ('D', 3)]]

        # BTOP strings
        self.btop_string = ['7AGAC39',
                            '7A-39',
                            '6-G-A41',
                            'AG8-GC-CTCT']
        self.btop_decoded_full = [[('=', 7), ('X', 2), ('=', 39)],
                                  [('=', 7), ('D', 1), ('=', 39)],
                                  [('=', 6), ('I', 2), ('=', 41)],
                                  [('X', 1), ('=', 8), ('I', 1), ('D', 1), ('X', 2)]]

        self.start = time.time()

    def tearDown(self):
        t = time.time() - self.start
        print('<- %s() %.3f' % (self._testMethodName, t))

    @unittest.skip("-> Receptor() skipped\n")
    def test_Receptor(self):
        print('IG>')
        print(self.ig_rec.getAlleleCalls(['v','d','j'], action='first'))
        print(self.ig_rec.getAlleleCalls(['v','j'], action='first'))
        print(self.ig_rec.getAlleleCalls(['j','v'], action='first'))
        print(self.ig_rec.getGeneCalls(['v','d','j'], action='first'))
        print(self.ig_rec.getGeneCalls(['v','j'], action='first'))
        print(self.ig_rec.getGeneCalls(['j','v'], action='first'))
        print(self.ig_rec.getFamilyCalls(['v','d','j'], action='first'))
        print(self.ig_rec.getFamilyCalls(['v','j'], action='first'))
        print(self.ig_rec.getFamilyCalls(['j','v'], action='first'))

        print(self.ig_rec.getAlleleCalls(['v','d','j'], action='set'))
        print(self.ig_rec.getAlleleCalls(['v','j'], action='set'))
        print(self.ig_rec.getAlleleCalls(['j','v'], action='set'))
        print(self.ig_rec.getGeneCalls(['v','d','j'], action='set'))
        print(self.ig_rec.getGeneCalls(['v','j'], action='set'))
        print(self.ig_rec.getGeneCalls(['j','v'], action='set'))
        print(self.ig_rec.getFamilyCalls(['v','d','j'], action='set'))
        print(self.ig_rec.getFamilyCalls(['v','j'], action='set'))
        print(self.ig_rec.getFamilyCalls(['j','v'], action='set'))

        print(self.ig_rec.getAlleleCalls(['v','d','j'], action='list'))
        print(self.ig_rec.getAlleleCalls(['v','j'], action='list'))
        print(self.ig_rec.getAlleleCalls(['j','v'], action='list'))
        print(self.ig_rec.getGeneCalls(['v','d','j'], action='list'))
        print(self.ig_rec.getGeneCalls(['v','j'], action='list'))
        print(self.ig_rec.getGeneCalls(['j','v'], action='list'))
        print(self.ig_rec.getFamilyCalls(['v','d','j'], action='list'))
        print(self.ig_rec.getFamilyCalls(['v','j'], action='list'))
        print(self.ig_rec.getFamilyCalls(['j','v'], action='list'))

        print('TR>')
        print(self.tr_rec.getAlleleCalls(['v','d','j'], action='first'))
        print(self.tr_rec.getAlleleCalls(['v','j'], action='first'))
        print(self.tr_rec.getAlleleCalls(['j','v'], action='first'))
        print(self.tr_rec.getGeneCalls(['v','d','j'], action='first'))
        print(self.tr_rec.getGeneCalls(['v','j'], action='first'))
        print(self.tr_rec.getGeneCalls(['j','v'], action='first'))
        print(self.tr_rec.getFamilyCalls(['v','d','j'], action='first'))
        print(self.tr_rec.getFamilyCalls(['v','j'], action='first'))
        print(self.tr_rec.getFamilyCalls(['j','v'], action='first'))

        print(self.tr_rec.getAlleleCalls(['v','d','j'], action='set'))
        print(self.tr_rec.getAlleleCalls(['v','j'], action='set'))
        print(self.tr_rec.getAlleleCalls(['j','v'], action='set'))
        print(self.tr_rec.getGeneCalls(['v','d','j'], action='set'))
        print(self.tr_rec.getGeneCalls(['v','j'], action='set'))
        print(self.tr_rec.getGeneCalls(['j','v'], action='set'))
        print(self.tr_rec.getFamilyCalls(['v','d','j'], action='set'))
        print(self.tr_rec.getFamilyCalls(['v','j'], action='set'))
        print(self.tr_rec.getFamilyCalls(['j','v'], action='set'))

        print(self.tr_rec.getAlleleCalls(['v','d','j'], action='list'))
        print(self.tr_rec.getAlleleCalls(['v','j'], action='list'))
        print(self.tr_rec.getAlleleCalls(['j','v'], action='list'))
        print(self.tr_rec.getGeneCalls(['v','d','j'], action='list'))
        print(self.tr_rec.getGeneCalls(['v','j'], action='list'))
        print(self.tr_rec.getGeneCalls(['j','v'], action='list'))
        print(self.tr_rec.getFamilyCalls(['v','d','j'], action='list'))
        print(self.tr_rec.getFamilyCalls(['v','j'], action='list'))
        print(self.tr_rec.getFamilyCalls(['j','v'], action='list'))

        print('JUNK>')
        print(self.bad_rec.getAlleleCalls(['v','d','j'], action='first'))
        print(self.bad_rec.getAlleleCalls(['v','j'], action='first'))
        print(self.bad_rec.getAlleleCalls(['j','v'], action='first'))
        print(self.bad_rec.getGeneCalls(['v','d','j'], action='first'))
        print(self.bad_rec.getGeneCalls(['v','j'], action='first'))
        print(self.bad_rec.getGeneCalls(['j','v'], action='first'))
        print(self.bad_rec.getFamilyCalls(['v','d','j'], action='first'))
        print(self.bad_rec.getFamilyCalls(['v','j'], action='first'))
        print(self.bad_rec.getFamilyCalls(['j','v'], action='first'))

        print(self.bad_rec.getAlleleCalls(['v','d','j'], action='set'))
        print(self.bad_rec.getAlleleCalls(['v','j'], action='set'))
        print(self.bad_rec.getAlleleCalls(['j','v'], action='set'))
        print(self.bad_rec.getGeneCalls(['v','d','j'], action='set'))
        print(self.bad_rec.getGeneCalls(['v','j'], action='set'))
        print(self.bad_rec.getGeneCalls(['j','v'], action='set'))
        print(self.bad_rec.getFamilyCalls(['v','d','j'], action='set'))
        print(self.bad_rec.getFamilyCalls(['v','j'], action='set'))
        print(self.bad_rec.getFamilyCalls(['j','v'], action='set'))

        print(self.bad_rec.getAlleleCalls(['v','d','j'], action='list'))
        print(self.bad_rec.getAlleleCalls(['v','j'], action='list'))
        print(self.bad_rec.getAlleleCalls(['j','v'], action='list'))
        print(self.bad_rec.getGeneCalls(['v','d','j'], action='list'))
        print(self.bad_rec.getGeneCalls(['v','j'], action='list'))
        print(self.bad_rec.getGeneCalls(['j','v'], action='list'))
        print(self.bad_rec.getFamilyCalls(['v','d','j'], action='list'))
        print(self.bad_rec.getFamilyCalls(['v','j'], action='list'))
        print(self.bad_rec.getFamilyCalls(['j','v'], action='list'))

        self.fail('TODO')

    #@unittest.skip("-> decodeCIGAR() skipped\n")
    def test_decodeCIGAR(self):
        for cigar, truth in zip(self.cigar_string, self.cigar_decoded):
            result = decodeCIGAR(cigar)
            print(result)
            self.assertListEqual(truth, result)

    #@unittest.skip("-> decodeBTOP() skipped\n")
    def test_decodeBTOP(self):
        for btop, truth in zip(self.btop_string, self.btop_decoded_full):
            result = decodeBTOP(btop)
            print('FULL> ', result)
            self.assertListEqual(truth, result)

    #@unittest.skip("-> encodeCIGAR() skipped\n")
    def test_encodeCIGAR(self):
        for align, truth in zip(self.cigar_decoded, self.cigar_string):
            result = encodeCIGAR(align)
            print(result)
            self.assertEqual(truth, result)

    #@unittest.skip("-> padAlignment() skipped\n")
    def test_padAlignment(self):
        cigar = self.cigar_decoded[0]
        for (s, r), truth in zip(self.cigar_pos, self.cigar_pad_1):
            #print('POS>', s, r)
            result = padAlignment(cigar, s, r)
            print('PAD>', '(%i, %i) =' % (s, r), result)
            self.assertEqual(truth, result)

        cigar = self.cigar_decoded[1]
        for (s, r), truth in zip(self.cigar_pos, self.cigar_pad_2):
            #print('POS>', s, r)
            result = padAlignment(cigar, s, r)
            print('PAD>', '(%i, %i) =' % (s, r), result)
            self.assertEqual(truth, result)

    #@unittest.skip("-> alignmentPositions() skipped\n")
    def test_alignmentPositions(self):
        for align, truth in zip(self.cigar_decoded, self.cigar_len):
            result = alignmentPositions(align)
            print('POS>', result)
            self.assertDictEqual(truth, result)


if __name__ == '__main__':
    unittest.main()