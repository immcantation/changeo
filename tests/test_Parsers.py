"""
Unit tests for MakeDb
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from changeo import __version__, __date__

# Imports
import os
import sys
import time
import unittest

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')

# Import script
sys.path.append(os.path.join(test_path, os.pardir, 'bin'))
import MakeDb
from changeo.IO import extractIMGT, readRepo
from changeo.Parsers import ChangeoReader, IgBLASTReader, IHMMuneReader, IMGTReader, \
                            decodeBTOP, decodeCIGAR, encodeCIGAR, padAlignment, alignmentPositions


class Test_MakeDb(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)

        # Germline files
        self.ig_repo_dir = '/usr/local/share/germlines/imgt/human/vdj'
        # Read files
        self.ig_read_file = os.path.join(data_path, 'reads_ig.fasta')
        # IMGT output
        self.ig_imgt_file = os.path.join(data_path, 'imgt_ig.txz')
        # IgBLAST output
        self.ig_igblast_file = os.path.join(data_path, 'igblast1.7_ig.fmt7')
        # iHMMune-Align output
        self.ig_ihmmune_file = os.path.join(data_path, 'ihmmune_ig.csv')
        # Change-O files
        self.ig_db_file = os.path.join(data_path, 'imgt_ig_db-pass.tsv')

        # CIGAR strings
        self.cigar_string = ['30M1I69M3D', '16P4H30M1I69M3D']
        self.cigar_decoded = [[('M', 30), ('I', 1), ('M', 69), ('D', 3)],
                              [('P', 16), ('H', 4), ('M', 30), ('I', 1), ('M', 69), ('D', 3)]]
        self.cigar_len = [{'q_start': 0, 'q_length': 30 + 1 + 69 + 0, 'r_start': 0, 'r_length': 30 + 0 + 69 + 3},
                          {'q_start': 4, 'q_length': 30 + 1 + 69 + 0, 'r_start': 16, 'r_length': 30 + 0 + 69 + 3}]

        self.cigar_pos = [(0, 0), (5, 0), (0, 3), (5, 3)]
        self.cigar_pad_1 = [[('M', 30), ('I', 1), ('M', 69), ('D', 3)],
                            [('H', 5), ('M', 30), ('I', 1), ('M', 69), ('D', 3)],
                            [('P', 3), ('M', 30), ('I', 1), ('M', 69), ('D', 3)],
                            [('P', 3), ('H', 5), ('M', 30), ('I', 1), ('M', 69), ('D', 3)]]
        self.cigar_pad_2 = [[('P', 16), ('H', 4), ('M', 30), ('I', 1), ('M', 69), ('D', 3)],
                            [('P', 16), ('H', 9), ('M', 30), ('I', 1), ('M', 69), ('D', 3)],
                            [('P', 19), ('H', 4), ('M', 30), ('I', 1), ('M', 69), ('D', 3)],
                            [('P', 19), ('H', 9), ('M', 30), ('I', 1), ('M', 69), ('D', 3)]]

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
        print("<- %s() %.3f" % (self._testMethodName, t))

    @unittest.skip("-> ChangeoReader() skipped\n")
    def test_ChangeoReader(self):
        # Parse
        with open(self.ig_db_file, 'r') as f:
            result = ChangeoReader(f, receptor=False)
            for x in result: print(x)

        with open(self.ig_db_file, 'r') as f:
            result = ChangeoReader(f, receptor=True)
            for x in result: print(x.toDict())

        self.fail('TODO')

    @unittest.skip("-> IMGTReader() skipped\n")
    def test_IMGTReader(self):
        # Extract IMGT files
        temp_dir, files = extractIMGT(self.ig_imgt_file)

        # Parse
        with open(files['summary'], 'r') as summary, \
                open(files['gapped'], 'r') as gapped, \
                open(files['ntseq'], 'r') as ntseq, \
                open(files['junction'], 'r') as junction:
            result = IMGTReader(summary, gapped, ntseq, junction, receptor=False)
            for x in result: print(x)

        # Remove IMGT temporary directory
        temp_dir.cleanup()

        self.fail('TODO')

    @unittest.skip("-> IgBLASTReader() skipped\n")
    def test_IgBLASTReader(self):
        # Load germlines and sequences
        seq_dict = MakeDb.getSeqDict(self.ig_read_file)
        repo_dict = readRepo([self.ig_repo_dir])

        # Parse
        with open(self.ig_igblast_file, 'r') as f:
            result = IgBLASTReader(f, seq_dict, repo_dict, receptor=False)
            for x in result: print(x)

        self.fail('TODO')

    @unittest.skip("-> IHMMReader() skipped\n")
    def test_IHMMReader(self):
        # Load germlines and sequences
        seq_dict = MakeDb.getSeqDict(self.ig_read_file)
        repo_dict = readRepo([self.ig_repo_dir])

        # Parse
        with open(self.ig_ihmmune_file, 'r') as f:
            result = IHMMuneReader(f, seq_dict, repo_dict, receptor=False)
            for x in result: print(x)

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