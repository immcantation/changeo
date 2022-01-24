"""
Receptor data structure
"""

# Info
__author__ = 'Jason Anthony Vander Heiden, Namita Gupta, Scott Christley'

# Imports
from collections import OrderedDict
from Bio.Seq import Seq
# import yaml
# from pkg_resources import resource_stream

# Presto and changeo imports
from presto.IO import printError, printWarning
from changeo.Gene import getAllele, getGene, getFamily, getAlleleNumber

# class Schema:
#     """
#     Schema for mapping Receptor attributes to column names
#     """
#     def __init__(self, schema):
#         """
#         Initializer
#
#         Arguments:
#           schema (str): name of schema to load.
#
#         Returns:
#           changeo.Receptor.Schema
#         """
#         with resource_stream(__name__, 'data/receptor.yaml') as f:
#             data = yaml.load(f, Loader=yaml.FullLoader)
#             receptor = {v[schema]: k for k, v in data['receptor'].items()}
#             definition = data[schema]
#
#         # Output extension
#         self.out_type = definition['out_type']
#
#         # Field sets
#         self.fields = list(receptor.keys())
#         self.required = definition['standard']
#         self.custom_fields = definition['custom']
#
#         # Mapping of schema column names to Receptor attributes
#         self._schema_map = {k : receptor[k] for k in self.fields}
#         self._receptor_map = {v: k for k, v in self._schema_map.items()}
#
#     def toReceptor(self, field):
#         """
#         Returns a Receptor attribute name from an Schema column name
#
#         Arguments:
#           field (str): schema column name.
#         Returns:
#           str: Receptor attribute name.
#         """
#         return self._schema_map.get(field, field)
#
#     def fromReceptor(self, field):
#         """
#         Returns a schema column name from a Receptor attribute name
#
#         Arguments:
#           field (str): Receptor attribute name.
#
#         Returns:
#           str: schema column name.
#         """
#         return self._receptor_map.get(field, field)
#
#
# AIRRSchema = Schema('airr')
# ChangeoSchema = Schema('changeo')


class AIRRSchema:
    """
    AIRR format to Receptor mappings
    """
    # Default file extension
    out_type = 'tsv'

    # Core fields
    required = ['sequence_id',
                'sequence',
                'sequence_alignment',
                'germline_alignment',
                'rev_comp',
                'productive',
                'stop_codon',
                'vj_in_frame',
                'locus',
                'v_call',
                'd_call',
                'j_call',
                'c_call',
                'junction',
                'junction_length',
                'junction_aa',
                'np1_length',
                'np2_length',
                'v_sequence_start',
                'v_sequence_end',
                'v_germline_start',
                'v_germline_end',
                'd_sequence_start',
                'd_sequence_end',
                'd_germline_start',
                'd_germline_end',
                'j_sequence_start',
                'j_sequence_end',
                'j_germline_start',
                'j_germline_end']

    # Mapping of AIRR column names to Receptor attributes
    _schema_map = OrderedDict([('sequence_id', 'sequence_id'),
                               ('sequence', 'sequence_input'),
                               ('sequence_alignment', 'sequence_imgt'),
                               ('germline_alignment', 'germline_imgt'),
                               ('sequence_aa', 'sequence_aa_input'),
                               ('sequence_aa_alignment', 'sequence_aa_imgt'),
                               ('germline_aa_alignment', 'germline_aa_imgt'),
                               ('rev_comp', 'rev_comp'),
                               ('productive', 'functional'),
                               ('stop_codon', 'stop'),
                               ('vj_in_frame', 'in_frame'),
                               ('v_frameshift', 'v_frameshift'),
                               ('locus', 'locus'),
                               ('v_call', 'v_call'),
                               ('d_call', 'd_call'),
                               ('j_call', 'j_call'),
                               ('c_call', 'c_call'),
                               ('junction', 'junction'),
                               ('junction_start', 'junction_start'),
                               ('junction_end', 'junction_end'),
                               ('junction_length', 'junction_length'),
                               ('junction_aa', 'junction_aa'),
                               ('junction_aa_length', 'junction_aa_length'),
                               ('np1_length', 'np1_length'),
                               ('np2_length', 'np2_length'),
                               ('np1_aa_length', 'np1_aa_length'),
                               ('np2_aa_length', 'np2_aa_length'),
                               ('v_sequence_start', 'v_seq_start'),
                               ('v_sequence_end', 'v_seq_end'),
                               ('v_sequence_length', 'v_seq_length'),
                               ('v_germline_start', 'v_germ_start_imgt'),
                               ('v_germline_end', 'v_germ_end_imgt'),
                               ('v_germline_length', 'v_germ_length_imgt'),
                               ('v_sequence_aa_start', 'v_seq_aa_start'),
                               ('v_sequence_aa_end', 'v_seq_aa_end'),
                               ('v_sequence_aa_length', 'v_seq_aa_length'),
                               ('v_germline_aa_start', 'v_germ_aa_start_imgt'),
                               ('v_germline_aa_end', 'v_germ_aa_end_imgt'),
                               ('v_germline_aa_length', 'v_germ_aa_length_imgt'),
                               ('d_sequence_start', 'd_seq_start'),
                               ('d_sequence_end', 'd_seq_end'),
                               ('d_sequence_length', 'd_seq_length'),
                               ('d_germline_start', 'd_germ_start'),
                               ('d_germline_end', 'd_germ_end'),
                               ('d_germline_length', 'd_germ_length'),
                               ('d_sequence_aa_start', 'd_seq_aa_start'),
                               ('d_sequence_aa_end', 'd_seq_aa_end'),
                               ('d_sequence_aa_length', 'd_seq_aa_length'),
                               ('d_germline_aa_start', 'd_germ_aa_start'),
                               ('d_germline_aa_end', 'd_germ_aa_end'),
                               ('d_germline_aa_length', 'd_germ_aa_length'),
                               ('j_sequence_start', 'j_seq_start'),
                               ('j_sequence_end', 'j_seq_end'),
                               ('j_sequence_length', 'j_seq_length'),
                               ('j_germline_start', 'j_germ_start'),
                               ('j_germline_end', 'j_germ_end'),
                               ('j_germline_length', 'j_germ_length'),
                               ('j_sequence_aa_start', 'j_seq_aa_start'),
                               ('j_sequence_aa_end', 'j_seq_aa_end'),
                               ('j_sequence_aa_length', 'j_seq_aa_length'),
                               ('j_germline_aa_start', 'j_germ_aa_start'),
                               ('j_germline_aa_end', 'j_germ_aa_end'),
                               ('j_germline_aa_length', 'j_germ_aa_length'),
                               ('germline_alignment_d_mask', 'germline_imgt_d_mask'),
                               ('v_score', 'v_score'),
                               ('v_identity', 'v_identity'),
                               ('v_support', 'v_evalue'),
                               ('v_cigar', 'v_cigar'),
                               ('d_score', 'd_score'),
                               ('d_identity', 'd_identity'),
                               ('d_support', 'd_evalue'),
                               ('d_cigar', 'd_cigar'),
                               ('j_score', 'j_score'),
                               ('j_identity', 'j_identity'),
                               ('j_support', 'j_evalue'),
                               ('j_cigar', 'j_cigar'),
                               ('vdj_score', 'vdj_score'),
                               ('cdr1', 'cdr1_imgt'),
                               ('cdr2', 'cdr2_imgt'),
                               ('cdr3', 'cdr3_imgt'),
                               ('fwr1', 'fwr1_imgt'),
                               ('fwr2', 'fwr2_imgt'),
                               ('fwr3', 'fwr3_imgt'),
                               ('fwr4', 'fwr4_imgt'),
                               ('cdr1_aa', 'cdr1_aa_imgt'),
                               ('cdr2_aa', 'cdr2_aa_imgt'),
                               ('cdr3_aa', 'cdr3_aa_imgt'),
                               ('fwr1_aa', 'fwr1_aa_imgt'),
                               ('fwr2_aa', 'fwr2_aa_imgt'),
                               ('fwr3_aa', 'fwr3_aa_imgt'),
                               ('fwr4_aa', 'fwr4_aa_imgt'),
                               ('cdr1_start', 'cdr1_start'),
                               ('cdr1_end', 'cdr1_end'),
                               ('cdr2_start', 'cdr2_start'),
                               ('cdr2_end', 'cdr2_end'),
                               ('cdr3_start', 'cdr3_start'),
                               ('cdr3_end', 'cdr3_end'),
                               ('fwr1_start', 'fwr1_start'),
                               ('fwr1_end', 'fwr1_end'),
                               ('fwr2_start', 'fwr2_start'),
                               ('fwr2_end', 'fwr2_end'),
                               ('fwr3_start', 'fwr3_start'),
                               ('fwr3_end', 'fwr3_end'),
                               ('fwr4_start', 'fwr4_start'),
                               ('fwr4_end', 'fwr4_end'),
                               ('n1_length', 'n1_length'),
                               ('n2_length', 'n2_length'),
                               ('p3v_length', 'p3v_length'),
                               ('p5d_length', 'p5d_length'),
                               ('p3d_length', 'p3d_length'),
                               ('p5j_length', 'p5j_length'),
                               ('d_frame', 'd_frame'),
                               ('cdr3_igblast', 'cdr3_igblast'),
                               ('cdr3_igblast_aa', 'cdr3_igblast_aa'),
                               ('duplicate_count', 'dupcount'),
                               ('consensus_count', 'conscount'),
                               ('umi_count', 'umicount'),
                               ('clone_id', 'clone'),
                               ('cell_id', 'cell')])

    # Mapping of Receptor attributes to AIRR column names
    _receptor_map = {v: k for k, v in _schema_map.items()}

    # All fields
    fields = list(_schema_map.keys())

    @staticmethod
    def toReceptor(field):
        """
        Returns a Receptor attribute name from an AIRR column name

        Arguments:
          field : AIRR column name.
        Returns:
          str: Receptor attribute name.
        """
        field = field.lower()
        return AIRRSchema._schema_map.get(field, field)

    @staticmethod
    def fromReceptor(field):
        """
        Returns an AIRR column name from a Receptor attribute name

        Arguments:
          field : Receptor attribute name.

        Returns:
          str: AIRR column name.
        """
        field = field.lower()
        return AIRRSchema._receptor_map.get(field, field)


class AIRRSchemaAA(AIRRSchema):
    """
    AIRR format to Receptor amino acid mappings
    """
    # Core fields
    required = ['sequence_id',
                'sequence',
                'sequence_alignment',
                'germline_alignment',
                'sequence_aa',
                'sequence_aa_alignment',
                'germline_aa_alignment',
                'rev_comp',
                'productive',
                'stop_codon',
                'locus',
                'v_call',
                'd_call',
                'j_call',
                'c_call',
                'junction',
                'junction_length',
                'junction_aa',
                'v_sequence_aa_start',
                'v_sequence_aa_end',
                'v_germline_aa_start',
                'v_germline_aa_end']


class ChangeoSchema:
    """
    Change-O to Receptor mappings
    """
    # Default file extension
    out_type = 'tab'

    # Standard fields
    required = ['SEQUENCE_ID',
                'SEQUENCE_INPUT',
                'FUNCTIONAL',
                'IN_FRAME',
                'STOP',
                'MUTATED_INVARIANT',
                'INDELS',
                'LOCUS',
                'V_CALL',
                'D_CALL',
                'J_CALL',
                'SEQUENCE_VDJ',
                'SEQUENCE_IMGT',
                'V_SEQ_START',
                'V_SEQ_LENGTH',
                'V_GERM_START_VDJ',
                'V_GERM_LENGTH_VDJ',
                'V_GERM_START_IMGT',
                'V_GERM_LENGTH_IMGT',
                'NP1_LENGTH',
                'D_SEQ_START',
                'D_SEQ_LENGTH',
                'D_GERM_START',
                'D_GERM_LENGTH',
                'NP2_LENGTH',
                'J_SEQ_START',
                'J_SEQ_LENGTH',
                'J_GERM_START',
                'J_GERM_LENGTH',
                'JUNCTION',
                'JUNCTION_LENGTH',
                'GERMLINE_IMGT']

    # Mapping of Change-O column names to Receptor attributes
    _schema_map = OrderedDict([('SEQUENCE_ID', 'sequence_id'),
                               ('SEQUENCE_INPUT', 'sequence_input'),
                               ('SEQUENCE_AA_INPUT', 'sequence_aa_input'),
                               ('FUNCTIONAL', 'functional'),
                               ('IN_FRAME', 'in_frame'),
                               ('STOP', 'stop'),
                               ('MUTATED_INVARIANT', 'mutated_invariant'),
                               ('INDELS', 'indels'),
                               ('V_FRAMESHIFT', 'v_frameshift'),
                               ('LOCUS', 'locus'),
                               ('V_CALL', 'v_call'),
                               ('D_CALL', 'd_call'),
                               ('J_CALL', 'j_call'),
                               ('C_CALL', 'c_call'),
                               ('SEQUENCE_VDJ', 'sequence_vdj'),
                               ('SEQUENCE_IMGT', 'sequence_imgt'),
                               ('SEQUENCE_AA_VDJ', 'sequence_aa_vdj'),
                               ('SEQUENCE_AA_IMGT', 'sequence_aa_imgt'),
                               ('V_SEQ_START', 'v_seq_start'),
                               ('V_SEQ_LENGTH', 'v_seq_length'),
                               ('V_GERM_START_VDJ', 'v_germ_start_vdj'),
                               ('V_GERM_LENGTH_VDJ', 'v_germ_length_vdj'),
                               ('V_GERM_START_IMGT', 'v_germ_start_imgt'),
                               ('V_GERM_LENGTH_IMGT', 'v_germ_length_imgt'),
                               ('V_SEQ_AA_START', 'v_seq_aa_start'),
                               ('V_SEQ_AA_LENGTH', 'v_seq_aa_length'),
                               ('V_GERM_AA_START_VDJ', 'v_germ_aa_start_vdj'),
                               ('V_GERM_AA_LENGTH_VDJ', 'v_germ_aa_length_vdj'),
                               ('V_GERM_AA_START_IMGT', 'v_germ_aa_start_imgt'),
                               ('V_GERM_AA_LENGTH_IMGT', 'v_germ_aa_length_imgt'),
                               ('NP1_LENGTH', 'np1_length'),
                               ('NP1_AA_LENGTH', 'np1_aa_length'),
                               ('D_SEQ_START', 'd_seq_start'),
                               ('D_SEQ_LENGTH', 'd_seq_length'),
                               ('D_GERM_START', 'd_germ_start'),
                               ('D_GERM_LENGTH', 'd_germ_length'),
                               ('D_SEQ_AA_START', 'd_seq_aa_start'),
                               ('D_SEQ_AA_LENGTH', 'd_seq_aa_length'),
                               ('D_GERM_AA_START', 'd_germ_aa_start'),
                               ('D_GERM_AA_LENGTH', 'd_germ_aa_length'),
                               ('NP2_LENGTH', 'np2_length'),
                               ('NP2_AA_LENGTH', 'np2_aa_length'),
                               ('J_SEQ_START', 'j_seq_start'),
                               ('J_SEQ_LENGTH', 'j_seq_length'),
                               ('J_GERM_START', 'j_germ_start'),
                               ('J_GERM_LENGTH', 'j_germ_length'),
                               ('J_SEQ_AA_START', 'j_seq_aa_start'),
                               ('J_SEQ_AA_LENGTH', 'j_seq_aa_length'),
                               ('J_GERM_AA_START', 'j_germ_aa_start'),
                               ('J_GERM_AA_LENGTH', 'j_germ_aa_length'),
                               ('JUNCTION', 'junction'),
                               ('JUNCTION_LENGTH', 'junction_length'),
                               ('GERMLINE_IMGT', 'germline_imgt'),
                               ('GERMLINE_AA_IMGT', 'germline_aa_imgt'),
                               ('JUNCTION_START', 'junction_start'),
                               ('V_SCORE', 'v_score'),
                               ('V_IDENTITY', 'v_identity'),
                               ('V_EVALUE', 'v_evalue'),
                               ('V_BTOP', 'v_btop'),
                               ('V_CIGAR', 'v_cigar'),
                               ('D_SCORE', 'd_score'),
                               ('D_IDENTITY', 'd_identity'),
                               ('D_EVALUE', 'd_evalue'),
                               ('D_BTOP', 'd_btop'),
                               ('D_CIGAR', 'd_cigar'),
                               ('J_SCORE', 'j_score'),
                               ('J_IDENTITY', 'j_identity'),
                               ('J_EVALUE', 'j_evalue'),
                               ('J_BTOP', 'j_btop'),
                               ('J_CIGAR', 'j_cigar'),
                               ('VDJ_SCORE', 'vdj_score'),
                               ('FWR1_IMGT', 'fwr1_imgt'),
                               ('FWR2_IMGT', 'fwr2_imgt'),
                               ('FWR3_IMGT', 'fwr3_imgt'),
                               ('FWR4_IMGT', 'fwr4_imgt'),
                               ('CDR1_IMGT', 'cdr1_imgt'),
                               ('CDR2_IMGT', 'cdr2_imgt'),
                               ('CDR3_IMGT', 'cdr3_imgt'),
                               ('FWR1_AA_IMGT', 'fwr1_aa_imgt'),
                               ('FWR2_AA_IMGT', 'fwr2_aa_imgt'),
                               ('FWR3_AA_IMGT', 'fwr3_aa_imgt'),
                               ('FWR4_AA_IMGT', 'fwr4_aa_imgt'),
                               ('CDR1_AA_IMGT', 'cdr1_aa_imgt'),
                               ('CDR2_AA_IMGT', 'cdr2_aa_imgt'),
                               ('CDR3_AA_IMGT', 'cdr3_aa_imgt'),
                               ('N1_LENGTH', 'n1_length'),
                               ('N2_LENGTH', 'n2_length'),
                               ('P3V_LENGTH', 'p3v_length'),
                               ('P5D_LENGTH', 'p5d_length'),
                               ('P3D_LENGTH', 'p3d_length'),
                               ('P5J_LENGTH', 'p5j_length'),
                               ('D_FRAME', 'd_frame'),
                               ('CDR3_IGBLAST', 'cdr3_igblast'),
                               ('CDR3_IGBLAST_AA', 'cdr3_igblast_aa'),
                               ('CONSCOUNT', 'conscount'),
                               ('DUPCOUNT', 'dupcount'),
                               ('UMICOUNT', 'umicount'),
                               ('CLONE', 'clone'),
                               ('CELL', 'cell')])

    # Mapping of Receptor attributes to Change-O column names
    _receptor_map = {v: k for k, v in _schema_map.items()}

    # All fields
    fields = list(_schema_map.keys())

    @staticmethod
    def toReceptor(field):
        """
        Returns a Receptor attribute name from a Change-O column name

        Arguments:
          field : Change-O column name.

        Returns:
          str: Receptor attribute name.
        """
        return ChangeoSchema._schema_map.get(field, field.lower())

    @staticmethod
    def fromReceptor(field):
        """
        Returns a Change-O column name from a Receptor attribute name

        Arguments:
          field : Receptor attribute name.

        Returns:
          str: Change-O column name.
        """
        return ChangeoSchema._receptor_map.get(field, field.upper())


class ChangeoSchemaAA(ChangeoSchema):
    """
    Change-O to Receptor amino acid mappings
    """
    # Standard fields
    required = ['SEQUENCE_ID',
                'SEQUENCE_AA_INPUT',
                'STOP',
                'INDELS',
                'LOCUS',
                'V_CALL',
                'SEQUENCE_AA_VDJ',
                'SEQUENCE_AA_IMGT',
                'V_SEQ_AA_START',
                'V_SEQ_AA_LENGTH',
                'V_GERM_AA_START_VDJ',
                'V_GERM_AA_LENGTH_VDJ',
                'V_GERM_AA_START_IMGT',
                'V_GERM_AA_LENGTH_IMGT',
                'GERMLINE_AA_IMGT']


class ReceptorData:
    """
    A class containing type conversion methods for Receptor data attributes

    Attributes:
      sequence_id (str): unique sequence identifier.

      rev_comp (bool): whether the alignment is relative to the reverse compliment of the input sequence.
      functional (bool): whether sample V(D)J sequence is predicted to be functional.
      in_frame (bool): whether junction region is in-frame.
      stop (bool): whether a stop codon is present in the V(D)J sequence.
      mutated_invariant (bool): whether the conserved amino acids are mutated in the V(D)J sequence.
      indels (bool): whether the V(D)J nucleotide sequence contains insertions and/or deletions.
      v_frameshift (bool): whether the V segment contains a frameshift

      sequence_input (Bio.Seq.Seq): input nucleotide sequence.
      sequence_vdj (Bio.Seq.Seq): Aligned V(D)J nucleotide sequence without IMGT-gaps.
      sequence_imgt (Bio.Seq.Seq): IMGT-gapped V(D)J nucleotide sequence.

      sequence_aa_input (Bio.Seq.Seq): input amino acid sequence.
      sequence_aa_vdj (Bio.Seq.Seq): Aligned V(D)J nucleotide sequence without IMGT-gaps.
      sequence_aa_imgt (Bio.Seq.Seq): IMGT-gapped V(D)J amino sequence.

      junction (Bio.Seq.Seq): ungapped junction region nucletide sequence.
      junction_aa (Bio.Seq.Seq): ungapped junction region amino acid sequence.
      junction_start (int): start positions of the junction in the input nucleotide sequence.
      junction_length (int): length of the junction in nucleotides.

      germline_vdj (Bio.Seq.Seq): full ungapped germline V(D)J nucleotide sequence.
      germline_vdj_d_mask (Bio.Seq.Seq): ungapped germline V(D)J nucleotides sequence with Ns masking the NP1-D-NP2 regions.
      germline_imgt (Bio.Seq.Seq): full IMGT-gapped germline V(D)J nucleotide sequence.
      germline_imgt_d_mask (Bio.Seq.Seq): IMGT-gapped germline V(D)J nucleotide sequence with ns masking the NP1-D-NP2 regions.

      germline_aa_vdj (Bio.Seq.Seq): full ungapped germline V(D)J amino acid sequence.
      germline_aa_imgt (Bio.Seq.Seq): full IMGT-gapped germline V(D)J amino acid sequence.

      v_call (str): V allele assignment(s).
      d_call (str): D allele assignment(s).
      j_call (str): J allele assignment(s).
      c_call (str): C region assignment.

      v_seq_start (int): position of the first V nucleotide in the input sequence (1-based).
      v_seq_length (int): number of V nucleotides in the input sequence.
      v_germ_start_imgt (int): position of the first V nucleotide in IMGT-gapped V germline sequence alignment (1-based).
      v_germ_length_imgt (int): length of the IMGT numbered germline V alignment.
      v_germ_start_vdj (int): position of the first nucleotide in ungapped V germline sequence alignment (1-based).
      v_germ_length_vdj (int): length of the ungapped germline V alignment.

      v_seq_aa_start (int): position of the first V amino acid in the amino acid input sequence (1-based).
      v_seq_aa_length (int): number of V amino acid in the amino acid input sequence.
      v_germ_aa_start_imgt (int): position of the first V amino acid in IMGT-gapped V germline amino acid alignment (1-based).
      v_germ_aa_length_imgt (int): length of the IMGT numbered germline V amino acid alignment.
      v_germ_aa_start_vdj (int): position of the first amino acid in ungapped V germline amino acid alignment (1-based).
      v_germ_aa_length_vdj (int): length of the ungapped germline V amino acid alignment.

      np1_start (int): position of the first untemplated nucleotide between the V and D segments in the input sequence (1-based).
      np1_length (int): number of untemplated nucleotides between the V and D segments.

      np1_aa_start (int): position of the first untemplated amino acid between the V and D segments in the input amino acid sequence (1-based).
      np1_aa_length (int): number of untemplated amino acids between the V and D segments.

      d_seq_start (int): position of the first D nucleotide in the input sequence (1-based).
      d_seq_length (int): number of D nucleotides in the input sequence.
      d_germ_start (int): position of the first nucleotide in D germline sequence alignment (1-based).
      d_germ_length (int): length of the germline D alignment.

      d_seq_aa_start (int): position of the first D amino acid in the input amino acidsequence (1-based).
      d_seq_aa_length (int): number of D amino acids in the input amino acid sequence.
      d_germ_aa_start (int): position of the first amino acid in D germline amino acid alignment (1-based).
      d_germ_aa_length (int): length of the germline D amino acid alignment.

      np2_start (int): position of the first untemplated nucleotide between the D and J segments in the input sequence (1-based).
      np2_length (int): number of untemplated nucleotides between the D and J segments.

      np2_aa_start (int): position of the first untemplated amino acid between the D and J segments in the input amino acid sequence (1-based).
      np2_aa_length (int): number of untemplated amino acid between the D and J segments.

      j_seq_start (int): position of the first J nucleotide in the input sequence (1-based).
      j_seq_length (int): number of J nucleotides in the input sequence.
      j_germ_start (int): position of the first nucleotide in J germline sequence alignment (1-based).
      j_germ_length (int): length of the germline J alignment.

      j_seq_aa_start (int): position of the first J amino acid in the input amino acidsequence (1-based).
      j_seq_aa_length (int): number of J amino acid in the input amino acidsequence.
      j_germ_aa_start (int): position of the first amino acid in J germline amino acid alignment (1-based).
      j_germ_aa_length (int): length of the germline J amino acid alignment.

      v_score (float): alignment score for the V.
      v_identity (float): alignment identity for the V.
      v_evalue (float): E-value for the alignment of the V.
      v_btop (str): BTOP for the alignment of the V.
      v_cigar (str): CIGAR for the alignment of the V.

      d_score (float): alignment score for the D.
      d_identity (float): alignment identity for the D.
      d_evalue (float): E-value for the alignment of the D.
      d_btop (str): BTOP for the alignment of the D.
      d_cigar (str): CIGAR for the alignment of the D.

      j_score (float): alignment score for the J.
      j_identity (float): alignment identity for the J.
      j_evalue (float): E-value for the alignment of the J.
      j_btop (str): BTOP for the alignment of the J.
      j_cigar (str): CIGAR for the alignment of the J.

      vdj_score (float): alignment score for the V(D)J.

      fwr1_imgt (Bio.Seq.Seq): IMGT-gapped FWR1 nucleotide sequence.
      fwr2_imgt (Bio.Seq.Seq): IMGT-gapped FWR2 nucleotide sequence.
      fwr3_imgt (Bio.Seq.Seq): IMGT-gapped FWR3 nucleotide sequence.
      fwr4_imgt (Bio.Seq.Seq): IMGT-gapped FWR4 nucleotide sequence.
      cdr1_imgt (Bio.Seq.Seq): IMGT-gapped CDR1 nucleotide sequence.
      cdr2_imgt (Bio.Seq.Seq): IMGT-gapped CDR2 nucleotide sequence.
      cdr3_imgt (Bio.Seq.Seq): IMGT-gapped CDR3 nucleotide sequence.
      cdr3_igblast (Bio.Seq.Seq): CDR3 nucleotide sequence assigned by IgBLAST.

      fwr1_aa_imgt (Bio.Seq.Seq): IMGT-gapped FWR1 amino acid sequence.
      fwr2_aa_imgt (Bio.Seq.Seq): IMGT-gapped FWR2 amino acid sequence.
      fwr3_aa_imgt (Bio.Seq.Seq): IMGT-gapped FWR3 amino acid sequence.
      fwr4_aa_imgt (Bio.Seq.Seq): IMGT-gapped FWR4 amino acid sequence.
      cdr1_aa_imgt (Bio.Seq.Seq): IMGT-gapped CDR1 amino acid sequence.
      cdr2_aa_imgt (Bio.Seq.Seq): IMGT-gapped CDR2 amino acid sequence.
      cdr3_aa_imgt (Bio.Seq.Seq): IMGT-gapped CDR3 amino acid sequence.
      cdr3_igblast_aa (Bio.Seq.Seq): CDR3 amino acid sequence assigned by IgBLAST.

      n1_length (int): M nucleotides 5' of the D segment.
      n2_length (int): nucleotides 3' of the D segment.
      p3v_length (int): palindromic nucleotides 3' of the V segment.
      p5d_length (int): palindromic nucleotides 5' of the D segment.
      p3d_length (int): palindromic nucleotides 3' of the D segment.
      p5j_length (int): palindromic nucleotides 5' of the J segment.
      d_frame (int): D segment reading frame.

      conscount (int): number of reads contributing to the UMI consensus sequence.
      dupcount (int): copy number of the sequence.
      umicount (int): number of UMIs representing the sequence.

      clone (str): clonal cluster identifier.
      cell (str): origin cell identifier.

      annotations (dict): dictionary containing all unknown fields.
    """
    #with resource_stream(__name__, 'data/receptor.yaml') as f:
    #    data = yaml.load(f, Loader=yaml.FullLoader)
    #
    # # Define type parsers
    # parsers = {k: v['type'] for k, v in data['receptor'].items()}
    #
    # # Define coordinate field sets
    # coordinates = {}
    # for k, v in data['receptor'].items():
    #     if 'coordinate' in v:
    #         position = {v['coordinate']['position']: k}
    #         group = coordinates.setdefault(v['coordinate']['group'], {})
    #         group.update(position)
    #
    # # Positional fields sets in the form {start: (length, end)}
    # self.start_fields = {x['start']: (x['length'], x['end']) for x in coordinates.values()}
    #
    # # Positional fields sets in the form {length: (start, end)}
    # self.length_fields = {x['length']: (x['start'], x['end']) for x in coordinates.values()}
    #
    # # Positional fields sets in the form {end: (start, length)}
    # self.end_fields = {x['end']: (x['start'], x['length']) for x in coordinates.values()}

    # Mapping of member variables to parsing functions
    parsers = {'sequence_id': 'identity',
               'rev_comp': 'logical',
               'functional': 'logical',
               'locus': 'identity',
               'in_frame': 'logical',
               'stop': 'logical',
               'mutated_invariant': 'logical',
               'indels': 'logical',
               'v_frameshift': 'logical',
               'sequence_input': 'nucleotide',
               'sequence_imgt': 'nucleotide',
               'sequence_vdj': 'nucleotide',
               'sequence_aa_input': 'aminoacid',
               'sequence_aa_imgt': 'aminoacid',
               'sequence_aa_vdj': 'aminoacid',
               'junction': 'nucleotide',
               'junction_aa': 'aminoacid',
               'junction_start': 'integer',
               'junction_length': 'integer',
               'germline_imgt': 'nucleotide',
               'germline_imgt_d_mask': 'nucleotide',
               'germline_vdj': 'nucleotide',
               'germline_vdj_d_mask': 'nucleotide',
               'germline_aa_imgt': 'aminoacid',
               'germline_aa_vdj': 'aminoacid',
               'v_call': 'identity',
               'd_call': 'identity',
               'j_call': 'identity',
               'c_call': 'identity',
               'v_seq_start': 'integer',
               'v_seq_length': 'integer',
               'v_germ_start_imgt': 'integer',
               'v_germ_length_imgt': 'integer',
               'v_germ_start_vdj': 'integer',
               'v_germ_length_vdj': 'integer',
               'v_seq_aa_start': 'integer',
               'v_seq_aa_length': 'integer',
               'v_germ_aa_start_imgt': 'integer',
               'v_germ_aa_length_imgt': 'integer',
               'v_germ_aa_start_vdj': 'integer',
               'v_germ_aa_length_vdj': 'integer',
               'np1_start': 'integer',
               'np1_length': 'integer',
               'np1_aa_start': 'integer',
               'np1_aa_length': 'integer',
               'd_seq_start': 'integer',
               'd_seq_length': 'integer',
               'd_germ_start': 'integer',
               'd_germ_length': 'integer',
               'd_seq_aa_start': 'integer',
               'd_seq_aa_length': 'integer',
               'd_germ_aa_start': 'integer',
               'd_germ_aa_length': 'integer',
               'np2_start': 'integer',
               'np2_length': 'integer',
               'np2_aa_start': 'integer',
               'np2_aa_length': 'integer',
               'j_seq_start': 'integer',
               'j_seq_length': 'integer',
               'j_germ_start': 'integer',
               'j_germ_length': 'integer',
               'j_seq_aa_start': 'integer',
               'j_seq_aa_length': 'integer',
               'j_germ_aa_start': 'integer',
               'j_germ_aa_length': 'integer',
               'v_score': 'double',
               'v_identity': 'double',
               'v_evalue': 'double',
               'v_btop': 'identity',
               'v_cigar': 'identity',
               'd_score': 'double',
               'd_identity': 'double',
               'd_evalue': 'double',
               'd_btop': 'identity',
               'd_cigar': 'identity',
               'j_score': 'double',
               'j_identity': 'double',
               'j_evalue': 'double',
               'j_btop': 'identity',
               'j_cigar': 'identity',
               'vdj_score': 'double',
               'fwr1_imgt': 'nucleotide',
               'fwr2_imgt': 'nucleotide',
               'fwr3_imgt': 'nucleotide',
               'fwr4_imgt': 'nucleotide',
               'cdr1_imgt': 'nucleotide',
               'cdr2_imgt': 'nucleotide',
               'cdr3_imgt': 'nucleotide',
               'fwr1_aa_imgt': 'aminoacid',
               'fwr2_aa_imgt': 'aminoacid',
               'fwr3_aa_imgt': 'aminoacid',
               'fwr4_aa_imgt': 'aminoacid',
               'cdr1_aa_imgt': 'aminoacid',
               'cdr2_aa_imgt': 'aminoacid',
               'cdr3_aa_imgt': 'aminoacid',
               'n1_length': 'integer',
               'n2_length': 'integer',
               'p3v_length': 'integer',
               'p5d_length': 'integer',
               'p3d_length': 'integer',
               'p5j_length': 'integer',
               'd_frame': 'integer',
               'cdr3_igblast': 'nucleotide',
               'cdr3_igblast_aa': 'aminoacid',
               'conscount': 'integer',
               'dupcount': 'integer',
               'umicount': 'integer',
               'clone': 'identity',
               'cell': 'identity'}

    # Positional fields sets in the form (start, length, end)
    _coordinate_map = [('v_seq_start', 'v_seq_length', 'v_seq_end'),
                       ('v_germ_start_imgt', 'v_germ_length_imgt', 'v_germ_end_imgt'),
                       ('v_germ_start_vdj', 'v_germ_length_vdj', 'v_germ_end_vdj'),
                       ('v_alignment_start', 'v_alignment_length', 'v_alignment_end'),
                       ('v_seq_aa_start', 'v_seq_aa_length', 'v_seq_aa_end'),
                       ('v_germ_aa_start_imgt', 'v_germ_aa_length_imgt', 'v_germ_aa_end_imgt'),
                       ('v_germ_aa_start_vdj', 'v_germ_aa_length_vdj', 'v_germ_aa_end_vdj'),
                       ('v_alignment_aa_start', 'v_alignment_aa_length', 'v_alignment_aa_end'),
                       ('d_seq_start', 'd_seq_length', 'd_seq_end'),
                       ('d_germ_start', 'd_germ_length', 'd_germ_end'),
                       ('d_seq_aa_start', 'd_seq_aa_length', 'd_seq_aa_end'),
                       ('d_germ_aa_start', 'd_germ_aa_length', 'd_germ_aa_end'),
                       ('j_seq_start', 'j_seq_length', 'j_seq_end'),
                       ('j_germ_start', 'j_germ_length', 'j_germ_end'),
                       ('j_seq_aa_start', 'j_seq_aa_length', 'j_seq_aa_end'),
                       ('j_germ_aa_start', 'j_germ_aa_length', 'j_germ_aa_end'),
                       ('junction_start', 'junction_length', 'junction_end'),
                       ('fwr1_start', 'fwr1_length', 'fwr1_end'),
                       ('fwr2_start', 'fwr2_length', 'fwr2_end'),
                       ('fwr3_start', 'fwr3_length', 'fwr3_end'),
                       ('fwr4_start', 'fwr4_length', 'fwr4_end'),
                       ('cdr1_start', 'cdr1_length', 'cdr1_end'),
                       ('cdr2_start', 'cdr2_length', 'cdr2_end'),
                       ('cdr3_start', 'cdr3_length', 'cdr3_end')]

    # Positional fields sets in the form {start: (length, end)}
    start_fields = {x[0]: (x[1], x[2]) for x in _coordinate_map}

    # Positional fields sets in the form {length: (start, end)}
    length_fields = {x[1]: (x[0], x[2]) for x in _coordinate_map}

    # Positional fields sets in the form {end: (start, length)}
    end_fields = {x[2]: (x[0], x[1]) for x in _coordinate_map}

    @staticmethod
    def identity(v, deparse=False):
        return v

    # Logical type conversion
    @staticmethod
    def logical(v, deparse=False):
        parse_map = {True: True, 'T': True, 'TRUE': True,
                     False: False, 'F': False, 'FALSE': False,
                     'NA': None, 'None': None, '': None}
        deparse_map = {False: 'F', True: 'T', None: ''}
        if not deparse:
            try:  return parse_map[v]
            except:  return None
        else:
            try:  return deparse_map[v]
            except:  return ''

    # Integer type conversion
    @staticmethod
    def integer(v, deparse=False):
        if not deparse:
            try:  return int(v)
            except:  return None
        else:
            return '' if v is None else str(v)

    # Float type conversion
    @staticmethod
    def double(v, deparse=False):
        if not deparse:
            try:  return float(v)
            except:  return None
        else:
            return '' if v is None else str(v)

    # Nucleotide sequence type conversion
    @staticmethod
    def nucleotide(v, deparse=False):
        if not deparse:
            try:
                #return '' if v in ('NA', 'None') else Seq(v, IUPAC.ambiguous_dna).upper()
                return '' if v in ('NA', 'None') else v.upper()
            except:
                return ''
        else:
            return '' if v in ('NA', 'None', None) else str(v)

    # Sequence type conversion
    @staticmethod
    def aminoacid(v, deparse=False):
        if not deparse:
            try:
                #return '' if v in ('NA', 'None') else Seq(v, IUPAC.extended_protein).upper()
                return '' if v in ('NA', 'None') else v.upper()
            except:
                return ''
        else:
            return '' if v in ('NA', 'None', None) else str(v)


class Receptor:
    """
    A class defining a V(D)J sequence and its annotations
    """
    # Mapping of derived properties to types
    _derived = {'v_seq_end': 'integer',
                'v_germ_end_vdj': 'integer',
                'v_germ_end_imgt': 'integer',
                'v_seq_aa_end': 'integer',
                'v_germ_aa_end_vdj': 'integer',
                'v_germ_aa_end_imgt': 'integer',
                'd_seq_end': 'integer',
                'd_germ_end': 'integer',
                'd_seq_aa_end': 'integer',
                'd_germ_aa_end': 'integer',
                'j_seq_end': 'integer',
                'j_germ_end': 'integer',
                'j_seq_aa_end': 'integer',
                'j_germ_aa_end': 'integer',
                'junction_end': 'integer'}

    def _junction_start(self):
        """
        Determine the position of the first junction nucleotide in the input sequence
        """
        try:
            x = self.v_germ_end_imgt - 310
            return self.v_seq_end - x if x >= 0 else None
        except TypeError:
            return None

    def __init__(self, data):
        """
        Initializer

        Arguments:
          data : dict of field/value data

        Returns:
          changeo.Receptor.Receptor
        """
        # Convert case of keys
        data = {k.lower(): v for k, v in data.items()}

        # Define known keys
        required_keys = ('sequence_id', )
        optional_keys = (x for x in ReceptorData.parsers if x not in required_keys)

        # Parse required fields
        try:
            for k in required_keys:
                f = getattr(ReceptorData, ReceptorData.parsers[k])
                setattr(self, k, f(data.pop(k)))
        except:
            printError('Input must contain valid %s values.' % ','.join(required_keys))

        # Parse optional known fields
        for k in optional_keys:
            f = getattr(ReceptorData, ReceptorData.parsers[k])
            setattr(self, k, f(data.pop(k, None)))

        # Derive junction_start if not provided
        if not hasattr(self, 'junction_start') or self.junction_start is None:
            setattr(self, 'junction_start', self._junction_start())

        # Add remaining elements as annotations dictionary
        self.annotations = data

    def setDict(self, data, parse=False):
        """
        Adds or updates multiple attributes and annotations

        Arguments:
          data : a dictionary of annotations to add or update.
          parse : if True pass values through string parsing functions for known fields.

        Returns:
          None : updates attribute values and the annotations attribute.
        """
        # Partition data
        attributes = {k.lower(): v for k, v in data.items() if k.lower() in ReceptorData.parsers}
        annotations = {k.lower(): v for k, v in data.items() if k.lower() not in attributes}

        # Update attributes
        for k, v in attributes.items():
            if parse:
                f = getattr(ReceptorData, ReceptorData.parsers[k])
                setattr(self, k, f(v))
            else:
                setattr(self, k, v)

        # Update annotations
        self.annotations.update(annotations)

    def setField(self, field, value, parse=False):
        """
        Set an attribute or annotation value

        Arguments:
          field : attribute name as a string
          value : value to assign
          parse : if True pass values through string parsing functions for known fields.

        Returns:
          None. Updates attribute or annotation.
        """
        field = field.lower()
        if field in ReceptorData.parsers and parse:
            f = getattr(ReceptorData, ReceptorData.parsers[field])
            setattr(self, field, f(value))
        elif field in ReceptorData.parsers:
            setattr(self, field, value)
        else:
            self.annotations[field] = value

    def getField(self, field):
        """
        Get an attribute or annotation value

        Arguments:
          field : attribute name as a string

        Returns:
          Value in the attribute. Returns None if the attribute cannot be found.
        """
        field = field.lower()

        if field in ReceptorData.parsers:
            return getattr(self, field)
        elif field in self.annotations:
            return self.annotations[field]
        else:
            return None

    def getSeq(self, field):
        """
        Get an attribute value converted to a Seq object

        Arguments:
          field : variable name as a string

        Returns:
          Bio.Seq.Seq : Value in the field as a Seq object
        """
        v = self.getField(field)

        if isinstance(v, Seq):
            return v
        elif isinstance(v, str):
            return Seq(v)
        else:
            return None

    def getAIRR(self, field, seq=False):
        """
        Get an attribute from an AIRR field name

        Arguments:
          field : AIRR column name as a string
          seq : if True return the attribute as a Seq object

        Returns:
          Value in the AIRR field. Returns None if the field cannot be found.
        """
        # Map to Receptor attribute
        field = AIRRSchema.toReceptor(field)

        if seq:
            return self.getSeq(field)
        else:
            return self.getField(field)

    def getChangeo(self, field, seq=False):
        """
        Get an attribute from a Change-O field name

        Arguments:
          field : Change-O column name as a string
          seq : if True return the attribute as a Seq object

        Returns:
          Value in the Change-O field. Returns None if the field cannot be found.
        """
        # Map to Receptor attribute
        field = ChangeoSchema.toReceptor(field)

        if seq:
            return self.getSeq(field)
        else:
            return self.getField(field)

    def toDict(self):
        """
        Convert the namespace to a dictionary

        Returns:
          dict : member fields with values converted to appropriate strings
        """
        d = {}
        n = self.__dict__
        # Parse attributes
        for k, v in n.items():
            if k == 'annotations':
                d.update(n['annotations'])
            else:
                f = getattr(ReceptorData, ReceptorData.parsers[k])
                d[k] = f(v, deparse=True)
        # Parse properties
        for k in Receptor._derived:
            f = getattr(ReceptorData, Receptor._derived[k])
            v = getattr(self, k)
            d[k] = f(v, deparse=True)
        return d

    def getAlleleCalls(self, calls, action='first'):
        """
        Get multiple allele calls

        Arguments:
          calls : iterable of calls to get; one or more of ('v','d','j')
          actions : One of ('first','set')

        Returns:
          list : List of requested calls in order
        """
        vdj = {'v': self.getVAllele(action),
               'd': self.getDAllele(action),
               'j': self.getJAllele(action)}

        return [vdj[k] for k in calls]

    def getGeneCalls(self, calls, action='first'):
        """
        Get multiple gene calls

        Arguments:
          calls : iterable of calls to get; one or more of ('v','d','j')
          actions : One of ('first','set')

        Returns:
          list : List of requested calls in order
        """
        vdj = {'v': self.getVGene(action),
               'd': self.getDGene(action),
               'j': self.getJGene(action)}

        return [vdj[k] for k in calls]

    def getFamilyCalls(self, calls, action='first'):
        """
        Get multiple family calls

        Arguments:
          calls : iterable of calls to get; one or more of ('v','d','j')
          actions : One of ('first','set')

        Returns:
          list : List of requested calls in order
        """
        vdj = {'v': self.getVFamily(action),
               'd': self.getDFamily(action),
               'j': self.getJFamily(action)}

        return [vdj[k] for k in calls]

    # TODO: this can't distinguish empty value ("") from missing field (no column)
    def getVAllele(self, action='first', field=None):
        """
        V segment allele getter

        Arguments:
          actions : One of 'first', 'set' or list'
          field : attribute or annotation name containing the V call. Use v_call attribute if None.

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        x = self.v_call if field is None else self.getField(field)
        return getAllele(x, action=action)

    def getDAllele(self, action='first', field=None):
        """
        D segment allele getter

        Arguments:
          actions : One of 'first', 'set' or 'list'
          field : attribute or annotation name containing the D call. Use d_call attribute if None.

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        x = self.d_call if field is None else self.getField(field)
        return getAllele(x, action=action)

    def getJAllele(self, action='first', field=None):
        """
        J segment allele getter

        Arguments:
          actions : One of 'first', 'set' or 'list'
          field : attribute or annotation name containing the J call. Use j_call attribute if None.

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        x = self.j_call if field is None else self.getField(field)
        return getAllele(x, action=action)

    def getVGene(self, action='first', field=None):
        """
        V segment gene getter

        Arguments:
          actions : One of 'first', 'set' or list'
          field : attribute or annotation name containing the V call. Use v_call attribute if None.

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        x = self.v_call if field is None else self.getField(field)
        return getGene(x, action=action)

    def getDGene(self, action='first', field=None):
        """
        D segment gene getter

        Arguments:
          actions : One of 'first', 'set' or list'
          field : attribute or annotation name containing the D call. Use d_call attribute if None.

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        x = self.d_call if field is None else self.getField(field)
        return getGene(x, action=action)

    def getJGene(self, action='first', field=None):
        """
        J segment gene getter

        Arguments:
          actions : One of 'first', 'set' or list'
          field : attribute or annotation name containing the J call. Use j_call attribute if None.

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        x = self.j_call if field is None else self.getField(field)
        return getGene(x, action=action)

    def getVFamily(self, action='first', field=None):
        """
        V segment family getter

        Arguments:
          actions : One of 'first', 'set' or list'
          field : attribute or annotation name containing the V call. Use v_call attribute if None.

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        x = self.v_call if field is None else self.getField(field)
        return getFamily(x, action=action)

    def getDFamily(self, action='first', field=None):
        """
        D segment family getter

        Arguments:
          actions : One of 'first', 'set' or list'
          field : attribute or annotation name containing the D call. Use d_call attribute if None.

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        x = self.d_call if field is None else self.getField(field)
        return getFamily(x, action=action)

    def getJFamily(self, action='first', field=None):
        """
        J segment family getter

        Arguments:
          actions : One of 'first', 'set' or list'
          field : attribute or annotation name containing the J call. Use j_call attribute if None.

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        x = self.j_call if field is None else self.getField(field)
        return getFamily(x, action=action)

    def getAlleleNumbers(self, calls, action='first'):
        """
        Get multiple allele numeric identifiers

        Arguments:
          calls : iterable of calls to get; one or more of ('v','d','j')
          actions : One of ('first','set')

        Returns:
          list : List of requested calls in order
        """
        vdj = {'v': self.getVAlleleNumber(action),
               'd': self.getDAlleleNumber(action),
               'j': self.getJAlleleNumber(action)}

        return [vdj[k] for k in calls]

    def getVAlleleNumber(self, action='first', field=None):
        """
        V segment allele number getter

        Arguments:
          actions : One of 'first', 'set' or list'
          field : attribute or annotation name containing the V call. Use v_call attribute if None.

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele numbers for 'set' or 'list' actions.
        """
        x = self.v_call if field is None else self.getField(field)
        return getAlleleNumber(x, action=action)

    def getDAlleleNumber(self, action='first', field=None):
        """
        D segment allele number getter

        Arguments:
          actions : One of 'first', 'set' or list'
          field : attribute or annotation name containing the D call. Use d_call attribute if None.

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele numbers for 'set' or 'list' actions.
        """
        x = self.d_call if field is None else self.getField(field)
        return getAlleleNumber(x, action=action)

    def getJAlleleNumber(self, action='first', field=None):
        """
        J segment allele number getter

        Arguments:
          actions : One of 'first', 'set' or list'
          field : attribute or annotation name containing the J call. Use j_call attribute if None.

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele numbers for 'set' or 'list' actions.
        """
        x = self.j_call if field is None else self.getField(field)
        return getAlleleNumber(x, action=action)

    @property
    def v_seq_end(self):
        """
        Position of the last V nucleotide in the input sequence
        """
        try:  return self.v_seq_start + self.v_seq_length - 1
        except TypeError:  return None

    @property
    def v_germ_end_imgt(self):
        """
        Position of the last nucleotide in the IMGT-gapped V germline sequence alignment
        """
        try:  return self.v_germ_start_imgt + self.v_germ_length_imgt - 1
        except TypeError:  return None

    @property
    def v_germ_end_vdj(self):
        """
        Position of the last nucleotide in the ungapped V germline sequence alignment
        """
        try:  return self.v_germ_start_vdj + self.v_germ_length_vdj - 1
        except TypeError:  return None

    @property
    def v_seq_aa_end(self):
        """
        Position of the last V nucleotide in the input sequence
        """
        try:  return self.v_seq_aa_start + self.v_seq_aa_length - 1
        except TypeError:  return None

    @property
    def v_germ_aa_end_imgt(self):
        """
        Position of the last nucleotide in the IMGT-gapped V germline sequence alignment
        """
        try:  return self.v_germ_aa_start_imgt + self.v_germ_aa_length_imgt - 1
        except TypeError:  return None

    @property
    def v_germ_aa_end_vdj(self):
        """
        Position of the last nucleotide in the ungapped V germline sequence alignment
        """
        try:  return self.v_germ_aa_start_vdj + self.v_germ_aa_length_vdj - 1
        except TypeError:  return None

    @property
    def d_seq_end(self):
        """
        Position of the last D nucleotide in the input sequence
        """
        try:  return self.d_seq_start + self.d_seq_length - 1
        except TypeError:  return None

    @property
    def d_germ_end(self):
        """
        Position of the last nucleotide in the D germline sequence alignment
        """
        try:  return self.d_germ_start + self.d_germ_length - 1
        except TypeError:  return None

    @property
    def d_seq_aa_end(self):
        """
        Position of the last D amino acid in the input amino acid sequence
        """
        try:  return self.d_seq_aa_start + self.d_seq_aa_length - 1
        except TypeError:  return None

    @property
    def d_germ_aa_end(self):
        """
        Position of the last amino acid in the D germline amino acid alignment
        """
        try:  return self.d_germ_aa_start + self.d_germ_aa_length - 1
        except TypeError:  return None

    @property
    def j_seq_end(self):
        """
        Position of the last J nucleotide in the input sequence
        """
        try:  return self.j_seq_start + self.j_seq_length - 1
        except TypeError:  return None

    @property
    def j_germ_end(self):
        """
        Position of the last nucleotide in the J germline sequence alignment
        """
        try:  return self.j_germ_start + self.j_germ_length - 1
        except TypeError:  return None

    @property
    def j_seq_aa_end(self):
        """
        Position of the last J amino acid in the input amino sequence
        """
        try:  return self.j_seq_aa_start + self.j_seq_aa_length - 1
        except TypeError:  return None

    @property
    def j_germ_aa_end(self):
        """
        Position of the last amino acid in the J germline amino acid alignment
        """
        try:  return self.j_germ_aa_start + self.j_germ_aa_length - 1
        except TypeError:  return None

    @property
    def junction_end(self):
        """
        Position of the last junction nucleotide in the input sequence
        """
        try:
            gaps = self.junction.count('.')
            return self.junction_start + self.junction_length - gaps - 1
        except TypeError:
            return None
