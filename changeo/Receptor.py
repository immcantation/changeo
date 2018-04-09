"""
Sequence manipulation and annotation functions
"""

# Info
__author__ = 'Jason Anthony Vander Heiden, Namita Gupta, Scott Christley'
from changeo import __version__, __date__

# Imports
import re
import sys
from collections import OrderedDict
from itertools import chain
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# Ig and TCR Regular expressions
allele_regex = re.compile(r'((IG[HLK]|TR[ABGD])([VDJ][A-Z0-9]+[-/\w]*[-\*][\.\w]+))')
gene_regex = re.compile(r'((IG[HLK]|TR[ABGD])([VDJ][A-Z0-9]+[-/\w]*))')
family_regex = re.compile(r'((IG[HLK]|TR[ABGD])([VDJ][A-Z0-9]+))')

v_allele_regex = re.compile(r'((IG[HLK]|TR[ABGD])V[A-Z0-9]+[-/\w]*[-\*][\.\w]+)')
d_allele_regex = re.compile(r'((IG[HLK]|TR[ABGD])D[A-Z0-9]+[-/\w]*[-\*][\.\w]+)')
j_allele_regex = re.compile(r'((IG[HLK]|TR[ABGD])J[A-Z0-9]+[-/\w]*[-\*][\.\w]+)')

allele_number_regex = re.compile(r'(?<=\*)([\.\w]+)')
c_gene_regex = re.compile(r'((IG[HLK]|TR[ABGD])([DMAGEC][P0-9]?[A-Z]?))')


class AIRRSchema:
    """
    AIRR format to Receptor mappings
    """
    # Core fields
    core = OrderedDict([('sequence_id', 'sequence_id'),
                        ('sequence', 'sequence_input'),
                        ('sequence_alignment', 'sequence_imgt'),
                        ('germline_alignment', 'germline_imgt'),
                        ('rev_comp', 'rev_comp'),
                        ('productive', 'functional'),
                        ('stop_codon', 'stop'),
                        ('vj_in_frame', 'in_frame'),
                        ('v_call', 'v_call'),
                        ('d_call', 'd_call'),
                        ('j_call', 'j_call'),
                        ('c_call', 'c_call'),
                        ('junction', 'junction'),
                        ('junction_length', 'junction_length'),
                        ('junction_aa', 'junction_aa'),
                        ('np1_length', 'np1_length'),
                        ('np2_length', 'np2_length'),
                        ('v_sequence_start', 'v_seq_start'),
                        ('v_sequence_end', 'v_seq_end'),
                        ('v_germline_start', 'v_germ_start_imgt'),
                        ('v_germline_end', 'v_germ_end_imgt'),
                        ('d_sequence_start', 'd_seq_start'),
                        ('d_sequence_end', 'd_seq_end'),
                        ('d_germline_start', 'd_germ_start'),
                        ('d_germline_end', 'd_germ_end'),
                        ('j_sequence_start', 'j_seq_start'),
                        ('j_sequence_end', 'j_seq_end'),
                        ('j_germline_start', 'j_germ_start'),
                        ('j_germline_end', 'j_germ_end')])
    core_fields = list(core.keys())

    # Alignment scoring fields
    imgt_score = OrderedDict([('v_score', 'v_score'),
                              ('v_identity', 'v_identity'),
                              ('d_score', 'd_score'),
                              ('d_identity', 'd_identity'),
                              ('j_score', 'j_score'),
                              ('j_identity', 'j_identity')])
    imgt_score_fields = list(imgt_score.keys())

    igblast_score = OrderedDict([('v_score', 'v_score'),
                                 ('v_identity', 'v_identity'),
                                 ('v_evalue', 'v_evalue'),
                                 ('v_cigar', 'v_cigar'),
                                 ('d_score', 'd_score'),
                                 ('d_identity', 'd_identity'),
                                 ('d_evalue', 'd_evalue'),
                                 ('d_cigar', 'd_cigar'),
                                 ('j_score', 'j_score'),
                                 ('j_identity', 'j_identity'),
                                 ('j_evalue', 'j_evalue'),
                                 ('j_cigar', 'j_cigar')])
    igblast_score_fields = list(igblast_score.keys())

    ihmm_score = OrderedDict([('vdj_score', 'hmm_score')])
    ihmm_score_fields = list(ihmm_score.keys())

    # FWR andd CDR fields
    region = OrderedDict([('cdr1', 'cdr1_imgt'),
                          ('cdr2', 'cdr2_imgt'),
                          ('cdr3', 'cdr3_imgt'),
                          ('fwr1', 'fwr1_imgt'),
                          ('fwr2', 'fwr2_imgt'),
                          ('fwr3', 'fwr3_imgt'),
                          ('fwr4', 'fwr4_imgt')])
    region_fields = list(region.keys())

    # CDR and FWR position fields
    region_position = OrderedDict([('cdr1_start', 'cdr1_start'),
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
                                   ('fwr4_end', 'fwr4_end')])
    region_position_fields = list(region_position.keys())

    # Detailed junction fields
    junction = OrderedDict([('n1_length', 'n1_length'),
                            ('n2_length', 'n2_length'),
                            ('p3v_length', 'p3v_length'),
                            ('p5d_length', 'p5d_length'),
                            ('p3d_length', 'p3d_length'),
                            ('p5j_length', 'p5j_length'),
                            ('d_frame', 'd_frame')])
    junction_fields = list(junction.keys())

    # IgBLAST CDR3 fields
    igblast_cdr3 = OrderedDict([('cdr3_igblast_nt', 'cdr3_igblast_nt'),
                                ('cdr3_igblast_aa', 'cdr3_igblast_aa'),
                                ('cdr3_igblast_start', 'cdr3_igblast_start'),
                                ('cdr3_igblast_end', 'cdr3_igblast_end')])
    igblast_cdr3_fields = list(igblast_cdr3.keys())

    # Grouping and counting fields
    count = OrderedDict([('duplicate_count', 'dupcount'),
                         ('consensus_count', 'conscount'),
                         ('clone_id', 'clone')])
    count_fields = list(count.keys())

    # Mapping of AIRR column names to Receptor attributes
    _airr = OrderedDict(chain(core.items(),
                              igblast_score.items(),
                              imgt_score.items(),
                              ihmm_score.items(),
                              region.items(),
                              junction.items(),
                              igblast_cdr3.items(),
                              count.items()))

    # Mapping of Receptor attributes to AIRR column names
    _receptor = {v: k for k, v in _airr.items()}

    # Set of starting Receptor positional fields
    _start = ['v_seq_start',
              'v_germ_start_imgt',
              'v_germ_start_vdj',
              'd_seq_start',
              'd_germ_start',
              'j_seq_start',
              'j_germ_start',
              'fwr1_start',
              'fwr2_start',
              'fwr3_start',
              'fwr4_start',
              'cdr1_start',
              'cdr2_start',
              'cdr3_start']

    # Positional fields in the form <Receptor end field>: (<Receptor start field>, <Receptor length field>)
    _end = {'v_seq_end': ('v_seq_start', 'v_seq_length'),
            'v_germ_end_imgt': ('v_germ_start_imgt', 'v_germ_length_imgt'),
            'v_germ_end_vdj': ('v_germ_start_vdj', 'v_germ_length_vdj'),
            'd_seq_end': ('d_seq_start', 'd_seq_length'),
            'd_germ_end': ('d_germ_start', 'd_germ_length'),
            'j_seq_end': ('j_seq_start', 'j_seq_length'),
            'j_germ_end': ('j_germ_start', 'j_germ_length'),
            'fwr1_end': ('fwr1_start', 'fwr1_length'),
            'fwr2_end': ('fwr2_start', 'fwr2_length'),
            'fwr3_end': ('fwr3_start', 'fwr3_length'),
            'fwr4_end': ('fwr4_start', 'fwr4_length'),
            'cdr1_end': ('cdr1_start', 'cdr1_length'),
            'cdr2_end': ('cdr2_start', 'cdr2_length'),
            'cdr3_end': ('cdr3_start', 'cdr3_length')}

    # Ordered list of known fields
    @staticmethod
    def fields(igblast_score=False, imgt_score=False, ihmm_score=False,
               region=False, junction=False, igblast_cdr3=False,
               count=False):
        """
        Returns a list of column names

        Arguments:
          igblast_score : if True include IgBLAST alignment scoring fields
          imgt_score : if True include IMGT alignment scoring fields
          ihmm_score : if True include iHMMune-Align alignment scoring fields
          region : if True include CDR and FWR region fields
          junction : if True include detailed junction fields
          igblast_cdr3 : if True include IgBLAST CDR3 assignment fields
          count : if True include count fields

        Returns:
          list : ordered column names
        """
        fields = AIRRSchema.core_fields[:]
        if igblast_score:  fields.extend(AIRRSchema.igblast_score_fields)
        if imgt_score:  fields.extend(AIRRSchema.imgt_score_fields)
        if ihmm_score:  fields.extend(AIRRSchema.ihmm_score_fields)
        if region:  fields.extend(AIRRSchema.region_fields)
        if junction:  fields.extend(AIRRSchema.junction_fields)
        if igblast_cdr3:  fields.extend(AIRRSchema.igblast_cdr3_fields)
        if count:  fields.extend(AIRRSchema.count_fields)

        return fields

    @staticmethod
    def asReceptor(field):
        """
        Returns a Receptor attribute name from an AIRR column name

        Arguments:
          field : AIRR column name
        Returns:
          str : Receptor attribute name
        """
        field = field.lower()
        return AIRRSchema._airr.get(field, field)

    @staticmethod
    def asAIRR(field, strict=False):
        """
        Returns an AIRR column name from a Receptor attribute name

        Arguments:
          field : Receptor attribute name
          strict : if True return None if the field is not in the AIRR definition;
                   otherwise return the input field name if it is undefined.
        Returns:
          str : AIRR column name
        """
        field = field.lower()
        if strict:
            return AIRRSchema._receptor.get(field, None)
        else:
            return AIRRSchema._receptor.get(field, field)

    @staticmethod
    def asChangeo(field):
        """
        Returns a Change-O column name from an AIRR column name

        Arguments:
          field : AIRR column name
        Returns:
          str : Change-O column name
        """
        return ChangeoSchema.asChangeo(AIRRSchema.asReceptor(field))


class ChangeoSchema:
    """
    Change-O to Receptor mappings
    """
    # Core fields
    core = OrderedDict([('SEQUENCE_ID', 'sequence_id'),
                        ('SEQUENCE_INPUT', 'sequence_input'),
                        ('REV_COMP', 'rev_comp'),
                        ('FUNCTIONAL', 'functional'),
                        ('IN_FRAME', 'in_frame'),
                        ('STOP', 'stop'),
                        ('MUTATED_INVARIANT', 'mutated_invariant'),
                        ('INDELS', 'indels'),
                        ('V_CALL', 'v_call'),
                        ('D_CALL', 'd_call'),
                        ('J_CALL', 'j_call'),
                        ('SEQUENCE_VDJ', 'sequence_vdj'),
                        ('SEQUENCE_IMGT', 'sequence_imgt'),
                        ('V_SEQ_START', 'v_seq_start'),
                        ('V_SEQ_LENGTH', 'v_seq_length'),
                        ('V_GERM_START_VDJ', 'v_germ_start_vdj'),
                        ('V_GERM_LENGTH_VDJ', 'v_germ_length_vdj'),
                        ('V_GERM_START_IMGT', 'v_germ_start_imgt'),
                        ('V_GERM_LENGTH_IMGT', 'v_germ_length_imgt'),
                        ('NP1_LENGTH', 'np1_length'),
                        ('D_SEQ_START', 'd_seq_start'),
                        ('D_SEQ_LENGTH', 'd_seq_length'),
                        ('D_GERM_START', 'd_germ_start'),
                        ('D_GERM_LENGTH', 'd_germ_length'),
                        ('NP2_LENGTH', 'np2_length'),
                        ('J_SEQ_START', 'j_seq_start'),
                        ('J_SEQ_LENGTH', 'j_seq_length'),
                        ('J_GERM_START', 'j_germ_start'),
                        ('J_GERM_LENGTH', 'j_germ_length'),
                        ('JUNCTION', 'junction'),
                        ('JUNCTION_LENGTH', 'junction_length')])
    core_fields = list(core.keys())

    # Alignment scoring fields
    imgt_score = OrderedDict([('V_SCORE', 'v_score'),
                              ('V_IDENTITY', 'v_identity'),
                              ('J_SCORE', 'j_score'),
                              ('J_IDENTITY', 'j_identity')])
    imgt_score_fields = list(imgt_score.keys())

    # IgBLAST scoring fields
    igblast_score = OrderedDict([('V_SCORE', 'v_score'),
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
                                 ('J_CIGAR', 'j_cigar')])
    igblast_score_fields = list(igblast_score.keys())

    ihmm_score = OrderedDict([('HMM_SCORE', 'hmm_score')])
    ihmm_score_fields = list(ihmm_score.keys())

    # Define default FWR amd CDR field ordering
    region = OrderedDict([('FWR1_IMGT', 'fwr1_imgt'),
                          ('FWR2_IMGT', 'fwr2_imgt'),
                          ('FWR3_IMGT', 'fwr3_imgt'),
                          ('FWR4_IMGT', 'fwr4_imgt'),
                          ('CDR1_IMGT', 'cdr1_imgt'),
                          ('CDR2_IMGT', 'cdr2_imgt'),
                          ('CDR3_IMGT', 'cdr3_imgt')])
    region_fields = list(region.keys())

    # Define default detailed junction field ordering
    junction = OrderedDict([('N1_LENGTH', 'n1_length'),
                            ('N2_LENGTH', 'n2_length'),
                            ('P3V_LENGTH', 'p3v_length'),
                            ('P5D_LENGTH', 'p5d_length'),
                            ('P3D_LENGTH', 'p3d_length'),
                            ('P5J_LENGTH', 'p5j_length'),
                            ('D_FRAME', 'd_frame')])
    junction_fields = list(junction.keys())

    # IgBLAST CDR3 field ordering
    igblast_cdr3 = OrderedDict([('CDR3_IGBLAST_NT', 'cdr3_igblast_nt'),
                                ('CDR3_IGBLAST_AA', 'cdr3_igblast_aa'),
                                ('CDR3_IGBLAST_START', 'cdr3_igblast_start'),
                                ('CDR3_IGBLAST_END', 'cdr3_igblast_end')])
    igblast_cdr3_fields = list(igblast_cdr3.keys())

    # Grouping and counting fields
    count = OrderedDict([('CONSCOUNT', 'conscount'),
                         ('DUPCOUNT', 'dupcount'),
                         ('CLONE', 'clone')])
    count_fields = list(count.keys())

    # Mapping of Change-O column names to Receptor attributes
    _changeo = OrderedDict(chain(core.items(),
                                 igblast_score.items(),
                                 imgt_score.items(),
                                 ihmm_score.items(),
                                 region.items(),
                                 junction.items(),
                                 igblast_cdr3.items(),
                                 count.items()))

    # Mapping of Receptor attributes to Change-O column names
    _receptor = {v: k for k, v in _changeo.items()}

    # Ordered list of known fields
    @staticmethod
    def fields(igblast_score=False, imgt_score=False, ihmm_score=False,
               region=False, junction=False, igblast_cdr3=False, count=False):
        """
        Returns a list of column names

        Arguments:
          igblast_score : if True include IgBLAST alignment scoring fields
          imgt_score : if True include IMGT alignment scoring fields
          ihmm_score : if True include iHMMune-Align alignment scoring fields
          region : if True include CDR and FWR region fields
          junction : if True include detailed junction fields
          igblast_cdr3 : if True include IgBLAST CDR3 assignment fields
          count : if True include count fields

        Returns:
          list : ordered column names
        """
        fields = ChangeoSchema.core_fields[:]
        if igblast_score:  fields.extend(ChangeoSchema.igblast_score_fields)
        if imgt_score:  fields.extend(ChangeoSchema.imgt_score_fields)
        if ihmm_score:  fields.extend(ChangeoSchema.ihmm_score_fields)
        if region:  fields.extend(ChangeoSchema.region_fields)
        if junction:  fields.extend(ChangeoSchema.junction_fields)
        if igblast_cdr3:  fields.extend(ChangeoSchema.igblast_cdr3_fields)
        if count:  fields.extend(ChangeoSchema.count_fields)

        return fields

    @staticmethod
    def asReceptor(field):
        """
        Returns a Receptor attribute name from a Change-O column name

        Arguments:
          field : Change-O column name
        Returns:
          str : Receptor attribute name
        """
        return ChangeoSchema._changeo.get(field, field.lower())

    @staticmethod
    def asChangeo(field):
        """
        Returns a Change-O column name from a Receptor attribute name

        Arguments:
          field : Receptor attribute name
        Returns:
          str : Change-O column name
        """
        return ChangeoSchema._receptor.get(field, field.upper())

    @staticmethod
    def asAIRR(field):
        """
        Returns an AIRR column name from a Change-O column name

        Arguments:
          field : Change-O column name
        Returns:
          str : AIRR column name
        """
        return AIRRSchema.asAIRR(ChangeoSchema.asReceptor(field))


class Receptor:
    """
    A class defining a V(D)J sequence and its annotations

    Attributes:
      sequence_id (str) : unique sequence identifier.
      v_call (str) : V allele assignment(s).
      d_call (str) : D allele assignment(s).
      j_call (str) : J allele assignment(s).
      c_call (str) : C region assignment.
      sequence_input (Bio.Seq.Seq) : input nucleotide sequence.
      sequence_imgt	(Bio.Seq.Seq) : IMGT-gapped V(D)J nucleotide sequence.
      sequence_vdj (Bio.Seq.Seq) : ungapped V(D)J nucleotide sequence.
      junction (Bio.Seq.Seq) : ungapped junction region nucletide sequence.
      junction_aa (Bio.Seq.Seq) : ungapped junction region amino acid sequence.
      junction_length (int) : length of the junction in nucleotides.

      functional (bool) : whether sample V(D)J sequence is predicted to be functional.
      rev_comp (bool) : whether the alignment is relative to the reverse compliment of the input sequence.
      in_frame (bool) : whether junction region is in-frame.
      stop (bool) : whether a stop codon is present in the V(D)J sequence.
      mutated_invariant (bool) : whether the conserved amino acids are mutated in the V(D)J sequence.
      indels (bool) : whether the V(D)J nucleotide sequence contains insertions and/or deletions.

      v_seq_start (int) : position of the first V nucleotide in the input sequence.
      v_seq_length (int) : number of V nucleotides in the input sequence.
      v_germ_start_imgt (int) : position of the first V nucleotide in IMGT-gapped V germline sequence alignment.
      v_germ_length_imgt (int) : length of the IMGT numbered germline V alignment.
      v_germ_start_vdj (int) : position of the first nucleotide in ungapped V germline sequence alignment.
      v_germ_length_vdj (int) : length of the ungapped germline V alignment.
      np1_start (int) : position of the first untemplated nucleotide between the V and D segments in the input sequence.
      np1_length (int) : number of untemplated nucleotides between the V and D segments.
      d_seq_start (int) : position of the first D nucleotide in the input sequence.
      d_seq_length (int) : number of D nucleotides in the input sequence.
      d_germ_start (int) : position of the first nucleotide in D germline sequence alignment.
      d_germ_length (int) : length of the germline D alignment.
      np2_start (int) : position of the first untemplated nucleotide between the D and J segments in the input sequence.
      np2_length (int) : number of untemplated nucleotides between the D and J segments.
      j_seq_start (int) : position of the first J nucleotide in the input sequence.
      j_seq_length (int) : number of J nucleotides in the input sequence.
      j_germ_start (int) : position of the first nucleotide in J germline sequence alignment.
      j_germ_length (int) : length of the germline J alignment.

      fwr1_imgt (Bio.Seq.Seq) : IMGT-gapped FWR1 nucleotide sequence.
      fwr2_imgt (Bio.Seq.Seq) : IMGT-gapped FWR2 nucleotide sequence.
      fwr3_imgt (Bio.Seq.Seq) : IMGT-gapped FWR3 nucleotide sequence.
      fwr4_imgt (Bio.Seq.Seq) : IMGT-gapped FWR4 nucleotide sequence.
      cdr1_imgt (Bio.Seq.Seq) : IMGT-gapped CDR1 nucleotide sequence.
      cdr2_imgt (Bio.Seq.Seq) : IMGT-gapped CDR2 nucleotide sequence.
      cdr3_imgt (Bio.Seq.Seq) : IMGT-gapped CDR3 nucleotide sequence.
      cdr3_igblast_nt (Bio.Seq.Seq) : CDR3 nucleotide sequence assigned by IgBLAST.
      cdr3_igblast_aa (Bio.Seq.Seq) : CDR3 amino acid sequence assigned by IgBLAST.

      n1_length (int) : M nucleotides 5' of the D segment.
      n2_length (int) : nucleotides 3' of the D segment.
      p3v_length (int) : palindromic nucleotides 3' of the V segment.
      p5d_length (int) : palindromic nucleotides 5' of the D segment.
      p3d_length (int) : palindromic nucleotides 3' of the D segment.
      p5j_length (int) : palindromic nucleotides 5' of the J segment.
      d_frame (int) : D segment reading frame.

      v_score	(float) : alignment score for the V.
      v_identity (float) : alignment identity for the V.
      v_evalue (float) : E-value for the alignment of the V.
      v_btop (str) : BTOP for the alignment of the V.
      v_cigar (str) : CIGAR for the alignment of the V.
      d_score	(float) : alignment score for the D.
      d_identity	(float) : alignment identity for the D.
      d_evalue (float) : E-value for the alignment of the D.
      d_btop (str) : BTOP for the alignment of the D.
      D_cigar (str) : CIGAR for the alignment of the D.
      j_score	(float) : alignment score for the J.
      j_identity	(float) : alignment identity for the J.
      j_evalue (float) : E-value for the alignment of the J.
      j_btop (str) : BTOP for the alignment of the J.
      j_cigar (str) : CIGAR for the alignment of the J.
      hmm_score (float) : alignment score for the V(D)J.

      germline_vdj (Bio.Seq.Seq) : full ungapped germline V(D)J nucleotide sequence.
      germline_vdj_d_mask (Bio.Seq.Seq) : ungapped germline V(D)J nucleotides sequence with Ns masking the NP1-D-NP2 regions.
      germline_imgt (Bio.Seq.Seq) : full IMGT-gapped germline V(D)J nucleotide sequence.
      germline_imgt_d_mask (Bio.Seq.Seq) : IMGT-gapped germline V(D)J nucleotide sequence with ns masking the NP1-D-NP2 regions.

      conscount (int) : number of reads contributing to the UMI consensus sequence.
      dupcount (int) : copy number of the sequence.
      clone (str): clonal cluster identifier.

      annotations (dict) : dictionary containing all unknown fields.
    """
    # Mapping of member variables to parsing functions
    _parsers = {'sequence_id': '_identity',
                'v_call': '_identity',
                'd_call': '_identity',
                'j_call': '_identity',
                'c_call': '_identity',
                'sequence_input': '_nucleotide',
                'sequence_imgt': '_nucleotide',
                'sequence_vdj': '_nucleotide',
                'junction': '_nucleotide',
                'junction_aa': '_aminoacid',
                'junction_length': '_integer',
                'rev_comp': '_logical',
                'functional': '_logical',
                'in_frame': '_logical',
                'stop': '_logical',
                'mutated_invariant': '_logical',
                'indels': '_logical',
                'v_seq_start': '_integer',
                'v_seq_length': '_integer',
                'v_germ_start_imgt': '_integer',
                'v_germ_length_imgt': '_integer',
                'v_germ_start_vdj': '_integer',
                'v_germ_length_vdj': '_integer',
                'np1_start': '_integer',
                'np1_length': '_integer',
                'd_seq_start': '_integer',
                'd_seq_length': '_integer',
                'd_germ_start': '_integer',
                'd_germ_length': '_integer',
                'np2_start': '_integer',
                'np2_length': '_integer',
                'j_seq_start': '_integer',
                'j_seq_length': '_integer',
                'j_germ_start': '_integer',
                'j_germ_length': '_integer',
                'v_score': '_float',
                'v_identity': '_float',
                'v_evalue': '_float',
                'v_btop': '_identity',
                'v_cigar': '_identity',
                'd_score': '_float',
                'd_identity': '_float',
                'd_evalue': '_float',
                'd_btop': '_identity',
                'd_cigar': '_identity',
                'j_score': '_float',
                'j_identity': '_float',
                'j_evalue': '_float',
                'j_btop': '_identity',
                'j_cigar': '_identity',
                'hmm_score': '_float',
                'fwr1_imgt': '_nucleotide',
                'fwr2_imgt': '_nucleotide',
                'fwr3_imgt': '_nucleotide',
                'fwr4_imgt': '_nucleotide',
                'cdr1_imgt': '_nucleotide',
                'cdr2_imgt': '_nucleotide',
                'cdr3_imgt': '_nucleotide',
                'germline_imgt': '_nucleotide',
                'germline_imgt_d_mask': '_nucleotide',
                'germline_vdj': '_nucleotide',
                'germline_vdj_d_mask': '_nucleotide',
                'n1_length': '_integer',
                'n2_length': '_integer',
                'p3v_length': '_integer',
                'p5d_length': '_integer',
                'p3d_length': '_integer',
                'p5j_length': '_integer',
                'd_frame': '_integer',
                'cdr3_igblast_nt': '_nucleotide',
                'cdr3_igblast_aa': '_aminoacid',
                'conscount': '_integer',
                'dupcount': '_integer',
                'clone': '_identity'}

    # Mapping of derived properties to parsing functions
    _derived = {'v_seq_end': '_integer',
                'v_germ_end_vdj': '_integer',
                'v_germ_end_imgt': '_integer',
                'j_seq_end': '_integer',
                'j_germ_end': '_integer',
                'd_seq_end': '_integer',
                'd_germ_end': '_integer'}

    @staticmethod
    def _identity(v, deparse=False):
        return v

    # Logical type conversion
    @staticmethod
    def _logical(v, deparse=False):
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
    def _integer(v, deparse=False):
        if not deparse:
            try:  return int(v)
            except:  return None
        else:
            return '' if v is None else str(v)

    # Float type conversion
    @staticmethod
    def _float(v, deparse=False):
        if not deparse:
            try:  return float(v)
            except:  return None
        else:
            return '' if v is None else str(v)

    # Nucleotide sequence type conversion
    @staticmethod
    def _nucleotide(v, deparse=False):
        if not deparse:
            try:  return '' if v in ('NA', 'None') else Seq(v, IUPAC.ambiguous_dna).upper()
            except:  return ''
        else:
            return '' if v in ('NA', 'None') else str(v)

    # Sequence type conversion
    @staticmethod
    def _aminoacid(v, deparse=False):
        if not deparse:
            try:  return '' if v in ('NA', 'None') else Seq(v, IUPAC.extended_protein).upper()
            except:  return ''
        else:
            return '' if v in ('NA', 'None') else str(v)

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
        required_keys = ('sequence_id',)
        optional_keys = (x for x in Receptor._parsers if x not in required_keys)

        # Parse required fields
        try:
            for k in required_keys:
                f = getattr(Receptor, Receptor._parsers[k])
                setattr(self, k, f(data.pop(k)))
        except:
            sys.exit('ERROR:  Input must contain valid %s values' % ','.join(required_keys))

        # Parse optional known fields
        for k in optional_keys:
            f = getattr(Receptor, Receptor._parsers[k])
            setattr(self, k, f(data.pop(k, None)))

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
        attributes = {k.lower(): v for k, v in data.items() if k.lower() in Receptor._parsers}
        annotations = {k.lower(): v for k, v in data.items() if k.lower() not in attributes}

        # Update attributes
        for k, v in attributes.items():
            if parse:
                f = getattr(Receptor, Receptor._parsers[k])
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
        if field in Receptor._parsers and parse:
            f = getattr(Receptor, Receptor._parsers[field])
            setattr(self, field, f(value))
        elif field in Receptor._parsers:
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

        if field in Receptor._parsers:
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
        value = self.getField(field)

        if isinstance(value, Seq):
            return value
        elif isinstance(value, str):
            return Seq(value, IUPAC.ambiguous_dna)
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
        field = AIRRSchema.asReceptor(field)

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
        field = ChangeoSchema.asReceptor(field)

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
                f = getattr(Receptor, Receptor._parsers[k])
                d[k] = f(v, deparse=True)
        # Parse properties
        for k in Receptor._derived:
            f = getattr(Receptor, Receptor._derived[k])
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
        return parseAllele(x, allele_regex, action)

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
        return parseAllele(x, allele_regex, action)

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
        return parseAllele(x, allele_regex, action)

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
        return parseAllele(x, gene_regex, action)

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
        return parseAllele(x, gene_regex, action)

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
        return parseAllele(x, gene_regex, action)

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
        return parseAllele(x, family_regex, action)

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
        return parseAllele(x, family_regex, action)

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
        return parseAllele(x, family_regex, action)

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
        return parseAllele(x, allele_number_regex, action)

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
        return parseAllele(x, allele_number_regex, action)

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
        return parseAllele(x, allele_number_regex, action)

    @property
    def v_seq_end(self):
        """
        position of the last V nucleotide in the input sequence.
        """
        try:  return self.v_seq_start + self.v_seq_length - 1
        except TypeError:  return None

    @property
    def v_germ_end_imgt(self):
        """
        position of the last nucleotide in the IMGT-gapped V germline sequence alignment.
        """
        try:  return self.v_germ_start_imgt + self.v_germ_length_imgt - 1
        except TypeError:  return None

    @property
    def v_germ_end_vdj(self):
        """
        position of the last nucleotide in the ungapped V germline sequence alignment.
        """
        try:  return self.v_germ_start_vdj + self.v_germ_length_vdj - 1
        except TypeError:  return None

    @property
    def d_seq_end(self):
        """
        position of the last D nucleotide in the input sequence.
        """
        try:  return self.d_seq_start + self.d_seq_length - 1
        except TypeError:  return None

    @property
    def d_germ_end(self):
        """
        position of the last nucleotide in the D germline sequence alignment.
        """
        try:  return self.d_germ_start + self.d_germ_length - 1
        except TypeError:  return None

    @property
    def j_seq_end(self):
        """
        position of the last J nucleotide in the input sequence.
        """
        try:  return self.j_seq_start + self.j_seq_length - 1
        except TypeError:  return None

    @property
    def j_germ_end(self):
        """
        position of the last nucleotide in the J germline sequence alignment.
        """
        try:  return self.j_germ_start + self.j_germ_length - 1
        except TypeError:  return None

    @property
    def junction_start(self):
        """
        position of the first junction nucleotide in the input sequence.
        """
        try:
            x = self.v_germ_end_imgt - 310
            return self.v_seq_end - x if x >= 0 else None
        except TypeError:
            return None

    @property
    def junction_end(self):
        """
        position of the last junction nucleotide in the input sequence.
        """
        try:
            gaps = self.junction.count('.')
            return self.junction_start + self.junction_length - gaps - 1
        except TypeError:
            return None


# TODO:  might be cleaner as getAllele(), getGene(), getFamily()
def parseAllele(alleles, regex, action='first'):
    """
    Extract alleles from strings

    Arguments:
      alleles : string with allele calls
      regex : compiled regular expression for allele match
      action : action to perform for multiple alleles;
               one of ('first', 'set', 'list').
    Returns:
      str : String of the allele when action is 'first';
      tuple : Tuple of allele calls for 'set' or 'list' actions.
    """
    try:
        match = [x.group(0) for x in regex.finditer(alleles)]
    except:
        match = None

    if action == 'first':
        return match[0] if match else None
    elif action == 'set':
        return tuple(sorted(set(match))) if match else None
    elif action == 'list':
        return tuple(sorted(match)) if match else None
    else:
        return None


def decodeBTOP(btop):
    """
    Parse a BTOP string into a list of tuples in CIGAR annotation.

    Arguments:
      btop : BTOP string.

    Returns:
      list : tuples of (operation, length) for each operation in the BTOP string using CIGAR annotation.
    """
    # Determine chunk type and length
    def _recode(m):
        if m.isdigit():  return ('=', int(m))
        elif m[0] == '-':  return ('I', len(m) // 2)
        elif m[1] == '-':  return ('D', len(m) // 2)
        else:  return ('X', len(m) // 2)

    # Split BTOP string into sections
    btop_split = re.sub(r'(\d+|[-A-Z]{2})', r'\1;', btop)
    # Parse each chunk of encoding
    matches = re.finditer(r'(\d+)|([A-Z]{2};)+|(-[A-Z];)+|([A-Z]-;)+', btop_split)

    return [_recode(m.group().replace(';', '')) for m in matches]


def decodeCIGAR(cigar):
    """
    Parse a CIGAR string into a list of tuples.

    Arguments:
      cigar : CIGAR string.

    Returns:
      list : tuples of (operation, length) for each operation in the CIGAR string.
    """
    matches = re.findall(r'(\d+)([A-Z])', cigar)

    return [(m[1], int(m[0])) for m in matches]


def encodeCIGAR(alignment):
    """
    Encodes a list of tuple with alignment information into a CIGAR string.

    Arguments:
      tuple : tuples of (type, length) for each alignment operation.

    Returns:
      str : CIGAR string.
    """
    return ''.join(['%i%s' % (x, s) for s, x in alignment])


def padAlignment(alignment, q_start, r_start):
    """
    Pads the start of an alignment based on query and reference positions.

    Arguments:
      alignment : tuples of (operation, length) for each alignment operation.
      q_start : query (input) start position (0-based)
      r_start : reference (subject) start position (0-based)

    Returns:
      list : updated list of tuples of (operation, length) for the alignment.
    """
    # Copy list to avoid weirdness
    result = alignment[:]

    # Add query deletions
    if result [0][0] == 'S':
        result[0] = ('S', result[0][1] + q_start)
    elif q_start > 0:
        result.insert(0, ('S', q_start))

    # Add reference padding if present
    if result[0][0] == 'N':
        result[0] = ('N', result[0][1] + r_start)
    elif result [0][0] == 'S' and result[1][0] == 'N':
        result[1] = ('N', result[1][1] + r_start)
    elif result[0][0] == 'S' and r_start > 0:
        result.insert(1, ('N', r_start))
    elif r_start > 0:
        result.insert(0, ('N', r_start))

    return result


def alignmentPositions(alignment):
    """
    Extracts start position and length from an alignment

    Arguments:
      alignment : tuples of (operation, length) for each alignment operation.

    Returns:
      dict : query (q) and reference (r) start and length information with keys
             {q_start, q_length, r_start, r_length}.
    """
    # Return object
    result = {'q_start': 0,
              'q_length': 0,
              'r_start': 0,
              'r_length': 0}

    # Query start
    if alignment[0][0] == 'S':
        result['q_start'] = alignment[0][1]

    # Reference start
    if alignment[0][0] == 'N':
        result['r_start'] = alignment[0][1]
    elif alignment[0][0] == 'S' and alignment[1][0] == 'N':
        result['r_start'] = alignment[1][1]

    # Reference length
    for x, i in alignment:
        if x in ('M', '=', 'X'):
            result['r_length'] += i
            result['q_length'] += i
        elif x == 'D':
            result['r_length'] += i
        elif x == 'I':
            result['q_length'] += i

    return result
