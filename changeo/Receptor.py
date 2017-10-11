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
    # TODO: Extra AIRR annotation columns currently not included
    # sample
    # locus
    # constant
    # isotype
    # rev_comp

    # Mapping of AIRR column names to Receptor attributes
    _airr = OrderedDict([('sequence_id', 'sequence_id'),
                         ('sequence', 'sequence_input'),
                         ('functional', 'functional'),
                         ('v_call', 'v_call'),
                         ('d_call', 'd_call'),
                         ('j_call', 'j_call'),
                         ('junction_nt', 'junction'),
                         ('junction_aa', 'junction_aa'),
                         ('v_score', 'v_score'),
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
                         ('j_cigar', 'j_cigar'),
                         ('vdj_score', 'vdj_score'),
                         ('vdj_identity', 'vdj_identity'),
                         ('vdj_evalue', 'vdj_evalue'),
                         ('vdj_cigar', 'vdj_cigar'),
                         ('np1_length', 'np1_length'),
                         ('np2_length', 'np2_length'),
                         ('n1_length', 'n1_length'),
                         ('n2_length', 'n2_length'),
                         ('p3v_length', 'p3v_length'),
                         ('p5d_length', 'p5d_length'),
                         ('p3d_length', 'p3d_length'),
                         ('p5j_length', 'p5j_length'),
                         ('v_start', 'v_seq_start'),
                         ('v_germ_start', 'v_germ_start_imgt'),
                         ('v_end', 'v_seq_end'),
                         ('v_germ_end', 'v_germ_end_imgt'),
                         ('d_start', 'd_seq_start'),
                         ('d_germ_start', 'd_germ_start_imgt'),
                         ('d_end', 'd_seq_end'),
                         ('d_end_start', 'd_germ_end_imgt'),
                         ('j_start', 'j_seq_start'),
                         ('j_germ_start', 'j_germ_start_imgt'),
                         ('j_end', 'j_seq_end'),
                         ('j_germ_end', 'j_germ_end_imgt'),
                         ('fwr1_start', 'fwr1_start'),
                         ('fwr1_end', 'fwr1_end'),
                         ('fwr2_start', 'fwr2_start'),
                         ('fwr2_end', 'fwr2_end'),
                         ('fwr3_start', 'fwr3_start'),
                         ('fwr3_end', 'fwr3_end'),
                         ('fwr4_start', 'fwr4_start'),
                         ('fwr4_end', 'fwr4_end'),
                         ('cdr1_start', 'cdr1_start'),
                         ('cdr1_end', 'cdr1_end'),
                         ('cdr2_start', 'cdr2_start'),
                         ('cdr2_end', 'cdr2_end'),
                         ('cdr3_start', 'cdr3_start'),
                         ('cdr3_end', 'cdr3_end'),
                         ('duplicate_count', 'dupcount'),
                         ('consensus_count', 'conscount')])

    # Mapping of Receptor attributes to Change-O column names
    _receptor = {v: k for k, v in _airr.items()}

    # Ordered list of known fields
    @staticmethod
    def fields():
        """
        Returns a list of column names

        Returns:
          list : ordered column names
        """
        return list(AIRRSchema._airr.keys())

    @staticmethod
    def asReceptor(field):
        """
        Returns a Receptor attribute name from an AIRR column name

        Arguments:
          field : AIRR column name
        Returns:
          str : Receptor attribute name
        """
        return AIRRSchema._airr.get(field.lower(), field.lower())

    @staticmethod
    def asAIRR(field):
        """
        Returns an AIRR column name from a Receptor attribute name

        Arguments:
          field : Receptor attribute name
        Returns:
          str : AIRR column name
        """
        return AIRRSchema._receptor.get(field.lower(), field.lower())


class ChangeoSchema:
    """
    Change-O to Receptor mappings
    """
    # Core fields
    core = OrderedDict([('SEQUENCE_ID', 'sequence_id'),
                        ('SEQUENCE_INPUT', 'sequence_input'),
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
                                ('CDR3_IGBLAST_AA', 'cdr3_igblast_aa')])
    igblast_cdr3_fields = list(igblast_cdr3.keys())

    # Mapping of Change-O column names to Receptor attributes
    _changeo = OrderedDict(chain(core.items(),
                                 igblast_score.items(),
                                 imgt_score.items(),
                                 ihmm_score.items(),
                                 region.items(),
                                 junction.items(),
                                 igblast_cdr3.items()))

    # Mapping of Receptor attributes to Change-O column names
    _receptor = {v: k for k, v in _changeo.items()}

    # Ordered list of known fields
    @staticmethod
    def fields(igblast_score=False, imgt_score=False, ihmm_score=False,
               region=False, junction=False, igblast_cdr3=False):
        """
        Returns a list of column names

        Arguments:
          igblast_score : if True include IgBLAST alignment scoring fields
          imgt_score : if True include IMGT alignment scoring fields
          ihmm_score : if True include iHMMune-Align alignment scoring fields
          region : if True include CDR and FWR region fields
          junction : if True include detailed junction fields
          igblast_cdr3 : if True include IgBLAST CDR3 assignment fields

        Returns:
          list : ordered column names
        """
        f = ChangeoSchema.core_fields[:]
        if igblast_score:  f.extend(ChangeoSchema.igblast_score_fields)
        if imgt_score:  f.extend(ChangeoSchema.imgt_score_fields)
        if ihmm_score:  f.extend(ChangeoSchema.ihmm_score_fields)
        if region:  f.extend(ChangeoSchema.region_fields)
        if junction:  f.extend(ChangeoSchema.junction_fields)
        if igblast_cdr3:  f.extend(ChangeoSchema.igblast_cdr3_fields)

        return f

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


class Receptor:
    """
    A class defining a V(D)J germline sequence alignment
    """
    # Mapping of member variables to parsing functions
    _parsers = {'sequence_id': '_identity',
                'v_call': '_identity',
                'v_call_genotyped': '_identity',
                'd_call': '_identity',
                'j_call': '_identity',
                'sequence_input': '_nucleotide',
                'sequence_vdj': '_nucleotide',
                'sequence_imgt': '_nucleotide',
                'junction': '_nucleotide',
                'junction_aa': '_aminoacid',
                'functional': '_logical',
                'in_frame': '_logical',
                'stop': '_logical',
                'mutated_invariant': '_logical',
                'indels': '_logical',
                'v_seq_start': '_integer',
                'v_seq_length': '_integer',
                'v_germ_start_vdj': '_integer',
                'v_germ_length_vdj': '_integer',
                'v_germ_start_imgt': '_integer',
                'v_germ_length_imgt': '_integer',
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
                'junction_length': '_integer',
                'l_seq_start': '_integer',
                'l_seq_length': '_integer',
                'c_seq_start': '_integer',
                'c_seq_length': '_integer',
                'v_score': '_float',
                'v_identity': '_float',
                'v_evalue': '_float',
                'v_btop': '_identity',
                'j_score': '_float',
                'j_identity': '_float',
                'j_evalue': '_float',
                'j_btop': '_identity',
                'hmm_score': '_float',
                'fwr1_imgt': '_nucleotide',
                'fwr2_imgt': '_nucleotide',
                'fwr3_imgt': '_nucleotide',
                'fwr4_imgt': '_nucleotide',
                'cdr1_imgt': '_nucleotide',
                'cdr2_imgt': '_nucleotide',
                'cdr3_imgt': '_nucleotide',
                'germline': '_nucleotide',
                'germline_d_mask': '_nucleotide',
                'n1_length': '_integer',
                'n2_length': '_integer',
                'p3v_length': '_integer',
                'p5d_length': '_integer',
                'p3d_length': '_integer',
                'p5j_length': '_integer',
                'd_frame': '_integer',
                'cdr3_igblast_nt': '_nucleotide',
                'cdr3_igblast_aa': '_aminoacid'}

    # Pass through type conversion
    @staticmethod
    def _identity(v, deparse=False):
        return v

    # Logical type conversion
    @staticmethod
    def _logical(v, deparse=False):
        parse_map = {'F': False, 'T': True, 'TRUE': True, 'FALSE': False,
                         'NA': None, 'None': None, '': None}
        deparse_map = {False: 'F', True: 'T', None: ''}
        if not deparse:
            try:
                return parse_map[v]
            except:
                return None
        else:
            try:
                return deparse_map[v]
            except:
                return ''

    # Integer type conversion
    @staticmethod
    def _integer(v, deparse=False):
        if not deparse:
            try:
                return int(v)
            except:
                return 0
        else:
            try:
                return str(v)
            except:
                return ''

    # Float type conversion
    @staticmethod
    def _float(v, deparse=False):
        if not deparse:
            try:
                return float(v)
            except:
                return 0.0
        else:
            try:
                return str(v)
            except:
                return ''

    # Nucleotide sequence type conversion
    @staticmethod
    def _nucleotide(v, deparse=False):
        if not deparse:
            try:
                if v in ('NA', 'None'):
                    return ''
                else:
                    return Seq(v, IUPAC.ambiguous_dna).upper()
            except:
                return ''
        else:
            try:
                if v in ('NA', 'None'):
                    return ''
                else:
                    return str(v)
            except:
                return ''

    # Sequence type conversion
    @staticmethod
    def _aminoacid(v, deparse=False):
        if not deparse:
            try:
                if v in ('NA', 'None'):
                    return ''
                else:
                    return Seq(v, IUPAC.extended_protein).upper()
            except:
                return ''
        else:
            try:
                if v in ('NA', 'None'):
                    return ''
                else:
                    return str(v)
            except:
                return ''

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

    def updateAnnotations(self, data):
        """
        Add entries to annotations

        Arguments:
          data : a dictionary of annotations to add

        Returns:
          None : updates the annotations attribute
        """
        # Convert case of keys
        data = {k.lower(): v for k, v in data.items()}
        self.annotations.update(data)

    def getField(self, field):
        """
        Get an attribute value

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
        for k, v in n.items():
            if k == 'annotations':
                #d.update({j.lower(): u for j, u  in n['annotations'].items()})
                d.update(n['annotations'])
            else:
                f = getattr(Receptor, Receptor._parsers[k])
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
    def getVAllele(self, action='first'):
        """
        V-region allele getter

        Arguments:
          actions : One of ('first', 'set')

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        x = self.v_call_genotyped if self.v_call_genotyped is not None else self.v_call
        return parseAllele(x, allele_regex, action)

    def getDAllele(self, action='first'):
        """
        D-region allele getter

        Arguments:
          actions : One of ('first', 'set')

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        return parseAllele(self.d_call, allele_regex, action)

    def getJAllele(self, action='first'):
        """
        J-region allele getter

        Arguments:
          actions : One of ('first', 'set')

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        return parseAllele(self.j_call, allele_regex, action)

    def getVGene(self, action='first'):
        """
        V-region gene getter

        Arguments:
          actions : One of ('first','set')

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        return parseAllele(self.v_call, gene_regex, action)

    def getDGene(self, action='first'):
        """
        D-region gene getter

        Arguments:
          actions : One of ('first', 'set')

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        return parseAllele(self.d_call, gene_regex, action)

    def getJGene(self, action='first'):
        """
        J-region gene getter

        Arguments:
          actions : One of ('first', 'set')

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        return parseAllele(self.j_call, gene_regex, action)

    def getVFamily(self, action='first'):
        """
        V-region family getter

        Arguments:
          actions : One of ('first', 'set')

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        return parseAllele(self.v_call, family_regex, action)

    def getDFamily(self, action='first'):
        """
        D-region family getter

        Arguments:
          actions : One of ('first', 'set')

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        return parseAllele(self.d_call, family_regex, action)

    def getJFamily(self, action='first'):
        """
        J-region family getter

        Arguments:
          actions : One of ('first', 'set')

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        return parseAllele(self.j_call, family_regex, action)

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

    def getVAlleleNumber(self, action='first'):
        """
        V-region allele number getter

        Arguments:
          actions : One of ('first','set')

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele numbers for 'set' or 'list' actions.
        """
        return parseAllele(self.v_call, allele_number_regex, action)

    def getDAlleleNumber(self, action='first'):
        """
        D-region allele number getter

        Arguments:
          actions : One of ('first','set')

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele numbers for 'set' or 'list' actions.
        """
        return parseAllele(self.d_call, allele_number_regex, action)

    def getJAlleleNumber(self, action='first'):
        """
        J-region allele number getter

        Arguments:
          actions : One of ('first','set')

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele numbers for 'set' or 'list' actions.
        """
        return parseAllele(self.j_call, allele_number_regex, action)

    @property
    def v_seq_end(self):
        return self.v_seq_start + self.v_seq_length - 1

    @property
    def v_germ_end_vdj(self):
        return self.v_germ_start_vdj + self.v_germ_length_vdj - 1

    @property
    def v_germ_end_imgt(self):
        return self.v_germ_start_imgt + self.v_germ_length_imgt - 1

    @property
    def d_seq_end(self):
        return self.d_seq_start + self.d_seq_length - 1

    @property
    def d_germ_end(self):
        return self.d_germ_start + self.d_germ_length - 1

    @property
    def j_seq_end(self):
        return self.j_seq_start + self.j_seq_length - 1

    @property
    def j_germ_end(self):
        return self.j_germ_start + self.j_germ_length - 1

    @property
    def junction_start(self):
        x = self.v_germ_end_imgt - 310
        return self.v_seq_end - x if x >= 0 else None

    @property
    def junction_end(self):
        gaps = self.junction.count('.')
        return self.junction_start + self.junction_length - gaps - 1


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
