"""
Sequence manipulation and annotation functions
"""
# Info
__author__ = 'Jason Anthony Vander Heiden, Namita Gupta'
from changeo import __version__, __date__

# Imports
import re
import sys
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# Ig and TCR Regular expressions
allele_regex = re.compile(r'((IG[HLK]|TR[ABGD])([VDJ][A-Z0-9]+[-/\w]*[-\*][\.\w]+))')
gene_regex = re.compile(r'((IG[HLK]|TR[ABGD])([VDJ][A-Z0-9]+[-/\w]*))')
family_regex = re.compile(r'((IG[HLK]|TR[ABGD])([VDJ][A-Z0-9]+))')

v_allele_regex = re.compile(r'((IG[HLK]|TR[ABGD])V[A-Z0-9]+[-/\w]*[-\*][\.\w]+)')
d_allele_regex = re.compile(r'((IG[HLK]|TR[ABGD])D[A-Z0-9]+[-/\w]*[-\*][\.\w]+)')
j_allele_regex = re.compile(r'((IG[HLK]|TR[ABGD])J[A-Z0-9]+[-/\w]*[-\*][\.\w]+)')

#allele_regex = re.compile(r'(IG[HLK][VDJ]\d+[-/\w]*[-\*][\.\w]+)')
#gene_regex = re.compile(r'(IG[HLK][VDJ]\d+[-/\w]*)')
#family_regex = re.compile(r'(IG[HLK][VDJ]\d+)')


# TODO:  might be better to just use the lower case column name as the member variable name. can use getattr and setattr.
class IgRecord:
    """
    A class defining a V(D)J germline sequence alignment
    """
    # Mapping of member variables to column names
    _key_map = {'id': 'SEQUENCE_ID',
                'v_call': 'V_CALL',
                'v_call_geno': 'V_CALL_GENOTYPED',
                'd_call': 'D_CALL',
                'j_call': 'J_CALL',
                'seq_input': 'SEQUENCE_INPUT',
                'seq_vdj': 'SEQUENCE_VDJ',
                'seq_imgt': 'SEQUENCE_IMGT',
                'junction': 'JUNCTION',
                'functional': 'FUNCTIONAL', 
                'in_frame': 'IN_FRAME', 
                'stop': 'STOP', 
                'mutated_invariant': 'MUTATED_INVARIANT', 
                'indels': 'INDELS',
                'v_seq_start': 'V_SEQ_START',
                'v_seq_length': 'V_SEQ_LENGTH',
                'v_germ_start_vdj': 'V_GERM_START_VDJ',
                'v_germ_length_vdj': 'V_GERM_LENGTH_VDJ',
                'v_germ_start_imgt': 'V_GERM_START_IMGT',
                'v_germ_length_imgt': 'V_GERM_LENGTH_IMGT',
                'np1_length': 'NP1_LENGTH',
                'd_seq_start': 'D_SEQ_START',
                'd_seq_length': 'D_SEQ_LENGTH',
                'd_germ_start': 'D_GERM_START',
                'd_germ_length': 'D_GERM_LENGTH',
                'np2_length': 'NP2_LENGTH',
                'j_seq_start': 'J_SEQ_START',
                'j_seq_length': 'J_SEQ_LENGTH',
                'j_germ_start': 'J_GERM_START',
                'j_germ_length': 'J_GERM_LENGTH',
                'junction_length': 'JUNCTION_LENGTH',
                'v_score': 'V_SCORE',
                'v_identity': 'V_IDENTITY',
                'v_evalue': 'V_EVALUE',
                'v_btop': 'V_BTOP',
                'j_score': 'J_SCORE',
                'j_identity': 'J_IDENTITY',
                'j_evalue': 'J_EVALUE',
                'j_btop': 'J_BTOP',
                'hmm_score': 'HMM_SCORE',
                'fwr1': 'FWR1_IMGT',
                'fwr2': 'FWR2_IMGT',
                'fwr3': 'FWR3_IMGT',
                'fwr4': 'FWR4_IMGT',
                'cdr1': 'CDR1_IMGT',
                'cdr2': 'CDR2_IMGT',
                'cdr3': 'CDR3_IMGT',
                'germline': 'GERMLINE',
                'germline_d_mask': 'GERMLINE_D_MASK',
                'n1_length': 'N1_LENGTH',
                'n2_length': 'N2_LENGTH',
                'p3v_length': 'P3V_LENGTH',
                'p5d_length': 'P5D_LENGTH',
                'p3d_length': 'P3D_LENGTH',
                'p5j_length': 'P5J_LENGTH',
                'd_frame': 'D_FRAME',
                'cdr3_igblast_nt': 'CDR3_IGBLAST_NT',
                'cdr3_igblast_aa': 'CDR3_IGBLAST_AA'}

    # Mapping of column names to member variables
    _field_map = {v: k for k, v in _key_map.items()}

    # Mapping of member variables to parsing functions
    _parse_map = {'id': '_identity',
                  'v_call': '_identity',
                  'v_call_geno': '_identity',
                  'd_call': '_identity',
                  'j_call': '_identity',
                  'seq_input': '_sequence',
                  'seq_vdj': '_sequence',
                  'seq_imgt': '_sequence',
                  'junction': '_sequence',
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
                  'np1_length': '_integer',
                  'd_seq_start': '_integer',
                  'd_seq_length': '_integer',
                  'd_germ_start': '_integer',
                  'd_germ_length': '_integer',
                  'np2_length': '_integer',
                  'j_seq_start': '_integer',
                  'j_seq_length': '_integer',
                  'j_germ_start': '_integer',
                  'j_germ_length': '_integer',
                  'junction_length': '_integer',
                  'v_score': '_float',
                  'v_identity': '_float',
                  'v_evalue': '_float',
                  'v_btop': '_identity',
                  'j_score': '_float',
                  'j_identity': '_float',
                  'j_evalue': '_float',
                  'j_btop': '_identity',
                  'hmm_score': '_float',
                  'fwr1': '_sequence',
                  'fwr2': '_sequence',
                  'fwr3': '_sequence',
                  'fwr4': '_sequence',
                  'cdr1': '_sequence',
                  'cdr2': '_sequence',
                  'cdr3': '_sequence',
                  'germline': '_sequence',
                  'germline_d_mask': '_sequence',
                  'n1_length': '_integer',
                  'n2_length': '_integer',
                  'p3v_length': '_integer',
                  'p5d_length': '_integer',
                  'p3d_length': '_integer',
                  'p5j_length': '_integer',
                  'd_frame': '_integer',
                  'cdr3_igblast_nt': '_sequence',
                  'cdr3_igblast_aa': '_sequence'}

    _logical_parse = {'F':False, 'T':True, 'TRUE':True, 'FALSE':False,
                      'NA':None, 'None':None, '':None}
    _logical_deparse = {False:'F', True:'T', None:''}

    # Private methods
    @staticmethod    
    def _identity(v, deparse=False):
        return v

    @staticmethod
    def _logical(v, deparse=False):
        if not deparse:
            try:  return IgRecord._logical_parse[v]
            except:  return None
        else:
            try:  return IgRecord._logical_deparse[v]
            except:  return ''

    @staticmethod
    def _integer(v, deparse=False):
        if not deparse:
            try:  return int(v)
            except:  return ''
        else:
            try:  return str(v)
            except:  return ''

    @staticmethod
    def _float(v, deparse=False):
        if not deparse:
            try:  return float(v)
            except:  return ''
        else:
            try:  return str(v)
            except:  return ''

    @staticmethod
    def _sequence(v, deparse=False):
        if not deparse:
            try:
                if v in ['NA','None']:
                    return ''
                else:
                    return Seq(v, IUPAC.ambiguous_dna).upper()
            except:  return ''
        else:
            try:
                if v in ['NA','None']:
                    return ''
                else:
                    return str(v)
            except:  return ''

    def __init__(self, row, genotyped=True):
        """
        Initializer

        Arguments:
          row : Dictionary of field/value data
          genotyped : If True assign v_call from genotyped field

        Returns:
          IgRecord
        """
        required_keys = ('id',)
        optional_keys = (x for x in IgRecord._parse_map if x not in required_keys)
        
        # Not ideal. Will place V_CALL_GENOTYPED in annotations
        if not genotyped and 'v_call_geno' in optional_keys:
            del optional_keys['v_call_geno']
            
        try:
            for k in required_keys:
                f = getattr(IgRecord, IgRecord._parse_map[k])
                setattr(self, k, f(row.pop(IgRecord._key_map[k])))
        except:
            sys.exit('ERROR:  Input must contain valid %s values' \
                     % ','.join([IgRecord._key_map[k] for k in required_keys]))

        # Defined optional logical values
        for k in optional_keys:
            f = getattr(IgRecord, IgRecord._parse_map[k])
            setattr(self, k, f(row.pop(IgRecord._key_map[k], None)))
            
        # Add remaining elements as annotations dictionary
        self.annotations = row
    
    def getField(self, field):
        """
        Get a field value by column name and return it as a string

        Arguments:
          field : Column name

        Returns:
          str : Value in the field as a string
        """
        if field in IgRecord._field_map:
            v = getattr(self, IgRecord._field_map[field])
        elif field in self.annotations:
            v = self.annotations[field]
        else:
            return None

        if isinstance(v, str):
            return v
        else:
            return str(v)

    def getSeqField(self, field):
        """
        Get a field value converted to a Seq object by column name

        Arguments:
          field : Column name

        Returns:
          Seq : Value in the field as a Seq object
        """
        if field in IgRecord._field_map:
            v = getattr(self, IgRecord._field_map[field])
        elif field in self.annotations:
            v = self.annotations[field]
        else:
            return None

        if isinstance(v, Seq):
            return v
        elif isinstance(v, str):
            return Seq(v, IUPAC.ambiguous_dna)
        else:
            return None

    # Returns: dictionary of the namespace
    def toDict(self):
        d = {}
        n = self.__dict__
        for k, v in n.items():
            if k == 'annotations':
                d.update({i.upper():j for i, j in n['annotations'].items()})
            else:
                f = getattr(IgRecord, IgRecord._parse_map[k])
                d[IgRecord._key_map[k]] = f(v, deparse=True)
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
        vdj = {'v':self.getVGene(action),
               'd':self.getDGene(action),
               'j':self.getJGene(action)}
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
        vdj = {'v':self.getVFamily(action),
               'd':self.getDFamily(action),
               'j':self.getJFamily(action)}
        return [vdj[k] for k in calls]

    # TODO: this can't distinguish empty value ("") from missing field (no column)
    def getVAllele(self, action='first'):
        """
        V-region allele getter

        Arguments:
          actions : One of ('first','set')

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        x = self.v_call_geno if self.v_call_geno is not None else self.v_call
        return parseAllele(x, allele_regex, action)

    def getDAllele(self, action='first'):
        """
        D-region allele getter

        Arguments:
          actions : One of ('first','set')

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        return parseAllele(self.d_call, allele_regex, action)

    def getJAllele(self, action='first'):
        """
        J-region allele getter

        Arguments:
          actions : One of ('first','set')

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
          actions : One of ('first','set')

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        return parseAllele(self.d_call, gene_regex, action)

    def getJGene(self, action='first'):
        """
        J-region gene getter

        Arguments:
          actions : One of ('first','set')

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        return parseAllele(self.j_call, gene_regex, action)
    
    def getVFamily(self, action='first'):
        """
        V-region family getter

        Arguments:
          actions : One of ('first','set')

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        return parseAllele(self.v_call, family_regex, action)

    def getDFamily(self, action='first'):
        """
        D-region family getter

        Arguments:
          actions : One of ('first','set')

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        return parseAllele(self.d_call, family_regex, action)

    def getJFamily(self, action='first'):
        """
        J-region family getter

        Arguments:
          actions : One of ('first','set')

        Returns:
          str : String of the allele when action is 'first';
          tuple : Tuple of allele calls for 'set' or 'list' actions.
        """
        return parseAllele(self.j_call, family_regex, action)


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
