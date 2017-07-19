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

# TODO: Add a method to extract attribute from changeo column name
class Receptor:
    """
    A class defining a V(D)J germline sequence alignment
    """
    # Mapping of member variables to parsing functions
    _parse_map = {'sequence_id': '_identity',
                  'v_call': '_identity',
                  'v_call_genotyped': '_identity',
                  'd_call': '_identity',
                  'j_call': '_identity',
                  'sequence_input': '_sequence',
                  'sequence_vdj': '_sequence',
                  'sequence_imgt': '_sequence',
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
                  'junction_start': '_integer',
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
                  'fwr1_imgt': '_sequence',
                  'fwr2_imgt': '_sequence',
                  'fwr3_imgt': '_sequence',
                  'fwr4_imgt': '_sequence',
                  'cdr1_imgt': '_sequence',
                  'cdr2_imgt': '_sequence',
                  'cdr3_imgt': '_sequence',
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

    # Pass through type conversion
    @staticmethod
    def _identity(v, deparse=False):
        return v

    # Logical type conversion
    @staticmethod
    def _logical(v, deparse=False):
        logical_parse = {'F': False, 'T': True, 'TRUE': True, 'FALSE': False,
                         'NA': None, 'None': None, '': None}
        logical_deparse = {False: 'F', True: 'T', None: ''}
        if not deparse:
            try:
                return logical_parse[v]
            except:
                return None
        else:
            try:
                return logical_deparse[v]
            except:
                return ''

    # Integer type conversion
    @staticmethod
    def _integer(v, deparse=False):
        if not deparse:
            try:
                return int(v)
            except:
                return ''
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
                return ''
        else:
            try:
                return str(v)
            except:
                return ''

    # Sequence type conversion
    @staticmethod
    def _sequence(v, deparse=False):
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

    def __init__(self, data):
        """
        Initializer

        Arguments:
          data : dict of field/value data

        Returns:
          changeo.Receptor.Receptor
        """
        # Convert case
        data = {k.lower(): v for k, v in data.items()}

        # Define known keys
        required_keys = ('sequence_id',)
        optional_keys = (x for x in Receptor._parse_map if x not in required_keys)

        # Parse required fields
        try:
            for k in required_keys:
                f = getattr(Receptor, Receptor._parse_map[k])
                setattr(self, k, f(data.pop(k)))
        except:
            sys.exit('ERROR:  Input must contain valid %s values' % ','.join(required_keys))

        # Parse optional known fields
        for k in optional_keys:
            f = getattr(Receptor, Receptor._parse_map[k])
            setattr(self, k, f(data.pop(k, None)))

        # Add remaining elements as annotations dictionary
        self.annotations = data

    def toSeq(self, field):
        """
        Get a field value converted to a Seq object

        Arguments:
          field : variable name as a string

        Returns:
          Seq : Value in the field as a Seq object
        """
        field = field.lower()

        if field in Receptor._parse_map:
            v = getattr(self, field)
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
                d.update(n['annotations'])
            else:
                f = getattr(Receptor, Receptor._parse_map[k])
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
