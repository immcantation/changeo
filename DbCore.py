#!/usr/bin/env python
"""
Core functions shared by Changeo modules
"""

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2013 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.0'
__date__      = '2013.12.1'

# Imports
import csv, os, re, sys
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# Defaults
default_id_field = 'ID'
default_seq_field = 'SEQ'


class IgRecord:
    """
    A class defining a V(D)J germline sequence alignment
    """
    # Private variables
    _key_map = {'id': 'SEQUENCE_ID',
                'v_call': 'V_CALL',
                'd_call': 'D_CALL',
                'j_call': 'J_CALL',
                'seq': 'SEQUENCE', 
                'seq_gap': 'SEQUENCE_GAP',
                'junction': 'JUNCTION',
                'functional': 'FUNCTIONAL', 
                'in_frame': 'IN_FRAME', 
                'stop': 'STOP', 
                'mutated_invariant': 'MUTATED_INVARIANT', 
                'indels': 'INDELS',
                'v_match': 'V_MATCH',
                'v_length': 'V_LENGTH',
                'j_match': 'J_MATCH',
                'j_length': 'J_LENGTH',
                'v_gap_length': 'V_GAP_LENGTH',
                'n1_length': 'N1_LENGTH',
                'd_5_trim': 'D_5_TRIM',
                'd_3_trim': 'D_3_TRIM',
                'n2_length': 'N2_LENGTH',
                'j_5_trim': 'J_5_TRIM',
                'j_gap_length': 'J_GAP_LENGTH',
                'junction_gap_length': 'JUNCTION_GAP_LENGTH'}
    
    _parse_map = {'id': '_identity',
                  'v_call': '_identity',
                  'd_call': '_identity',
                  'j_call': '_identity',
                  'seq': '_sequence', 
                  'seq_gap': '_sequence',
                  'junction': '_sequence',
                  'functional': '_logical', 
                  'in_frame': '_logical', 
                  'stop': '_logical', 
                  'mutated_invariant': '_logical', 
                  'indels': '_logical',
                  'v_match': '_integer',
                  'v_length': '_integer',
                  'j_match': '_integer',
                  'j_length': '_integer',
                  'v_gap_length': '_integer',
                  'n1_length': '_integer',
                  'd_5_trim': '_integer',
                  'd_3_trim': '_integer',
                  'n2_length': '_integer',
                  'j_5_trim': '_integer',
                  'j_gap_length': '_integer',
                  'junction_gap_length': '_integer'}

    _logical_parse = {'F':False, 'T':True}
    _logical_deparse = {False:'F', True:'T'}
    
    # Public variables
    allele_regex = re.compile(r'(IG[HLK][VDJ]\d+[-/\w]*[-\*][\.\w]+)')
    gene_regex = re.compile(r'(IG[HLK][VDJ]\d+[-/\w]*)')
    family_regex = re.compile(r'(IG[HLK][VDJ]\d+)')

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
            except:  return None
        else:
            try:  return str(v)
            except:  return ''
            
    @staticmethod
    def _sequence(v, deparse=False):
        if not deparse:
            try:  return Seq(v, IUPAC.ambiguous_dna)
            except:  return None
        else:
            try:  return str(v)
            except:  return ''
    
    @staticmethod
    def _parseAllele(alleles, regex, action='first'):
        x = regex.findall(alleles)
        if action == 'first':
            return x[0] if x else None
        elif action == 'set':
            return tuple(sorted(set(x))) if x else None
        elif action == 'list':
            return tuple(sorted(x)) if x else None
        else:
            return None

    # Initializer
    def __init__(self, row):
        required_keys = ('id', 'seq', 'v_call', 'd_call', 'j_call')
        optional_keys = (x for x in IgRecord._parse_map if x not in required_keys)
        
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
    
    # Return a dictionary of the namespace
    def toDict(self):
        d = {}
        n = self.__dict__
        for k, v in n.iteritems():
            if k == 'annotations':
                d.update({i.upper():j for i, j in n['annotations'].iteritems()})
            else:
                f = getattr(IgRecord, IgRecord._parse_map[k])
                d[IgRecord._key_map[k]] = f(v, deparse=True)
        return d
    
    # Allele, gene and family getter functions
    def getVAllele(self, action='first'):
        return IgRecord._parseAllele(self.v_call, self.allele_regex, action)

    def getDAllele(self, action='first'):
        return IgRecord._parseAllele(self.d_call, self.allele_regex, action)

    def getJAllele(self, action='first'):
        return IgRecord._parseAllele(self.j_call, self.allele_regex, action)
    
    def getVGene(self, action='first'):
        return IgRecord._parseAllele(self.v_call, self.gene_regex, action)

    def getDGene(self, action='first'):
        return IgRecord._parseAllele(self.d_call, self.gene_regex, action)

    def getJGene(self, action='first'):
        return IgRecord._parseAllele(self.j_call, self.gene_regex, action)
    
    def getVFamily(self, action='first'):
        return IgRecord._parseAllele(self.v_call, self.family_regex, action)

    def getDFamily(self, action='first'):
        return IgRecord._parseAllele(self.d_call, self.family_regex, action)

    def getJFamily(self, action='first'):
        return IgRecord._parseAllele(self.j_call, self.family_regex, action)


def getDbWriter(out_handle, in_file=None, add_fields=None, exclude_fields=None):
    """
    Opens a writer object for an output database file
    
    Arguments: 
    out_handle = the file handle to write to
    in_file = the input filename to determine output fields from;
              if None do not define output fields from input file
    add_fields = a list of fields added to the writer not present in the in_file;
                 if None do not add fields
    exclude_fields = a list of fields in the in_file excluded from the writer;
                     if None do not exclude fields
    
    Returns:
    a writer object
    """
    # Get output field names from input file
    if in_file is not None:
        fields = (readDbFile(in_file, ig=False)).fieldnames
    else:
        fields = []
    # Add extra fields
    if add_fields is not None:
        if not isinstance(add_fields, list):  add_fields = [add_fields]
        fields.extend([f for f in add_fields if f not in fields])
    # Remove unwanted fields
    if exclude_fields is not None:
        if not isinstance(exclude_fields, list):  exclude_fields = [exclude_fields]
        fields = [f for f in fields if f not in exclude_fields]

    # Create writer
    try:
        db_writer = csv.DictWriter(out_handle, fieldnames=fields, dialect='excel-tab')
        db_writer.writeheader()
    except:
        sys.exit('ERROR:  File %s cannot be written' % out_handle.name)

    return db_writer


# >>> HOW TO CLOSE db_handle?
def readDbFile(db_file, ig=True):
    """
    Reads database files

    Arguments: 
    db_file = a tab delimited database file
    ig = if True convert fields to an IgRecord
    
    Returns: 
    a database record iterator
    """
    # Read and check file
    try:
        db_handle = open(db_file, 'rb')
        db_reader = csv.DictReader(db_handle, dialect='excel-tab')
        if ig:  
            db_iter = (IgRecord(r) for r in db_reader)
        else:  
            db_iter = db_reader
    except IOError:
        sys.exit('ERROR:  File %s cannot be read' % db_file)
    except:
        sys.exit('ERROR:  File %s is invalid' % db_file)
    
    return db_iter


def countDbFile(db_file):
    """
    Counts the records in database files

    Arguments: 
    db_file = a tab delimited database file

    Returns: 
    the count of records in the database file
    """
    # Count records and check file
    try:
        with open(db_file) as db_handle:
            db_records = csv.reader(db_handle, dialect='excel-tab') 
            for i, __ in enumerate(db_records):  pass
        db_count = i
    except IOError:
        sys.exit('ERROR:  File %s cannot be read' % db_file)
    except:
        sys.exit('ERROR:  File %s is invalid' % db_file)
    else:
        if db_count == 0:  sys.exit('ERROR:  File %s is empty' % db_file)
        
    return db_count


if __name__ == '__main__':
    """
    Print module information
    """
    print 'Version: %s %s %s' % (os.path.basename(__file__), __version__, __date__)
    print 'Location: %s' % os.path.dirname(os.path.realpath(__file__))
    #print 'Parent Dir: %s' % path.join(path.dirname(path.realpath(__file__)), path.pardir)
    