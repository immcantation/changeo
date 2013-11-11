#!/usr/bin/env python
"""
Core functions shared by Changeo modules
"""

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2013 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.0'
__date__      = '2013.10.12'

# Imports
import csv, math, os, re, sys
import multiprocessing as mp
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from itertools import izip, izip_longest, product
from collections import OrderedDict
from time import time, strftime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

# Defaults
default_id_field = 'ID'
default_seq_field = 'SEQ'


class IgRecord:
    """
    A class defining a V(D)J germline sequence alignment
    """
    # Private static methods
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
            try:  return str(v.seq)
            except:  return ''
            
    # Private variables
    _parse_map = {'SEQUENCE_ID': ('id', '_identity'),
                  'V_CALL': ('v_call', '_identity'),
                  'D_CALL': ('d_call', '_identity'),
                  'J_CALL': ('j_call', '_identity'),
                  'SEQUENCE': ('seq', '_sequence'), 
                  'SEQUENCE_GAP': ('seq_gap', '_sequence'),
                  'JUNCTION': ('junction', '_sequence'),
                  'FUNCTIONAL': ('functional', '_logical'), 
                  'IN_FRAME': ('in_frame', '_logical'), 
                  'STOP': ('stop', '_logical'), 
                  'MUTATED_INVARIANT': ('mutated_invariant', '_logical'), 
                  'INDELS': ('indels', '_logical'),
                  'V_MATCH': ('v_match', '_integer'),
                  'V_LENGTH': ('v_length', '_integer'),
                  'J_MATCH': ('j_match', '_integer'),
                  'J_LENGTH': ('j_length', '_integer'),
                  'V_GAP_LENGTH': ('v_gap_length', '_integer'),
                  'N1_LENGTH': ('n1_length', '_integer'),
                  'D_5_TRIM': ('d_5_trim', '_integer'),
                  'D_3_TRIM': ('d_3_trim', '_integer'),
                  'N2_LENGTH': ('n2_length', '_integer'),
                  'J_5_TRIM': ('j_5_trim', '_integer'),
                  'J_GAP_LENGTH': ('j_gap_length', '_integer'),
                  'JUNCTION_GAP_LENGTH': ('junction_gap_length', '_integer')}
    _deparse_map = {v[0]: (k, v[1]) for k, v in _parse_map.iteritems()}
    
    _logical_parse = {'F':False, 'T':True}
    _logical_deparse = {False:'F', True:'T'}
    
    # Public variables
    allele_regex = re.compile(r'(IG[HLK][VDJ]\d+[-/\w]*[-\*][\.\w]+)')
    gene_regex = re.compile(r'(IG[HLK][VDJ]\d+[-/\w]*)')
    family_regex = re.compile(r'(IG[HLK][VDJ]\d+)')

    def _parseAllele(self, alleles, regex, action='first'):
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
        try:
            self.id = row.pop('SEQUENCE_ID')
            self.seq = IgRecord._sequence(row.pop('SEQUENCE'))
            self.v_call = row.pop('V_CALL')
            self.d_call = row.pop('D_CALL')
            self.j_call = row.pop('J_CALL')
        except:
            sys.exit('ERROR:  Input must contain valid ID,SEQUENCE,V_CALL,D_CALL,J_CALL values')

        # Defined optional logical values
        self.functional = IgRecord._logical(row.pop('FUNCTIONAL', None))
        self.in_frame = IgRecord._logical(row.pop('IN_FRAME', None))
        self.stop = IgRecord._logical(row.pop('STOP', None))
        self.mutated_invariant = IgRecord._logical(row.pop('MUTATED_INVARIANT', None))
        self.indels = IgRecord._logical(row.pop('INDELS', None))
        
        # Defined optional integer values
        self.v_match = IgRecord._integer(row.pop('V_MATCH', None))
        self.v_length = IgRecord._integer(row.pop('V_LENGTH', None))
        self.j_match = IgRecord._integer(row.pop('J_MATCH', None))
        self.j_length = IgRecord._integer(row.pop('J_LENGTH', None))
        self.v_gap_length = IgRecord._integer(row.pop('V_GAP_LENGTH', None))
        self.n1_length = IgRecord._integer(row.pop('N1_LENGTH', None))
        self.d_5_trim = IgRecord._integer(row.pop('D_5_TRIM', None))
        self.d_3_trim = IgRecord._integer(row.pop('D_3_TRIM', None))
        self.n2_length = IgRecord._integer(row.pop('N2_LENGTH', None))
        self.j_5_trim = IgRecord._integer(row.pop('J_5_TRIM', None))
        self.j_gap_length = IgRecord._integer(row.pop('J_GAP_LENGTH', None))
        self.junction_gap_length = IgRecord._integer(row.pop('JUNCTION_GAP_LENGTH', None))
        
        # Defined option Seq records
        self.junction = IgRecord._sequence(row.pop('JUNCTION', None))
        self.seq_gap = IgRecord._sequence(row.pop('SEQUENCE_GAP', None))
            
        # Add remaining elements as annotations dictionary
        self.annotations = row
    
    # Return a dictionary of the namespace
    def to_dict(self):
        d = {}
        n = self.__dict__
        if 'annotations' in n:
            d.update({k.upper():v for k, v in n['annotations'].iteritems()})
            del n['annotations']
        for k, v in n.iteritems():
            p = IgRecord._deparse_map[k]
            f = getattr(IgRecord, p[1])
            d[p[0]] = f(v, deparse=True)

        return d
        
        
    # Public methods
    def getVAllele(self, action='first'):
        return self._parseAllele(self.v_call, self.allele_regex, action)

    def getDAllele(self, action='first'):
        return self._parseAllele(self.d_call, self.allele_regex, action)

    def getJAllele(self, action='first'):
        return self._parseAllele(self.j_call, self.allele_regex, action)
    
    def getVGene(self, action='first'):
        return self._parseAllele(self.v_call, self.gene_regex, action)

    def getDGene(self, action='first'):
        return self._parseAllele(self.d_call, self.gene_regex, action)

    def getJGene(self, action='first'):
        return self._parseAllele(self.j_call, self.gene_regex, action)
    
    def getVFamily(self, action='first'):
        return self._parseAllele(self.v_call, self.family_regex, action)

    def getDFamily(self, action='first'):
        return self._parseAllele(self.d_call, self.family_regex, action)

    def getJFamily(self, action='first'):
        return self._parseAllele(self.j_call, self.family_regex, action)


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
    