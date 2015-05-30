#!/usr/bin/env python
"""
Core functions shared by Change-O modules
"""

__author__    = 'Jason Anthony Vander Heiden, Namita Gupta'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.0'
__date__      = '2015.05.30'

# Imports
import csv, os, re, sys
import pandas as pd
from itertools import product
from collections import OrderedDict
from time import time
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# Core imports
from IgCore import getScoreDict, scoreDNA, scoreAA
from IgCore import getOutputHandle, getFileType
from IgCore import printLog, printProgress

# Defaults
#default_repo = 'germlines'

# Regular expression globals
allele_regex = re.compile(r'((IG[HLK]|TR[ABGD])([VDJ]\d+[-/\w]*[-\*][\.\w]+))')
gene_regex = re.compile(r'((IG[HLK]|TR[ABGD])([VDJ]\d+[-/\w]*))')
family_regex = re.compile(r'((IG[HLK]|TR[ABGD])([VDJ]\d+))')

v_allele_regex = re.compile(r'((IG[HLK]|TR[ABGD])V\d+[-/\w]*[-\*][\.\w]+)')
d_allele_regex = re.compile(r'((IG[HLK]|TR[ABGD])D\d+[-/\w]*[-\*][\.\w]+)')
j_allele_regex = re.compile(r'((IG[HLK]|TR[ABGD])J\d+[-/\w]*[-\*][\.\w]+)')

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
                'v_germ_start': 'V_GERM_START',
                'v_germ_length': 'V_GERM_LENGTH',
                'n1_length': 'N1_LENGTH',
                'd_seq_start': 'D_SEQ_START',
                'd_seq_length': 'D_SEQ_LENGTH',
                'd_germ_start': 'D_GERM_START',
                'd_germ_length': 'D_GERM_LENGTH',
                'n2_length': 'N2_LENGTH',
                'j_seq_start': 'J_SEQ_START',
                'j_seq_length': 'J_SEQ_LENGTH',
                'j_germ_start': 'J_GERM_START',
                'j_germ_length': 'J_GERM_LENGTH',
                'junction_length': 'JUNCTION_LENGTH'}

    # Mapping of column names to member variables
    _field_map = {v: k for k, v in _key_map.iteritems()}

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
                  'v_germ_start': '_integer',
                  'v_germ_length': '_integer',
                  'n1_length': '_integer',
                  'd_seq_start': '_integer',
                  'd_seq_length': '_integer',
                  'd_germ_start': '_integer',
                  'd_germ_length': '_integer',
                  'n2_length': '_integer',
                  'j_seq_start': '_integer',
                  'j_seq_length': '_integer',
                  'j_germ_start': '_integer',
                  'j_germ_length': '_integer',
                  'junction_length': '_integer'}

    _logical_parse = {'F':False, 'T':True, 'TRUE':True, 'FALSE':False, 'NA':None, 'None':None}
    _logical_deparse = {False:'F', True:'T', None:'None'}

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

    # Initializer
    #
    # Arguments:  row = dictionary of {field:value} data
    #             genotyped = if True assign v_call from genotyped field
    # Returns:    IgRecord
    def __init__(self, row, genotyped=True):
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
    
    # Get a field value by column name
    #
    # Arguments:  field = column name
    # Returns:    value in the field
    def getField(self, field):
        if field in IgRecord._field_map:
            return getattr(self, IgRecord._field_map[field])
        elif field in self.annotations:
            return self.annotations[field]
        else:
            return None

    # Returns: dictionary of the namespace
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
    

    # Methods to get multiple allele, gene and family calls
    #
    # Arguments:  calls = iterable of calls to get; one or more of ('v','d','j')
    #             actions = one of ('first','set')
    # Returns:    list of requested calls in order
    def getAlleleCalls(self, calls, action='first'):
        vdj = {'v': self.getVAllele(action),
               'd': self.getDAllele(action),
               'j': self.getJAllele(action)}
        return [vdj[k] for k in calls]

    def getGeneCalls(self, calls, action='first'):
        vdj = {'v':self.getVGene(action),
               'd':self.getDGene(action),
               'j':self.getJGene(action)}
        return [vdj[k] for k in calls]

    def getFamilyCalls(self, calls, action='first'):
        vdj = {'v':self.getVFamily(action),
               'd':self.getDFamily(action),
               'j':self.getJFamily(action)}
        return [vdj[k] for k in calls]

    # Individual allele, gene and family getter methods
    #
    # Arguments:  actions = one of ('first','set')
    # Returns:    call as a string
    def getVAllele(self, action='first'):
        # TODO: this can't distinguish empty value ("") from missing field (no column)
        x = self.v_call_geno if self.v_call_geno is not None else self.v_call
        return parseAllele(x, allele_regex, action)

    def getDAllele(self, action='first'):
        return parseAllele(self.d_call, allele_regex, action)

    def getJAllele(self, action='first'):
        return parseAllele(self.j_call, allele_regex, action)
    
    def getVGene(self, action='first'):
        return parseAllele(self.v_call, gene_regex, action)

    def getDGene(self, action='first'):
        return parseAllele(self.d_call, gene_regex, action)

    def getJGene(self, action='first'):
        return parseAllele(self.j_call, gene_regex, action)
    
    def getVFamily(self, action='first'):
        return parseAllele(self.v_call, family_regex, action)

    def getDFamily(self, action='first'):
        return parseAllele(self.d_call, family_regex, action)

    def getJFamily(self, action='first'):
        return parseAllele(self.j_call, family_regex, action)


class DbData:
    """
    A class defining IgRecord data objects for worker processes
    """
    # Instantiation
    def __init__(self, key, records):
        self.id = key
        self.data = records
        self.valid = (key is not None and records is not None)

    # Boolean evaluation
    def __nonzero__(self):
        return self.valid

    # Length evaluation
    def __len__(self):
        if isinstance(self.data, IgRecord):
            return 1
        elif self.data is None:
            return 0
        else:
            return len(self.data)


class DbResult:
    """
    A class defining IgRecord result objects for collector processes
    """
    # Instantiation
    def __init__(self, key, records):
        self.id = key
        self.data = records
        self.results = None
        self.valid = False
        self.log = OrderedDict([('ID', key)])
        #if isinstance(values, list):
        #    for v in values:  setattr(self, v, None)
        #else:
        #    setattr(self, values, None)

    # Boolean evaluation
    def __nonzero__(self):
        return self.valid

    # Length evaluation
    def __len__(self):
        if isinstance(self.results, IgRecord):
            return 1
        elif self.results is None:
            return 0
        else:
            return len(self.results)

    # Set data_count to number of data records
    @property
    def data_count(self):
        if isinstance(self.data, IgRecord):
            return 1
        elif self.data is None:
            return 0
        else:
            return len(self.data)


# TODO:  might be cleaner as getAllele(), getGene(), getFamily()
def parseAllele(alleles, regex, action='first'):
    """
    Extract alleles from strings

    Arguments:  alleles = string with allele calls
                regex = compiled regular expression for allele match
                action = action to perform for multiple alleles;
                         one of ('first', 'set', 'list').
    Returns:    string of the allele for action='first';
                tuple of allele calls for 'set' or 'list' actions.
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


# TODO:  change to require output fields rather than in_file? probably better that way.
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
        # >>> THIS NEEDS TO BE FIXED, extrasaction='ignore' IS A WORKAROUND FOR ADDITIONS TO IgRecord
        db_writer = csv.DictWriter(out_handle, fieldnames=fields, dialect='excel-tab', extrasaction='ignore')
        db_writer.writeheader()
    except:
        sys.exit('ERROR:  File %s cannot be written' % out_handle.name)

    return db_writer


# TODO:  Need to close db_handle?
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


def getDistMat(mat=None, n_score=0, gap_score=0, alphabet='dna'):
    """
    Generates a distance matrix

    Arguments:
    mat = input distance matrix to extend to full alphabet;
          if unspecified, creates Hamming distance matrix that incorporates IUPAC equivalencies
    n_score = score for all matches against an N character
    gap_score = score for all matches against a [-, .] character
    alphabet = the type of score dictionary to generate;
               one of [dna, aa] for DNA and amino acid characters

    Returns:
    a distance matrix (pandas DataFrame)
    """
    if alphabet=='dna':
        IUPAC_chars = list('-.ACGTRYSWKMBDHVN')
        n = 'N'
        score_func = scoreDNA
    elif alphabet=='aa':
        IUPAC_chars = list('-.*ABCDEFGHIJKLMNOPQRSTUVWXYZ')
        n = 'X'
        score_func = scoreAA
    else:
        sys.stderr.write('ERROR:  The alphabet %s is not a recognized type.\n' % alphabet)

    # Default matrix to inf
    dist_mat = pd.DataFrame(float('inf'), index=IUPAC_chars, columns=IUPAC_chars, dtype=float)
    # Set gap score
    for c in '-.':
        dist_mat.loc[c] = dist_mat.loc[:,c] = gap_score
    # Set n score
    dist_mat.loc[n] = dist_mat.loc[:,n] = n_score
    # Fill in provided distances from input matrix
    if mat is not None:
        for i,j in product(mat.index, mat.columns):
            dist_mat.loc[i,j] = mat.loc[i,j]
    # If no input matrix, create IUPAC-defined Hamming distance
    else:
        for i,j in product(dist_mat.index, dist_mat.columns):
            dist_mat.loc[i,j] = 1 - score_func(i, j, n_score=1-n_score, gap_score=1-gap_score)

    return dist_mat


def feedDbQueue(alive, data_queue, db_file, group_func=None, group_args={}):
    """
    Feeds the data queue with Ig records

    Arguments:
    alive = a multiprocessing.Value boolean controlling whether processing continues
            if False exit process
    data_queue = a multiprocessing.Queue to hold data for processing
    db_file = the Ig record database file
    group_func = the function to use for grouping records
    group_args = a dictionary of arguments to pass to group_func

    Returns:
    None
    """
    # Open input file and perform grouping
    try:
        # Iterate over Ig records and assign groups
        db_iter = readDbFile(db_file)
        if group_func is not None:
            group_dict = group_func(db_iter, **group_args)
            group_iter = group_dict.iteritems()
        else:
            group_iter = ((r.id, r) for r in db_iter)
    except:
        alive.value = False
        raise

    # Add groups to data queue
    try:
        # Iterate over groups and feed data queue
        while alive.value:
            # Get data from queue
            if data_queue.full():  continue
            else:  data = next(group_iter, None)
            # Exit upon reaching end of iterator
            if data is None:  break

            # Feed queue
            data_queue.put(DbData(*data))
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
    except:
        #sys.stderr.write('Exception in feeder queue feeding step\n')
        alive.value = False
        raise

    return None


def processDbQueue(alive, data_queue, result_queue, process_func, process_args={}):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments:
    alive = a multiprocessing.Value boolean controlling whether processing
            continues; when False function returns
    data_queue = a multiprocessing.Queue holding data to process
    result_queue = a multiprocessing.Queue to hold processed results
    process_func = the function to use for filtering sequences
    process_args = a dictionary of arguments to pass to process_func

    Returns:
    None
    """
    try:
        # Iterator over data queue until sentinel object reached
        while alive.value:
            # Get data from queue
            if data_queue.empty():  continue
            else:  data = data_queue.get()
            # Exit upon reaching sentinel
            if data is None:  break

            # Perform work
            result = process_func(data, **process_args)

            # Feed results to result queue
            result_queue.put(result)
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
    except:
        alive.value = False
        sys.stderr.write('Error processing sequence with ID: %s.\n' % data.id)
        raise

    return None


def collectDbQueue(alive, result_queue, collect_queue, db_file, task_label, out_args,
                   add_fields=None):
    """
    Pulls from results queue, assembles results and manages log and file IO

    Arguments:
    alive = a multiprocessing.Value boolean controlling whether processing
            continues; when False function returns
    result_queue = a multiprocessing.Queue holding worker results
    collect_queue = a multiprocessing.Queue to store collector return values
    db_file = the database file name
    task_label = the task label used to tag the output files
    out_args = common output argument dictionary from parseCommonArgs
    add_fields = a list of fields added to the writer not present in the in_file;
                 if None do not add fields

    Returns:
    None
    (adds a dictionary of {log: log object, out_files: output file names} to collect_queue)
    """
    try:
        result_count = countDbFile(db_file)

        # Define output format
        out_type = getFileType(db_file) if out_args['out_type'] is None \
                   else out_args['out_type']

        # Defined valid alignment output handle
        pass_handle = getOutputHandle(db_file,
                                      '%s-pass' % task_label,
                                      out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'],
                                      out_type=out_type)
        pass_writer = getDbWriter(pass_handle, db_file, add_fields=add_fields)
        # Defined failed alignment output handle
        if out_args['failed']:
            fail_handle = getOutputHandle(db_file,
                                          '%s-fail'  % task_label,
                                          out_dir=out_args['out_dir'],
                                          out_name=out_args['out_name'],
                                          out_type=out_type)
            fail_writer = getDbWriter(fail_handle, db_file)
        else:
            fail_handle = None

        # Define log handle
        if out_args['log_file'] is None:
            log_handle = None
        else:
            log_handle = open(out_args['log_file'], 'w')
    except:
        alive.value = False
        raise

    try:
        # Iterator over results queue until sentinel object reached
        start_time = time()
        set_count = rec_count = pass_count = fail_count = 0
        while alive.value:
            # Get result from queue
            if result_queue.empty():  continue
            else:  result = result_queue.get()
            # Exit upon reaching sentinel
            if result is None:  break

            # Print progress for previous iteration
            printProgress(pass_count, result_count, 0.05, start_time)

            # Update counts for current iteration
            set_count += 1
            rec_count += result.data_count

            # Write log
            printLog(result.log, handle=log_handle)

            # Write alignments
            if result:
                pass_count += result.data_count
                for rec in result.results:
                    pass_writer.writerow(rec.toDict())
            else:
                fail_count += result.data_count
                if fail_handle is not None:
                    for rec in result.data:
                        pass_writer.writerow(rec.toDict())
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None

        # Print total counts
        printProgress(pass_count, result_count, 0.05, start_time)

        # Update return values
        log = OrderedDict()
        log['OUTPUT'] = os.path.basename(pass_handle.name)
        log['RECORDS'] = rec_count
        log['GROUPS'] = set_count
        log['PASS'] = pass_count
        log['FAIL'] = fail_count
        collect_dict = {'log':log, 'out_files': [pass_handle.name]}
        collect_queue.put(collect_dict)

        # Close file handles
        pass_handle.close()
        if fail_handle is not None:  fail_handle.close()
        if log_handle is not None:  log_handle.close()
    except:
        alive.value = False
        raise

    return None


if __name__ == '__main__':
    """
    Print module information
    """
    print 'Version: %s %s %s' % (os.path.basename(__file__), __version__, __date__)
    print 'Location: %s' % os.path.dirname(os.path.realpath(__file__))
    #print 'Parent Dir: %s' % path.join(path.dirname(path.realpath(__file__)), path.pardir)

