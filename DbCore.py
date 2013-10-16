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
default_delimiter = ('|', '=', ',')
default_separator = default_delimiter[2]
default_action_choices = ('min', 'max', 'sum', 'first', 'last', 'set')
default_coord_type = 'presto'
default_missing_chars = ('-', '.', 'N')
default_out_args = {'log_file':None, 
                    'delimiter':default_delimiter,
                    'separator':default_separator,
                    'out_dir':None,
                    'out_name':None,
                    'out_type':None,
                    'clean':False}



class IgRecord:
    """
    A class defining a V(D)J germline sequence alignment
    """
    # Public variables
    allele_regex = re.compile(r'(IG[HLK][VDJ]\d+[-/\w]*[-\*][\.\w]+)')
    gene_regex = re.compile(r'(IG[HLK][VDJ]\d+[-/\w]*)')
    family_regex = re.compile(r'(IG[HLK][VDJ]\d+)')
    
    # Private methods
    def _logical(self, v):
        trans = {'F':False, 'T':True}
        try:  return trans[v]
        except:  return None

    def _integer(self, v):
        try:  return int(v)
        except:  return None
        
    def _sequence(self, v):
        try:  return Seq(v, IUPAC.ambiguous_dna)
        except:  return None

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
            self.seq = Seq(row.pop('SEQUENCE'), IUPAC.ambiguous_dna)
            self.v_call = row.pop('V_CALL')
            self.d_call = row.pop('D_CALL')
            self.j_call = row.pop('J_CALL')
        except:
            sys.exit('ERROR:  Input must contain valid ID,SEQUENCE,V_CALL,D_CALL,J_CALL values')

        # Defined optional logical values
        self.functional = self._logical(row.pop('FUNCTIONAL', None))
        self.in_frame = self._logical(row.pop('IN_FRAME', None))
        self.stop = self._logical(row.pop('STOP', None))
        self.mutated_invariant = self._logical(row.pop('MUTATED_INVARIANT', None))
        self.indels = self._logical(row.pop('INDELS', None))
        
        # Defined optional integer values
        self.v_match = self._integer(row.pop('V_MATCH', None))
        self.v_length = self._integer(row.pop('V_LENGTH', None))
        self.j_match = self._integer(row.pop('J_MATCH', None))
        self.j_length = self._integer(row.pop('J_LENGTH', None))
        self.v_gap_length = self._integer(row.pop('V_GAP_LENGTH', None))
        self.n1_length = self._integer(row.pop('N1_LENGTH', None))
        self.d_5_trim = self._integer(row.pop('D_5_TRIM', None))
        self.d_3_trim = self._integer(row.pop('D_3_TRIM', None))
        self.n2_length = self._integer(row.pop('N2_LENGTH', None))
        self.j_5_trim = self._integer(row.pop('J_5_TRIM', None))
        self.j_gap_length = self._integer(row.pop('J_GAP_LENGTH', None))
        self.junction_gap_length = self._integer(row.pop('JUNCTION_GAP_LENGTH', None))
        
        # Defined option Seq records
        self.junction = self._sequence(row.pop('JUNCTION', None))
        self.seq_gap = self._sequence(row.pop('SEQUENCE_GAP', None))
            
        # Add remaining elements as annotations dictionary
        self.annotations = row
    
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
    
    
def getFileType(filename):
    """
    Determines the type of a file by file extension

    Arguments: 
    filename = a filename
    
    Returns: 
    a string defining the sequence type for SeqIO operations
    """
    # Read and check file
    try:
        file_type = os.path.splitext(filename)[1].lower().lstrip('.')
    except IOError:
        sys.exit('ERROR:  File %s cannot be read' % filename)
    except:
        sys.exit('ERROR:  File %s is invalid' % filename)
    else:
        if file_type not in ['fasta', 'fastq', 'embl', 'gb', 'tab']:
            sys.exit('ERROR:  File %s has an unrecognized type' % filename)
    
    return file_type


def readDbFile(db_file):
    """
    Reads database files

    Arguments: 
    db_file = a tab delimited database file
    
    Returns: 
    a tuple of database record iterator
    """
    # Read and check file
    try:
        db_handle = open(db_file, 'rb')
        db_records = csv.DictReader(db_handle, dialect='excel-tab')
        db_iter = (IgRecord(r) for r in db_records)
    except IOError:
        sys.exit('ERROR:  File %s cannot be read' % db_file)
    except:
        sys.exit('ERROR:  File %s is invalid' % db_file)
    
    return db_iter


def countDbRecords(db_file):
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


def getOutputHandle(in_file, out_label=None, out_dir=None, out_name=None, out_type=None):
    """
    Opens an output file handle
    
    Arguments: 
    in_file = the input filename
    out_label = text to be inserted before the file extension;
                if None do not add a label
    out_type = the file extension of the output file; 
               if None use input file extension
    out_dir = the output directory;
              if None use directory of input file
    out_name = the short file anme to use for the output file; 
               if None use input file short name
    
    Returns:
    a file handle
    """
    # Get in_file components
    dir_name, file_name = os.path.split(in_file)
    short_name, ext_name = os.path.splitext(file_name)
    
    # Define output directory
    if out_dir is None:
        out_dir = dir_name
    else:
        out_dir = os.path.abspath(out_dir)
        if not os.path.exists(out_dir):  os.mkdir(out_dir)
    # Define output file prefix
    if out_name is None:  out_name = short_name
    # Define output file extension
    if out_type is None:  out_type = ext_name.lstrip('.')

    # Define output file name
    if out_label is None:  
        out_file = os.path.join(out_dir, '%s.%s' % (out_name, out_type))
    else:
        out_file = os.path.join(out_dir, '%s_%s.%s' % (out_name, out_label, out_type))
    
    return open(out_file, 'w')


def printLog(record, handle=sys.stdout, inset=None):
    """
    Formats a dictionary into an IgPipeline log string

    Arguments: 
    record = a dict or OrderedDict of {field names: values}
    handle = the file handle to write the log to;
             if None do not write to file
    inset = minimum field name inset;
            if None automatically space field names
                    
    Returns:
    a formatted multi-line string in IgPipeline log format
    """
    # Return empty string if empty dictionary passed
    if not record:
        return ''
    
    # Determine inset
    max_len = max(map(len, record))
    inset = max(max_len, inset)
    
    # Assemble log string
    record_str = ''
    if isinstance(record, OrderedDict):
        key_list = record.keys()
    else:
        key_list = sorted(record)
    for key in key_list:
        record_str += '%s%s> %s\n' % (' ' * (inset - len(key)), key, record[key])
    
    # Write log record
    if handle is not None:
        try:
            handle.write('%s\n' % record_str)
        except IOError as e:
            sys.stderr.write('I/O error writing to log file: %s\n' % e.strerror)

    return record_str


def printProgress(current, total=None, step=None, start_time=None, end=False):
    """
    Prints a progress bar to standard out
    
    Arguments:
    current = the count of completed tasks 
    total = the total task count;
            if None do not print percentage
    step = a float defining the fractional progress increment to print if total is defined;
           an int defining the progress increment to print at if total is not defined;
           if None always output the progress
    start_time = task start time returned by time.time();
                 if None do not add run time to progress
    end = if True print final log (add newline)
    
    Returns:
    None
    """
    try:
        # Check update condition
        if total is None:
            update = (current % step == 0)
        else:
            update = (current % math.ceil(step*total) == 0)
    except:
        # Return on modulo by zero error
        return None
    else:
        # Check end condition and return if no update needed
        if current == total:
            end = True
        if not end and not update:
            return None
        
    # Define progress bar
    if total is not None and total != 0:
        p = float(current) / total
        c = format(current, "%i,d" % len(format(total, ",d")))
        bar = 'PROGRESS> %s [%-20s] %3.0f%% (%s)' \
              % (strftime('%H:%M:%S'), '#' * int(p*20), p*100, c)
    else:
        bar = 'PROGRESS> %s (%s)' % (strftime('%H:%M:%S'), current)
        
    # Add run time to bar if start_time is specified
    if start_time is not None:
        bar = '%s %.1f min' % (bar, (time() - start_time)/60)
    
    # Print progress bar    
    if current == 0:
        print '%s' % bar,
        sys.stdout.flush()
    elif end:
        print '\r%s\n' % bar
        sys.stdout.flush()
    else:
        print '\r%s' % bar,
        sys.stdout.flush()
    
    return None


def getCommonParser(seq_in=True, seq_out=True, paired=False, db_in=False, 
                    annotation=True, log=True, multiproc=False):
    """
    Defines an ArgumentParser object with common pRESTO arguments

    Arguments: 
    seq_in = if True include sequence input arguments
    seq_out = if True include sequence output arguments
    paired = if True defined paired-end sequence input and output arguments
    db_in = if True include tab delimited database input arguments
    annotation = if True include annotation arguments
    log = if True include log arguments
    multiproc = if True include multiprocessing arguments
    
    Returns:
    an ArgumentParser object
    """
    parser = ArgumentParser(add_help=False, formatter_class=ArgumentDefaultsHelpFormatter)

    # Database arguments
    if db_in:
        parser.add_argument('-d', nargs='+', action='store', dest='db_files', required=True,
                        help='Tab delimited database files containing Ig records')
            
    # Sequence arguments
    if seq_in and not paired:
        parser.add_argument('-s', nargs='+', action='store', dest='seq_files', required=True,
                            help='List of FASTA/FASTQ files containing sequences')
    elif seq_in and paired:
        parser.add_argument('-1', nargs='+', action='store', dest='seq_files_1', required=True,
                            help='Ordered list of FASTA/FASTQ files containing head/primary sequences')
        parser.add_argument('-2', nargs='+', action='store', dest='seq_files_2', required=True,
                            help='Ordered list of FASTA/FASTQ files containing tail/secondary sequences')
    if seq_out:
        parser.add_argument('--fasta', action='store_const', dest='out_type', const='fasta',
                            help='Specify to force output as FASTA rather than FASTQ')
        parser.add_argument('--clean', action='store_true', dest='clean', 
                            help='If specified do not create files of sequences that fail processing')
            
    # Annotation arguments
    if annotation:
        parser.add_argument('--delim', nargs=3, action='store', dest='delimiter', 
                            type=str, default=default_delimiter, 
                            help='The three delimiters that separate annotation blocks, \
                                  field names and values, and values within a list, respectively')

    # Log arguments
    if log:
        parser.add_argument('--log', action='store', dest='log_file', default=None,
                            help='Specify to write verbose logging to a file')    
    
    # Multiprocessing arguments
    if multiproc:
        parser.add_argument('--nproc', action='store', dest='nproc', type=int, default=mp.cpu_count(),
                            help='The number of computational simultaneous computational processes to execute')
    
    # Universal arguments
    parser.add_argument('--outdir', action='store', dest='out_dir', default=None,
                        help='Changes the output directory to the path specified \
                              Uses input file directory if not specified')
    parser.add_argument('--outname', action='store', dest='out_name', default=None,
                        help='Changes the prefix of the successfully processed output file to the string specified')
    
    return parser


def parseCommonArgs(args, file_args=None):
    """
    Checks common arguments from getCommonParser and transforms output options to a dictionary

    Arguments: 
    args = argument Namespace defined by ArgumentParser.parse_args
    file_args = a list of non-common file arguments to verify; 
                ['seq_files', '1', '2', 'primer_file', 'out_dir'] 
                are checked by default
                    
    Returns:
    a dictionary copy of args with output arguments embedded in the dictionary out_args
    """ 
    db_types = ['.tab']
    seq_types = ['.fasta', '.fastq']
    primer_types = ['.fasta', '.regex']
    args_dict = args.__dict__.copy()
    
    # Verify database files
    db_files = []
    if 'db_files' in args_dict:
        db_files.extend(args_dict['db_files'])
    for f in db_files:
        if not os.path.isfile(f):
            sys.exit('ERROR:  database file %s does not exist' % f)
        if os.path.splitext(f)[-1].lower() not in db_types:
            sys.exit('ERROR:  database file %s is not a supported type. Must be one: %s' \
                     % (','.join(db_types), f))
                
    # Verify sequence files
    seq_files = []
    if 'seq_files' in args_dict:
        seq_files.extend(args_dict['seq_files'])
    elif 'seq_files_1' and 'seq_files_2' in args_dict:
        if len(args_dict['seq_files_1']) != len(args_dict['seq_files_2']):
            sys.exit('ERROR:  The -1 and -2 arguments must contain the same number of files')
        for f1, f2 in zip(args_dict['seq_files_1'], args_dict['seq_files_2']):
            if os.path.splitext(f1)[-1].lower() != os.path.splitext(f2)[-1].lower():
                sys.exit('ERROR:  Each pair of files in the -1 and -2 arguments must be the same file type')
        seq_files.extend(args_dict['seq_files_1'])
        seq_files.extend(args_dict['seq_files_2'])
    for f in seq_files:
        if not os.path.isfile(f):
            sys.exit('ERROR:  sequence file %s does not exist' % f)
        if os.path.splitext(f)[-1].lower() not in seq_types:
            sys.exit('ERROR:  sequence file %s is not a supported type. Must be one: %s' \
                     % (','.join(seq_types), f))

    # Verify primer files
    if 'primer_file' in args_dict:
        primer_files = args_dict['primer_file'] if isinstance(args_dict['primer_file'], list) \
                       else [args_dict['primer_file']]
        for f in primer_files:
            if not os.path.isfile(f):
                sys.exit('ERROR:  primer file %s does not exist' % f)
            if os.path.splitext(f)[-1].lower() not in primer_types:
                sys.exit('ERROR:  primer file %s is not a supported type. Must be one: %s' \
                         % (','.join(primer_types), f))
            
    # Verify output directory
    if 'out_dir' in args_dict and args_dict['out_dir'] is not None:
        if os.path.exists(args_dict['out_dir']) and not os.path.isdir(args_dict['out_dir']):
            sys.exit('ERROR:  directory %s exists and is not a directory' % args_dict['out_dir'])
    
    # Verify optional files
    if file_args is not None:
        if not isinstance(file_args, list):  file_args = [file_args]
        for arg in file_args:
            files = args_dict[arg] if isinstance(args_dict[arg], list) \
                    else [args_dict[arg]]
            for f in files:
                if not os.path.isfile(f):
                    sys.exit('ERROR:  file %s does not exist' % f)
    
    # Redefine common output options as out_args dictionary
    out_args = ['log_file', 'delimiter', 'out_dir', 'out_name', 'out_type', 'clean']
    args_dict['out_args'] = {k:args_dict.setdefault(k, None) for k in out_args}
    for k in out_args: del args_dict[k]
    
    return args_dict


if __name__ == '__main__':
    """
    Print module information
    """
    print 'Version: %s %s %s' % (os.path.basename(__file__), __version__, __date__)
    print 'Location: %s' % os.path.dirname(os.path.realpath(__file__))
    #print 'Parent Dir: %s' % path.join(path.dirname(path.realpath(__file__)), path.pardir)
    