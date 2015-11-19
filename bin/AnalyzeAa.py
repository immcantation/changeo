#!/usr/bin/env python3
"""
Performs amino acid analysis of Ig sequences
"""
# Info
__author__ = 'Namita Gupta, Daniel Gadala-Maria'
from changeo import __version__, __date__

# Imports
import re
from argparse import ArgumentParser
from collections import OrderedDict
from os import path
from textwrap import dedent
from time import time
from Bio.Seq import Seq

# Presto and changeo imports
from presto.Defaults import default_out_args
from presto.Commandline import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from presto.IO import getOutputHandle, printLog, printProgress
from changeo.IO import getDbWriter, readDbFile, countDbFile

# Defaults
default_seq_field = 'JUNCTION'


def gravy(aa_seq):
    """
    Calculate GRAVY (Grand Average of Hydropathy) index for amino acid sequence
    (http://web.expasy.org/tools/protparam/protparam-doc.html)

    Arguments:
    aa_seq = amino acid sequence for which index is to be calculated

    Returns:
    GRAVY index
    """
    hydropath = {"A":1.8, "R":-4.5, "N":-3.5, "D":-3.5, "C":2.5, "Q":-3.5, 
                 "E":-3.5, "G":-0.4, "H":-3.2, "I":4.5, "L":3.8, "K":-3.9,
                 "M":1.9, "F":2.8, "P":-1.6, "S":-0.8, "T":-0.7, "W":-0.9,
                 "Y":-1.3, "V":4.2}
    g = sum([hydropath.get(c,0) for c in aa_seq])
    return g/len(aa_seq)
    

def AaProperties(seq):
    """
    Calculate amino acid properties of a sequence

    Arguments:
    seq = input nucleotide sequence

    Returns:
    dictionary with amino acid properties
    """
    # TODO: needs a better solution to the gap character problem at some point
    seq = re.sub('\.|-', 'N', seq)
    aa_seq = str(Seq(seq).translate())

    props = {}
    
    # Calculate length
    props['_AA_LENGTH'] = len(aa_seq)
    
    # Count the percent of aa that are positively charged
    props['_AA_POSITIVE'] = round(100*float(len(re.findall("[RK]", aa_seq)))/props['_AA_LENGTH'], 2)
    
    # Count percent of aa that are negatively charged
    props['_AA_NEGATIVE'] = round(100*float(len(re.findall("[DE]", aa_seq)))/props['_AA_LENGTH'], 2)
    
    # Count the percent of aa that are Arg
    props['_ARGININE'] = round(100*float(len(re.findall("[R]", aa_seq)))/props['_AA_LENGTH'], 2)
    
    # Count percent of aa that are His
    props['_HISTIDINE'] = round(100*float(len(re.findall("[H]", aa_seq)))/props['_AA_LENGTH'], 2)
    
    # Count the percent of aa that are Lys
    props['_LYSINE'] = round(100*float(len(re.findall("[K]", aa_seq)))/props['_AA_LENGTH'], 2)
    
    # Count percent of aa that are Tyr
    props['_TYROSINE'] = round(100*float(len(re.findall("[Y]", aa_seq)))/props['_AA_LENGTH'], 2)
    
    # Aliphatic index
    # Some documentation: http://web.expasy.org/tools/protparam/protparam-doc.html
    nAla    = len(re.findall("[A]", aa_seq))
    nVal    = len(re.findall("[V]", aa_seq))
    nLeuIle = len(re.findall("[LI]", aa_seq))
    a = 2.9
    b = 3.9
    props['_ALIPHATIC'] = round(100*float(nAla + a*nVal + b*nLeuIle)/props['_AA_LENGTH'], 2)
    
    # Percent aromatic AAs
    props['_AROMATIC'] = round(100*float(len(re.findall("[FWHY]", aa_seq)))/props['_AA_LENGTH'], 2)
    
    # GRAVY (Grand Average of Hydropathy) index
    # Some documentation: http://web.expasy.org/tools/protparam/protparam-doc.html
    props['_GRAVY'] = round(gravy(aa_seq), 2)
    
    return props


def analyzeAa(db_file, seq_field=default_seq_field, out_args=default_out_args):
    """
    Calculate amino acid properties for specified regions and add to database file

    Arguments:
    db_file = input tab-delimited database file
    seq_field = sequence field for which amino acid properties are analyzed
    out_args = arguments for output preferences

    Returns:
    None
    """
    log = OrderedDict()
    log['START'] = 'AnalyzeAa'
    log['FILE'] = path.basename(db_file)
    log['SEQ_FIELD'] = seq_field
    printLog(log)

    # Define reader instance for input file
    reader = readDbFile(db_file, ig=False)
    # Define passed output handle and writer
    pass_handle = getOutputHandle(db_file,
                                  'aa-pass',
                                  out_dir=out_args['out_dir'],
                                  out_name=out_args['out_name'],
                                  out_type=out_args['out_type'])
    prop_fields = [seq_field + p for p in ['_AA_LENGTH', '_AA_POSITIVE',
                                           '_AA_NEGATIVE','_ARGININE',
                                           '_HISTIDINE', '_LYSINE', '_TYROSINE',
                                           '_ALIPHATIC', '_AROMATIC', '_GRAVY']]
    pass_writer = getDbWriter(pass_handle, db_file,
                         add_fields=prop_fields)
    # Define failed output handle and writer
    if out_args['failed']:
        fail_handle = getOutputHandle(db_file,
                                      out_label='aa-fail',
                                      out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'],
                                      out_type=out_args['out_type'])
        fail_writer = getDbWriter(fail_handle, db_file)
    else:
        fail_handle = None
        fail_writer = None

    # Initialize time and total count for progress bar
    start_time = time()
    rec_count = countDbFile(db_file)

    # Iterate over rows
    pass_count = 0
    fail_count = 0
    for i,row in enumerate(reader):
        # Print progress bar
        printProgress(i, rec_count, 0.05, start_time)

        # Check that sequence field is not empty and has length a multiple of three
        if(row[seq_field] != '' and len(row[seq_field])%3 == 0):
            # Calculate amino acid properties
            aaprops = AaProperties(row[seq_field])
            for k,v in aaprops.items(): row[seq_field + k] = v

            pass_count += 1
            # Write row to pass file
            pass_writer.writerow(row)
        else:
            fail_count += 1
            # Write row to fail file
            if fail_writer is not None:
                fail_writer.writerow(row)
        
    # Print log    
    printProgress(i+1, rec_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT'] = pass_handle.name
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    log['END'] = 'AnalyzeAa'
    printLog(log)

    # Close file handles
    pass_handle.close()
    if fail_handle is not None: fail_handle.close()
    
    
def getArgParser():
    """
    Defines the ArgumentParser

    Arguments: 
    None
                      
    Returns: 
    an ArgumentParser object
    """
    # Define input and output field help message
    fields = dedent(
             '''
             output files:
               aa-pass             database with amino acid properties.
               aa-fail             database with records failing analysis.

             required fields:
               <user defined>      sequence field specified by the --sf parameter
                
             output fields:
               <user defined>_AA_LENGTH
               <user defined>_AA_POSITIVE
               <user defined>_AA_NEGATIVE
               <user defined>_ARGININE
               <user defined>_HISTIDINE
               <user defined>_LYSINE
               <user defined>_TYROSINE
               <user defined>_ALIPHATIC
               <user defined>_AROMATIC
               <user defined>_GRAVY
              ''')
                  
    # Parent parser    
    parser_parent = getCommonArgParser(seq_in=False, seq_out=False, db_in=True,
                                       failed=True, annotation=False, log=False)
    # Define argument parser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            parents=[parser_parent], 
                            formatter_class=CommonHelpFormatter)
    parser.add_argument('--version', action='version',
                        version='%(prog)s:' + ' %s-%s' %(__version__, __date__))
    parser.add_argument('--sf', action='store', dest='seq_field',
                        default=default_seq_field,
                        help='The name of the field to be analyzed')
    
    return parser


if __name__ == "__main__":
    """
    Parses command line arguments and calls main
    """
    # Parse arguments
    parser = getArgParser()
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)
    args_dict['seq_field'] = args_dict['seq_field'].upper()
    
    del args_dict['db_files']
    
    for f in args.__dict__['db_files']:
        args_dict['db_file'] = f
        analyzeAa(**args_dict)
