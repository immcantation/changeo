#!/usr/bin/env python
"""
Performs amino acid analysis of Ig sequences
"""

__author__    = 'Namita Gupta, Daniel Gadala-Maria'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.2.0'
__date__      = '2015.05.30'

# Imports
import re, textwrap
from os import path
from argparse import ArgumentParser
from Bio.Seq import Seq
from collections import OrderedDict
from time import time

# IgCore and DbCore imports 
from sys import path as syspath
syspath.append(path.dirname(path.realpath(__file__)))
from IgCore import default_out_args
from IgCore import getOutputHandle, printLog, printProgress
from IgCore import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from DbCore import getDbWriter, countDbFile, readDbFile

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
    

def cdr3Properties(junc, out_args):
    """
    Calculate amino acid properties of CDR3

    Arguments:
    junc = input junction nucleotide sequence
    out_args = arguments for output preferences

    Returns:
    dictionary with CDR3 amino acid properties
    """    
    # Trim junction to CDR3
    cdr3_in = junc[3:len(junc)-3].upper()
    # TODO: needs a better solution to the gap character problem at some point
    cdr3_in = re.sub('\.|-', 'N', cdr3_in)

    # Remove sequences that are too short to translate
    not_empty = cdr3_in if len(cdr3_in) > 2 else ''
    cdr3_aa = str(Seq(not_empty).translate())
    

    cdr3 = {}
    
    # Calculate CDR3 Lengths
    cdr3['CDR3_AA_LENGTH'] = len(cdr3_aa)
    
    # Count the percent of aa that are positively charged
    cdr3['CDR3_AA_POSITIVE'] = round(100*float(len(re.findall("[RK]", cdr3_aa)))/cdr3['CDR3_AA_LENGTH'], 2)
    
    # Count percent of aa that are negatively charged
    cdr3['CDR3_AA_NEGATIVE'] = round(100*float(len(re.findall("[DE]", cdr3_aa)))/cdr3['CDR3_AA_LENGTH'], 2)
    
    # Count the percent of aa that are Arg
    cdr3['CDR3_ARGININE'] = round(100*float(len(re.findall("[R]", cdr3_aa)))/cdr3['CDR3_AA_LENGTH'], 2)
    
    # Count percent of aa that are His
    cdr3['CDR3_HISTIDINE'] = round(100*float(len(re.findall("[H]", cdr3_aa)))/cdr3['CDR3_AA_LENGTH'], 2)
    
    # Count the percent of aa that are Lys
    cdr3['CDR3_LYSINE'] = round(100*float(len(re.findall("[K]", cdr3_aa)))/cdr3['CDR3_AA_LENGTH'], 2)
    
    # Count percent of aa that are Tyr
    cdr3['CDR3_TYROSINE'] = round(100*float(len(re.findall("[Y]", cdr3_aa)))/cdr3['CDR3_AA_LENGTH'], 2)
    
    # Aliphatic index
    # Some documentation: http://web.expasy.org/tools/protparam/protparam-doc.html
    nAla    = len(re.findall("[A]", cdr3_aa))
    nVal    = len(re.findall("[V]", cdr3_aa))
    nLeuIle = len(re.findall("[LI]", cdr3_aa))
    a = 2.9
    b = 3.9
    cdr3['CDR3_ALIPHATIC'] = round(100*float(nAla + a*nVal + b*nLeuIle)/cdr3['CDR3_AA_LENGTH'], 2)
    
    # Percent CDR3 Aromatic AAs
    cdr3['CDR3_AROMATIC'] = round(100*float(len(re.findall("[FWHY]", cdr3_aa)))/cdr3['CDR3_AA_LENGTH'], 2)
    
    # GRAVY (Grand Average of Hydropathy) index
    # Some documentation: http://web.expasy.org/tools/protparam/protparam-doc.html
    cdr3['CDR3_GRAVY'] = round(gravy(cdr3_aa), 2)
    
    return cdr3


def analyzeAa(db_file, seq_field=default_seq_field, out_args=default_out_args):
    """
    Calculate amino acid properties for specified regions and add to tab-delimited database

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
    
    # Create reader instance for input file
    reader = readDbFile(db_file, ig=False)
    # Create passed output handle and writer
    pass_handle = getOutputHandle(db_file,
                                 'aa',
                                 out_dir=out_args['out_dir'], 
                                 out_name=out_args['out_name'], 
                                 out_type=out_args['out_type'])
    pass_writer = getDbWriter(pass_handle, db_file,
                         add_fields=['CDR3_AA_LENGTH', 'CDR3_AA_POSITIVE', 'CDR3_AA_NEGATIVE',
                                     'CDR3_ARGININE', 'CDR3_HISTIDINE', 'CDR3_LYSINE',
                                     'CDR3_TYROSINE', 'CDR3_ALIPHATIC', 'CDR3_AROMATIC', 'CDR3_GRAVY'])

    # Defined failed output handle and writer
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
            aacdr3 = cdr3Properties(row[seq_field], out_args)
            for k,v in aacdr3.iteritems(): row[k] = v

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
    if fail_handle is not None:  fail_handle.close()
    
    
def getArgParser():
    """
    Defines the ArgumentParser

    Arguments: 
    None
                      
    Returns: 
    an ArgumentParser object
    """
    # Define input and output field help message
    fields = textwrap.dedent(
             '''
             required fields:
                 JUNCTION
                
              output fields:
                 CDR3_AA_LENGTH
                 CDR3_AA_POSITIVE
                 CDR3_AA_NEGATIVE 
                 CDR3_ARGININE
                 CDR3_HISTIDINE
                 CDR3_LYSINE
                 CDR3_TYROSINE
                 CDR3_ALIPHATIC
                 CDR3_AROMATIC
                 CDR3_GRAVY
              ''')
                  
    # Parent parser    
    parser_parent = getCommonArgParser(seq_in=False, seq_out=False, db_in=True,
                                       failed=True, annotation=False, log=False)
    # Define argument parser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            version='%(prog)s:' + ' v%s-%s' %(__version__, __date__),
                            parents=[parser_parent], 
                            formatter_class=CommonHelpFormatter)

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
    
    del args_dict['db_files']
    
    for f in args.__dict__['db_files']:
        args_dict['db_file'] = f
        analyzeAa(**args_dict)
