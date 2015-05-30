#!/usr/bin/env python
"""
Performs amino acid analysis of Ig sequences
"""

__author__    = 'Namita Gupta, Daniel Gadala-Maria'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.0'
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


def analyzeAa(db_file, cdr3=True, out_args=default_out_args):
    """
    Calculate amino acid properties for specified regions and add to tab-delimited database

    Arguments:
    db_file = input tab-delimited database file
    cdr3 = true if cdr3 amino acid properties are desired
    out_args = arguments for output preferences

    Returns:
    None
    """
    log = OrderedDict()
    log['START'] = 'AAnalysis'
    log['FILE'] = path.basename(db_file)
    #log['CDR3'] = 'True'
    printLog(log)
    
    # Create reader instance for input file
    reader = readDbFile(db_file, ig=False)
    # Create output handle
    out_handle = getOutputHandle(db_file, 
                                 'aaproperties', 
                                 out_dir=out_args['out_dir'], 
                                 out_name=out_args['out_name'], 
                                 out_type=out_args['out_type'])

    if(cdr3):
        # Create writer instance
        writer = getDbWriter(out_handle, db_file, 
                             add_fields=['CDR3_AA_LENGTH', 'CDR3_AA_POSITIVE', 'CDR3_AA_NEGATIVE', 
                                         'CDR3_ARGININE', 'CDR3_HISTIDINE', 'CDR3_LYSINE', 
                                         'CDR3_TYROSINE', 'CDR3_ALIPHATIC', 'CDR3_AROMATIC', 'CDR3_GRAVY'])
    
    # Initialize time and total count for progress bar
    start_time = time()
    rec_count = countDbFile(db_file)
    
    # Iterate over rows
    for i,row in enumerate(reader):
        # Print progress bar
        printProgress(i, rec_count, 0.05, start_time)
        
        # Calculate CDR3 amino acid properties
        if(cdr3): 
            aacdr3 = cdr3Properties(row['JUNCTION'], out_args)
            for k,v in aacdr3.iteritems(): row[k] = v
            
        # Write row to output
        writer.writerow(row)
        
    # Print log    
    printProgress(i+1, rec_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT'] = out_handle.name
    log['PASS'] = i+1
    log['END'] = 'AAnalysis'
    printLog(log)
    
    # Close handles
    out_handle.close()
    
    
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
                                       failed=False, annotation=False, log=False)
    # Define argument parser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            version='%(prog)s:' + ' v%s-%s' %(__version__, __date__),
                            parents=[parser_parent], 
                            formatter_class=CommonHelpFormatter)
    # parser.add_argument('--cdr3', action='store_true', dest='cdr3', default=True,
    
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
