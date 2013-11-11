#!/usr/bin/env python
"""
Parses records in the console log of pRESTO modules
"""

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2013 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.1'
__date__      = '2013.11.10'

# Imports
import os, sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import OrderedDict
from time import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

# IgCore imports
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from IgCore import flattenAnnotation 
from IgCore import getOutputHandle, printLog, printProgress
from IgCore import default_delimiter, default_out_args
from IgCore import getCommonArgParser, parseCommonArgs
from DbCore import default_id_field, default_seq_field
from DbCore import countDbFile, readDbFile


def getDbSeqRecord(db_record, id_field, seq_field, meta_fields=None, 
                   delimiter=default_delimiter):
    """
    Parses a database record into a SeqRecord

    Arguments: 
    db_record = a dictionary containing a database record
    id_field = the field containing identifiers
    seq_field = the field containing sequences
    meta_fields = a list of fields to add to sequence annotations
    delimiter = a tuple of delimiters for (fields, values, value lists) 

    Returns: 
    a SeqRecord
    """
    # Return None if ID or sequence fields are empty
    if not db_record[id_field] or not db_record[seq_field]:
        return None
    
    # Create description string
    desc_dict = OrderedDict([('ID', db_record[id_field])])
    if meta_fields is not None:
        desc_dict.update([(f, db_record[f]) for f in meta_fields if f in db_record]) 
    desc_str = flattenAnnotation(desc_dict, delimiter=delimiter)
    
    # Create SeqRecord
    seq_record = SeqRecord(Seq(db_record[seq_field], IUPAC.ambiguous_dna),
                           id=desc_str, name=desc_str, description='')
        
    return seq_record


def parseDbSeq(db_file, id_field=default_id_field, seq_field=default_seq_field, 
               meta_fields=None, out_args=default_out_args):
    """
    Builds fasta files from database records

    Arguments: 
    db_file = the database file name
    id_field = the field containing identifiers
    seq_field = the field containing sequences
    meta_fields = a list of fields to add to sequence annotations
    out_args = common output argument dictionary from parseCommonArgs
                    
    Returns: 
    the output file name
    """
    log = OrderedDict()
    log['START'] = 'ParseDb'
    log['COMMAND'] = 'fasta'
    log['FILE'] = os.path.basename(db_file)
    log['ID_FIELD'] = id_field
    log['SEQ_FIELD'] = seq_field
    if meta_fields is not None:  log['META_FIELDS'] = ','.join(meta_fields)
    printLog(log)
    
    # Open file handles
    out_type = 'fasta'
    db_iter = readDbFile(db_file, ig=False)
    pass_handle = getOutputHandle(db_file, out_label='sequences', out_dir=out_args['out_dir'], 
                                  out_name=out_args['out_name'], out_type=out_type)
    # Count records
    result_count = countDbFile(db_file)
    
    # Iterate over records
    start_time = time()
    rec_count = pass_count = fail_count = 0
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time)
        rec_count += 1

        # Get SeqRecord
        seq = getDbSeqRecord(rec, id_field, seq_field, meta_fields, out_args['delimiter'])

        # Write sequences
        if seq is not None:
            pass_count += 1
            SeqIO.write(seq, pass_handle, out_type)
        else:
            fail_count += 1
        
    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    log['END'] = 'ParseDb'
    printLog(log)

    # Close file handles
    pass_handle.close()
 
    return pass_handle.name


def getArgParser():
    """
    Defines the ArgumentParser

    Arguments: 
    None
                      
    Returns: 
    an ArgumentParser object
    """
    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, version='%(prog)s:' + ' v%s-%s' %(__version__, __date__), 
                            formatter_class=ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers()
    
    # Subparser to add header fields
    parser_parent = getCommonArgParser(seq_in=False, seq_out=False, db_in=True, log=False)
    parser_seq = subparsers.add_parser('seq', parents=[parser_parent], 
                                       formatter_class=ArgumentDefaultsHelpFormatter,
                                       help='Creates a fasta file from database records')
    parser_seq.add_argument('--if', action='store', dest='id_field', 
                              default=default_id_field,
                              help='The name of the field containing identifiers')
    parser_seq.add_argument('--sf', action='store', dest='seq_field', 
                              default=default_seq_field,
                              help='The name of the field containing sequences')
    parser_seq.add_argument('--mf', nargs='+', action='store', dest='meta_fields',
                              help='List of annotation fields to add to the sequence description')
    parser_seq.set_defaults(func=parseDbSeq)
    
    return parser

    
if __name__ == '__main__':
    """
    Parses command line arguments and calls main function
    """
    # Parse arguments
    parser = getArgParser()
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)
    # Convert case of fields
    args_dict['id_field'] = args_dict['id_field'].upper()
    args_dict['seq_field'] = args_dict['seq_field'].upper() 
    if args_dict['meta_fields']:  args_dict['meta_fields'] = map(str.upper, args_dict['meta_fields']) 

    # Call parser function for each database file
    del args_dict['func']
    del args_dict['db_files']
    for f in args.__dict__['db_files']:
        args_dict['db_file'] = f
        args.func(**args_dict)
 