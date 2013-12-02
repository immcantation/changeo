#!/usr/bin/env python
"""
Parses records in the console log of pRESTO modules
"""

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2013 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.1'
__date__      = '2013.12.2'

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
from DbCore import countDbFile, readDbFile, getDbWriter

# Defaults
default_id_field = 'SEQUENCE_ID'
default_seq_field = 'SEQUENCE_GAP'
default_germ_field = 'GERMLINE_GAP_D_MASK'


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


# >>> SHOULD ALLOW FOR UNSORTED CLUSTER COLUMN
# >>> SHOULD ALLOW FOR GROUPING FIELDS
def convertDbClip(db_file, id_field=default_id_field, seq_field=default_seq_field, 
                  germ_field=default_germ_field, cluster_field=None, 
                  meta_fields=None, out_args=default_out_args):
    """
    Builds fasta files from database records

    Arguments: 
    db_file = the database file name
    id_field = the field containing identifiers
    seq_field = the field containing sample sequences
    germ_field = the field containing germline sequences
    cluster_field = the field containing clonal groupings
                    if None write the germline for each record
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
    log['GERM_FIELD'] = germ_field
    log['CLUSTER_FIELD'] = cluster_field
    if meta_fields is not None:  log['META_FIELDS'] = ','.join(meta_fields)
    printLog(log)
    
    # Open file handles
    db_iter = readDbFile(db_file, ig=False)
    pass_handle = getOutputHandle(db_file, out_label='sequences', out_dir=out_args['out_dir'], 
                                  out_name=out_args['out_name'], out_type='clip')
    # Count records
    result_count = countDbFile(db_file)
    
    # Iterate over records
    start_time = time()
    rec_count = germ_count = pass_count = fail_count = 0
    cluster_last = None
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time)
        rec_count += 1
        
        # Update cluster ID
        cluster = rec.get(cluster_field, None)
        
        # Get germline SeqRecord when needed
        if cluster_field is None:
            germ = getDbSeqRecord(rec, id_field, germ_field, meta_fields, 
                                  delimiter=out_args['delimiter'])
            germ.id = '>' + germ.id
        elif cluster != cluster_last:
            germ = getDbSeqRecord(rec, cluster_field, germ_field, 
                                  delimiter=out_args['delimiter'])
            germ.id = '>' + germ.id            
        else:
            germ = None

        # Get read SeqRecord
        seq = getDbSeqRecord(rec, id_field, seq_field, meta_fields, 
                             delimiter=out_args['delimiter'])
        
        # Write germline
        if germ is not None:
            germ_count += 1
            SeqIO.write(germ, pass_handle, 'fasta')
        
        # Write sequences
        if seq is not None:
            pass_count += 1
            SeqIO.write(seq, pass_handle, 'fasta')
        else:
            fail_count += 1
        
        # Set last cluster ID
        cluster_last = cluster
        
    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['GERMLINES'] = germ_count
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    log['END'] = 'ParseDb'
    printLog(log)

    # Close file handles
    pass_handle.close()
 
    return pass_handle.name


def convertDbSeq(db_file, id_field=default_id_field, seq_field=default_seq_field, 
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


def deleteDbRecords(db_file, fields, values, out_args=default_out_args):
    """
    Builds fasta files from database records

    Arguments: 
    db_file = the database file name
    fields = a list of fields to check for deletion criteria
    values = a list of values defining deletion targets
    out_args = common output argument dictionary from parseCommonArgs
                    
    Returns: 
    the output file name
    """
    log = OrderedDict()
    log['START'] = 'ParseDb'
    log['COMMAND'] = 'delete'
    log['FILE'] = os.path.basename(db_file)
    log['FIELDS'] = ','.join(fields)
    log['VALUES'] = ','.join(values)
    printLog(log)
    
    # Open file handles
    db_iter = readDbFile(db_file, ig=False)
    pass_handle = getOutputHandle(db_file, out_label='parse-delete', out_dir=out_args['out_dir'], 
                                  out_name=out_args['out_name'], out_type='tab')
    pass_writer = getDbWriter(pass_handle, db_file)
    # Count records
    result_count = countDbFile(db_file)
    
    # Iterate over records
    start_time = time()
    rec_count = pass_count = fail_count = 0
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time)
        rec_count += 1

        # Check for deletion values in all fields
        delete = any([rec.get(f, False) in values for f in fields])
        
        # Write sequences
        if not delete:
            pass_count += 1
            pass_writer.writerow(rec)
        else:
            fail_count += 1
        
    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['KEPT'] = pass_count
    log['DELETED'] = fail_count
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
    
    # Define parent parser
    parser_parent = getCommonArgParser(seq_in=False, seq_out=False, db_in=True, log=False)

    # Subparser to convert database entries to sequence file
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
    parser_seq.set_defaults(func=convertDbSeq)
    
    # Subparser to convert database entries to clip-fasta file
    parser_clip = subparsers.add_parser('clip', parents=[parser_parent], 
                                        formatter_class=ArgumentDefaultsHelpFormatter,
                                        help='Creates a fasta file from database records')
    parser_clip.add_argument('--if', action='store', dest='id_field', 
                             default=default_id_field,
                             help='The name of the field containing identifiers')
    parser_clip.add_argument('--sf', action='store', dest='seq_field',
                             default=default_seq_field,
                             help='The name of the field containing reads')
    parser_clip.add_argument('--gf', action='store', dest='germ_field',
                             default=default_germ_field,
                             help='The name of the field containing germline sequences')
    parser_clip.add_argument('--cf', action='store', dest='cluster_field', default=None,
                             help='The name of the field containing containing sorted clone IDs')
    parser_clip.add_argument('--mf', nargs='+', action='store', dest='meta_fields',
                             help='List of annotation fields to add to the sequence description')
    parser_clip.set_defaults(func=convertDbClip)

    # Subparser to delete records
    parser_delete = subparsers.add_parser('delete', parents=[parser_parent], 
                                       formatter_class=ArgumentDefaultsHelpFormatter,
                                       help='Deletes database records')
    parser_delete.add_argument('-f', nargs='+', action='store', dest='fields', required=True,
                               help='The name of the fields to check for deletion criteria')
    parser_delete.add_argument('-u', nargs='+', action='store', dest='values', default=['', 'NA'],
                               help='The values defining with records to delete')
    parser_delete.set_defaults(func=deleteDbRecords)
        
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
    #args_dict['id_field'] = args_dict['id_field'].upper()
    #args_dict['seq_field'] = args_dict['seq_field'].upper() 
    #if args_dict['meta_fields']:  args_dict['meta_fields'] = map(str.upper, args_dict['meta_fields']) 

    # Call parser function for each database file
    del args_dict['func']
    del args_dict['db_files']
    for f in args.__dict__['db_files']:
        args_dict['db_file'] = f
        args.func(**args_dict)
 