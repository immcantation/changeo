#!/usr/bin/env python
"""
Parses tab delimited database files
"""

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.0'
__date__      = '2014.11.26'

# Imports
import os, sys, textwrap
from argparse import ArgumentParser
from itertools import izip
from collections import OrderedDict
from time import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

# IgCore imports
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from IgCore import default_delimiter, default_out_args
from IgCore import flattenAnnotation 
from IgCore import getOutputHandle, printLog, printProgress
from IgCore import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from DbCore import countDbFile, readDbFile, getDbWriter

# Defaults
default_id_field = 'SEQUENCE_ID'
default_seq_field = 'SEQUENCE_IMGT'
default_germ_field = 'GERMLINE_IMGT_D_MASK'

# TODO:  add regex support for values (with flag)
# TODO:  add partial match support for values (with flag)
# TODO:  convert SQL-ish operations to modify_func() as per ParseHeaders

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


# TODO:  SHOULD ALLOW FOR UNSORTED CLUSTER COLUMN
# TODO:  SHOULD ALLOW FOR GROUPING FIELDS
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


def convertDbFasta(db_file, id_field=default_id_field, seq_field=default_seq_field,
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


def addDbRecords(db_file, fields, values, out_args=default_out_args):
    """
    Adds field and value pairs to a database file

    Arguments:
    db_file = the database file name
    fields = a list of fields to add
    values = a list of values to assign to all rows of each field
    out_args = common output argument dictionary from parseCommonArgs

    Returns:
    the output file name
    """
    log = OrderedDict()
    log['START'] = 'ParseDb'
    log['COMMAND'] = 'add'
    log['FILE'] = os.path.basename(db_file)
    log['FIELDS'] = ','.join(fields)
    log['VALUES'] = ','.join(values)
    printLog(log)

    # Open file handles
    db_iter = readDbFile(db_file, ig=False)
    pass_handle = getOutputHandle(db_file, out_label='parse-add', out_dir=out_args['out_dir'],
                                  out_name=out_args['out_name'], out_type='tab')
    pass_writer = getDbWriter(pass_handle, db_file, add_fields=fields)
    # Count records
    result_count = countDbFile(db_file)

    # Define fields and values to append
    add_dict = {k:v for k,v in izip(fields, values) if k not in db_iter.fieldnames}

    # Iterate over records
    start_time = time()
    rec_count = 0
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time)
        rec_count += 1
        # Write updated row
        rec.update(add_dict)
        pass_writer.writerow(rec)

    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['END'] = 'ParseDb'
    printLog(log)

    # Close file handles
    pass_handle.close()

    return pass_handle.name


def dropDbRecords(db_file, fields, out_args=default_out_args):
    """
    Deletes entire fields from a database file

    Arguments:
    db_file = the database file name
    fields = a list of fields to drop
    out_args = common output argument dictionary from parseCommonArgs

    Returns:
    the output file name
    """
    log = OrderedDict()
    log['START'] = 'ParseDb'
    log['COMMAND'] = 'add'
    log['FILE'] = os.path.basename(db_file)
    log['FIELDS'] = ','.join(fields)
    printLog(log)

    # Open file handles
    db_iter = readDbFile(db_file, ig=False)
    pass_handle = getOutputHandle(db_file, out_label='parse-drop', out_dir=out_args['out_dir'],
                                  out_name=out_args['out_name'], out_type='tab')
    pass_writer = getDbWriter(pass_handle, db_file, exclude_fields=fields)
    # Count records
    result_count = countDbFile(db_file)

    # Iterate over records
    start_time = time()
    rec_count = 0
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time)
        rec_count += 1
        # Write row
        pass_writer.writerow(rec)

    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['END'] = 'ParseDb'
    printLog(log)

    # Close file handles
    pass_handle.close()

    return pass_handle.name


def deleteDbRecords(db_file, fields, values, out_args=default_out_args):
    """
    Deletes records from a database file

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


def selectDbRecords(db_file, fields, values, out_args=default_out_args):
    """
    Selects records from a database file

    Arguments:
    db_file = the database file name
    fields = a list of fields to check for selection criteria
    values = a list of values defining selection targets
    out_args = common output argument dictionary from parseCommonArgs

    Returns:
    the output file name
    """
    log = OrderedDict()
    log['START'] = 'ParseDb'
    log['COMMAND'] = 'select'
    log['FILE'] = os.path.basename(db_file)
    log['FIELDS'] = ','.join(fields)
    log['VALUES'] = ','.join(values)
    printLog(log)

    # Open file handles
    db_iter = readDbFile(db_file, ig=False)
    pass_handle = getOutputHandle(db_file, out_label='parse-select', out_dir=out_args['out_dir'],
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
        select = any([rec.get(f, False) in values for f in fields])

        # Write sequences
        if select:
            pass_count += 1
            pass_writer.writerow(rec)
        else:
            fail_count += 1

    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['SELECTED'] = pass_count
    log['DISCARDED'] = fail_count
    log['END'] = 'ParseDb'
    printLog(log)

    # Close file handles
    pass_handle.close()

    return pass_handle.name


def updateDbRecords(db_file, field, values, updates, out_args=default_out_args):
    """
    Updates field and value pairs to a database file

    Arguments:
    db_file = the database file name
    field = the field to update
    values = a list of values to specifying which rows to update
    updates = a list of values to update each value with
    out_args = common output argument dictionary from parseCommonArgs

    Returns:
    the output file name
    """
    log = OrderedDict()
    log['START'] = 'ParseDb'
    log['COMMAND'] = 'update'
    log['FILE'] = os.path.basename(db_file)
    log['FIELD'] = ','.join(field)
    log['VALUES'] = ','.join(values)
    log['UPDATES'] = ','.join(updates)
    printLog(log)

    # Open file handles
    db_iter = readDbFile(db_file, ig=False)
    pass_handle = getOutputHandle(db_file, out_label='parse-update', out_dir=out_args['out_dir'],
                                  out_name=out_args['out_name'], out_type='tab')
    pass_writer = getDbWriter(pass_handle, db_file)
    # Count records
    result_count = countDbFile(db_file)

    # Iterate over records
    start_time = time()
    rec_count = pass_count = 0
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time)
        rec_count += 1

        # Updated values if found
        for x, y in izip(values, updates):
            if rec[field] == x:
                rec[field] = y
                pass_count += 1

        # Write records
        pass_writer.writerow(rec)

    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['UPDATED'] = pass_count
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
    # Define input and output field help message
    fields = textwrap.dedent(
             '''
             required fields:
                 SEQUENCE_ID
                 
              optional fields:
                 JUNCTION
                 SEQUENCE_IMGT
                 GERMLINE_IMGT
                 GERMLINE_IMGT_D_MASK
                 GERMLINE_IMGT_V_REGION
                 SEQUENCE_VDJ
                 GERMLINE_VDJ
                 GERMLINE_VDJ_D_MASK
                 GERMLINE_VDJ_V_REGION
                
              output fields:
                 None
              ''')
    
    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            version='%(prog)s:' + ' v%s-%s' %(__version__, __date__), 
                            formatter_class=CommonHelpFormatter)
    subparsers = parser.add_subparsers(title='subcommands', dest='command', metavar='',
                                       help='Database operation')

    # Define parent parser
    parser_parent = getCommonArgParser(seq_in=False, seq_out=False, db_in=True,
                                       failed=False, log=False)

    # Subparser to convert database entries to sequence file
    parser_seq = subparsers.add_parser('fasta', parents=[parser_parent],
                                       formatter_class=CommonHelpFormatter,
                                       help='Creates a fasta file from database records')
    parser_seq.add_argument('--if', action='store', dest='id_field', 
                            default=default_id_field,
                            help='The name of the field containing identifiers')
    parser_seq.add_argument('--sf', action='store', dest='seq_field', 
                            default=default_seq_field,
                            help='The name of the field containing sequences')
    parser_seq.add_argument('--mf', nargs='+', action='store', dest='meta_fields',
                            help='List of annotation fields to add to the sequence description')
    parser_seq.set_defaults(func=convertDbFasta)
    
    # Subparser to convert database entries to clip-fasta file
    parser_clip = subparsers.add_parser('clip', parents=[parser_parent], 
                                        formatter_class=CommonHelpFormatter,
                                        help='''Creates a clip-fasta file from database
                                             records, wherein germline sequences precede
                                             each clone and are denoted by ">>" headers.''')
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

    # Subparser to add records
    parser_add = subparsers.add_parser('add', parents=[parser_parent],
                                       formatter_class=CommonHelpFormatter,
                                       help='Adds field and value pairs')
    parser_add.add_argument('-f', nargs='+', action='store', dest='fields', required=True,
                               help='The name of the fields to add.')
    parser_add.add_argument('-u', nargs='+', action='store', dest='values', required=True,
                               help='The value to assign to all rows for each field.')
    parser_add.set_defaults(func=addDbRecords)

    # Subparser to delete records
    parser_delete = subparsers.add_parser('delete', parents=[parser_parent], 
                                          formatter_class=CommonHelpFormatter,
                                          help='Deletes specific records')
    parser_delete.add_argument('-f', nargs='+', action='store', dest='fields', required=True,
                               help='The name of the fields to check for deletion criteria.')
    parser_delete.add_argument('-u', nargs='+', action='store', dest='values', default=['', 'NA'],
                               help='''The values defining which records to delete. A value
                                    may appear in any of the fields specified with -f.''')
    parser_delete.set_defaults(func=deleteDbRecords)

    # Subparser to drop fields
    parser_drop = subparsers.add_parser('drop', parents=[parser_parent],
                                        formatter_class=CommonHelpFormatter,
                                        help='Deletes entire fields')
    parser_drop.add_argument('-f', nargs='+', action='store', dest='fields', required=True,
                               help='The name of the fields to delete from the database.')
    parser_drop.set_defaults(func=dropDbRecords)

    # Subparser to select records
    parser_select = subparsers.add_parser('select', parents=[parser_parent],
                                          formatter_class=CommonHelpFormatter,
                                          help='Selects specific records')
    parser_select.add_argument('-f', nargs='+', action='store', dest='fields', required=True,
                               help='The name of the fields to check for selection criteria.')
    parser_select.add_argument('-u', nargs='+', action='store', dest='values', required=True,
                               help='''The values defining with records to select. A value
                                    may appear in any of the fields specified with -f.''')
    parser_select.set_defaults(func=selectDbRecords)

    # Subparser to update records
    parser_update = subparsers.add_parser('update', parents=[parser_parent],
                                       formatter_class=CommonHelpFormatter,
                                       help='Updates field and value pairs')
    parser_update.add_argument('-f', action='store', dest='field', required=True,
                               help='The name of the field to update.')
    parser_update.add_argument('-u', nargs='+', action='store', dest='values', required=True,
                               help='The values that will be replaced.')
    parser_update.add_argument('-t', nargs='+', action='store', dest='updates', required=True,
                               help='''The new value to assign to each selected row.''')
    parser_update.set_defaults(func=updateDbRecords)

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

    # Check modify_args arguments
    if args.command == 'add' and len(args_dict['fields']) != len(args_dict['values']):
        parser.error('You must specify exactly one value (-u) per field (-f)')
    elif args.command == 'update' and len(args_dict['values']) != len(args_dict['updates']):
        parser.error('You must specify exactly one value (-u) per replacement (-t)')

    # Call parser function for each database file
    del args_dict['command']
    del args_dict['func']
    del args_dict['db_files']
    for f in args.__dict__['db_files']:
        args_dict['db_file'] = f
        args.func(**args_dict)
 