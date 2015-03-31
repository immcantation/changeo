#!/usr/bin/env python
"""
Sorts, samples and splits tab-delimited database files
"""

__author__    = 'Namita Gupta, Jason Vander Heiden'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.0'
__date__      = '2014.10.2'

# Imports
import os, sys, textwrap
from argparse import ArgumentParser
from collections import OrderedDict
from time import time

# IgCore and DbCore imports 
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from IgCore import default_out_args
from IgCore import getOutputHandle, printLog, printProgress
from IgCore import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from DbCore import readDbFile, getDbWriter, countDbFile
 

def downsizeDbFile(db_file, max_count, out_args=default_out_args):
    """
    Splits a tab-delimited database file into segments with a limited number of records

    Arguments: 
    db_file = filename of the tab-delimited database file to split
    max_count = number of records in each output file
    out_args = common output argument dictionary from parseCommonArgs

    Returns: 
    a list of output file names
    """
    log = OrderedDict()
    log['START'] = 'SplitTab'
    log['COMMAND'] = 'count'
    log['FILE'] = os.path.basename(db_file) 
    log['MAX_COUNT'] = max_count
    printLog(log)
    
    # Determine total numbers of records
    rec_count = countDbFile(db_file)
    
    # Loop through iterator writing each record
    # Open new output file handle as needed
    start_time = time()
    count, part_num = 0, 1
    reader = readDbFile(db_file)
    out_handle = getOutputHandle(db_file, 'part%06i' % part_num, 
                                 out_type = out_args['out_type'],
                                 out_name = out_args['out_name'], 
                                 out_dir = out_args['out_dir'])
    writer = getDbWriter(out_handle, db_file)
    out_files = [out_handle.name]
    for row in reader:
        # Print progress
        printProgress(count, rec_count, 0.05, start_time)
        
        # Update count
        count += 1
        
        # Write IgRecord row
        writer.writerow(row.toDict())
        # Check if total records reached to avoid extra empty file
        if count == rec_count:
            break
        
        # Open new file if needed
        if count % max_count == 0:
            out_handle.close()
            part_num += 1
            out_handle = getOutputHandle(db_file, 'part%06i' % part_num, 
                                         out_type = out_args['out_type'], 
                                         out_name = out_args['out_name'],
                                         out_dir = out_args['out_dir'])
            writer = getDbWriter(out_handle, db_file)
            out_files.append(out_handle.name)
    
    # Print log
    printProgress(count, rec_count, 0.05, start_time)
    log = OrderedDict()
    for i, f in enumerate(out_files): 
        log['OUTPUT%i' % (i + 1)] = os.path.basename(f)
    log['RECORDS'] = rec_count
    log['PARTS'] = len(out_files)
    log['END'] = 'SplitTab'
    printLog(log)
    
    # Close file handles
    out_handle.close()

    return out_files


def groupDbFile(db_file, field, threshold, out_args=default_out_args):
    """
    Divides a tab-delimited database file into segments by description tags

    Arguments: 
    db_file = filename of the tab-delimited database file to split
    field = The field name by which to split db_file
    threshold = The numerical threshold by which to group sequences;
                if None treat field as textual
    out_args = common output argument dictionary from parseCommonArgs

    Returns: 
    a list of output file names
    """
    log = OrderedDict()
    log['START'] = 'SplitTab'
    log['COMMAND'] = 'group'
    log['FILE'] = os.path.basename(db_file) 
    log['FIELD'] = field
    log['THRESHOLD'] = threshold
    printLog(log)
    
    # Open IgRecord reader iter object
    reader = readDbFile(db_file, ig=False) 

    # Determine total numbers of records
    rec_count = countDbFile(db_file)
    
    start_time = time()
    count = 0
    # Sort records into files based on textual field
    if threshold is None:
        # Create set of unique field tags
        tmp_iter = readDbFile(db_file, ig=False)
        tag_list = list(set([row[field] for row in tmp_iter]))
        
        # Forbidden characters in filename and replacements
        noGood = {'\/':'f','\\':'b','?':'q','\%':'p','*':'s',':':'c',
                  '\|':'pi','\"':'dq','\'':'sq','<':'gt','>':'lt',' ':'_'}
        # replace forbidden characters in tag_list
        tag_dict = {}
        for tag in tag_list:
            for c,r in noGood.iteritems():
                tag_dict[tag] = (tag_dict.get(tag,tag).replace(c,r) if c in tag else tag_dict.get(tag,tag)) 
        
        # Create output handles
        handles_dict = {tag:getOutputHandle(db_file, 
                                            label, 
                                            out_type = out_args['out_type'],
                                            out_name = out_args['out_name'], 
                                            out_dir = out_args['out_dir'])
                        for tag,label in tag_dict.iteritems()}
        
        # Create Db writer instances
        writers_dict = {tag:getDbWriter(handles_dict[tag], db_file)
                        for tag,label in tag_dict.iteritems()}
        
        # Iterate over IgRecords
        for row in reader:
            printProgress(count, rec_count, 0.05, start_time)
            count += 1
            # Write row to appropriate file
            tag = row[field]
            writers_dict[tag].writerow(row)
            
    # Sort records into files based on numeric threshold    
    else:
        threshold = float(threshold)
        
        # Create output handles
        handles_dict = {'under':getOutputHandle(db_file, 
                                                'under-%.1f' % threshold, 
                                                out_type = out_args['out_type'], 
                                                out_name = out_args['out_name'], 
                                                out_dir = out_args['out_dir']),
                        'atleast':getOutputHandle(db_file, 
                                                  'atleast-%.1f' % threshold, 
                                                out_type = out_args['out_type'], 
                                                out_name = out_args['out_name'], 
                                                out_dir = out_args['out_dir'])}
        
        # Create Db writer instances
        writers_dict = {'under':getDbWriter(handles_dict['under'], db_file),
                        'atleast':getDbWriter(handles_dict['atleast'], db_file)}

        # Iterate over IgRecords
        for row in reader:
            printProgress(count, rec_count, 0.05, start_time)
            count += 1
            tag = row[field]
            tag = 'under' if float(tag) < threshold else 'atleast'
            writers_dict[tag].writerow(row)
    
    # Write log
    printProgress(count, rec_count, 0.05, start_time)
    log = OrderedDict()
    for i, k in enumerate(handles_dict): 
        log['OUTPUT%i' % (i + 1)] = os.path.basename(handles_dict[k].name)
    log['RECORDS'] = rec_count
    log['PARTS'] = len(handles_dict)
    log['END'] = 'SplitTab'
    printLog(log)
    
    # Close output file handles
    for t in handles_dict: handles_dict[t].close()

    return [handles_dict[t].name for t in handles_dict]


def sortDbFile(db_file, field, numeric, max_count, out_args=default_out_args):
    """
    Sorts a tab-delimited database file by annotation fields

    Arguments: 
    db_file = filename of the tab-delimited database file to sort
    field = the field name to sort by
    numeric = if True sort field numerically;
              if False sort field alphabetically
    max_count = maximum number of records in each output file
                if None do not create multiple files
    out_args = common output argument dictionary from parseCommonArgs
    
    Returns: 
    a list of output file names
    """
    log = OrderedDict()
    log['START'] = 'SplitTab'
    log['COMMAND'] = 'sort'
    log['FILE'] = os.path.basename(db_file)
    log['FIELD'] = field
    log['NUMERIC'] = numeric
    log['MAX_COUNT'] = max_count
    printLog(log)
    
    # Save all IgRecords in file as dictionary
    reader = readDbFile(db_file, ig=False)
    db_dict = {}
    for row in reader:
        db_dict[row['SEQUENCE_ID']] = row
        
    # Get fields and sort db_dict by annotation values
    tag_dict = {k:v[field] for k,v in db_dict.iteritems()}
    if numeric:  tag_dict = {k:float(v or 0) for k,v in tag_dict.iteritems()}
    sorted_ids = sorted(tag_dict, key=tag_dict.get)

    # Determine total numbers of records
    rec_count = len(db_dict)
    if max_count >= rec_count: max_count = None
    
    # Loop through sorted sequence IDs writing each record
    # Open initial output file handle and Db writer
    part_num = 1
    if max_count is None: out_label = 'sorted'
    else: out_label = 'sorted-part%06i' % part_num
    out_handle = getOutputHandle(db_file, 
                                 out_label, 
                                 out_dir=out_args['out_dir'], 
                                 out_name=out_args['out_name'], 
                                 out_type=out_args['out_type'])
    writer = getDbWriter(out_handle, db_file)
    out_files = [out_handle.name]
    
    # Loop through sorted IgRecord dictionary keys  
    start_time = time()  
    saved_tag, saved_ids = None, []
    count = chunk_count = 0
    for seq_id in sorted_ids:
        # Print progress and update count
        printProgress(count, rec_count, 0.05, start_time)
        count += 1
        
        # Write saved group of sequences when tag changes or end reached
        if saved_tag is not None and tag_dict[seq_id] != saved_tag:
            if max_count is not None and chunk_count + len(saved_ids) > max_count:
                # Update partition counts
                part_num += 1
                # Open new file handle and add name to output file list
                out_handle.close()
                out_handle = getOutputHandle(db_file, 
                                 out_label, 
                                 out_dir=out_args['out_dir'], 
                                 out_name=out_args['out_name'], 
                                 out_type=out_args['out_type'])
                writer = getDbWriter(out_handle, db_file)
                out_files.append(out_handle.name)
            
            # Write saved records  
            for k in saved_ids:
                chunk_count += 1
                writer.writerow(db_dict[k])
            # Reset saved keys to current key only
            saved_ids = [seq_id]
            
        # Update list of saved tags if tag is unchanged
        else:
            saved_ids.append(seq_id)
            
        # Check if total records reached, write all saved IDs, and exit loop
        if count == rec_count:
            for k in saved_ids:
                chunk_count += 1
                writer.writerow(db_dict[k])
            out_handle.close()
            break
        
        # Update tag tracker
        saved_tag = tag_dict[seq_id]
        
    # Print log
    printProgress(count, rec_count, 0.05, start_time)
    log = OrderedDict()
    for i, f in enumerate(out_files): 
        log['OUTPUT%i' % (i + 1)] = os.path.basename(f)
    log['RECORDS'] = rec_count
    log['PARTS'] = len(out_files)
    log['END'] = 'SplitTab'
    printLog(log)
    
    # Close file handles
    out_handle.close()
    
    return out_files


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
                
              output fields:
                 None
              ''')
    
    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            version='%(prog)s:' + ' v%s-%s' %(__version__, __date__), 
                            formatter_class=CommonHelpFormatter)
    subparsers = parser.add_subparsers(title='subcommands', dest='command', 
                                       help='Splitting operation mode', metavar='')
    
    # Parent parser    
    parser_parent = getCommonArgParser(seq_in=False, seq_out=False, db_in=True, 
                                       failed=False, annotation=False, log=False)


    # Subparser to downsize files to a maximum count
    parser_count = subparsers.add_parser('count', parents=[parser_parent],
                                         formatter_class=CommonHelpFormatter,
                                         help='Splits sequences files by number of records')
    parser_count.add_argument('-c', action='store', dest='max_count', type=int, required=True,
                              help='Maximum number of records in each new file')
    parser_count.set_defaults(func=downsizeDbFile)
    
    # Subparser to partition files by annotation
    parser_group = subparsers.add_parser('group', parents=[parser_parent],
                                         formatter_class=CommonHelpFormatter,
                                         help='Splits tab-delimited database files by annotation')
    parser_group.add_argument('-f', action='store', dest='field', type=str, required=True,
                              help='Annotation field by which to split database files')
    parser_group.add_argument('--num', action='store', dest='threshold', type=float, default=None, 
                              help='Specify to define the split column as numeric and group records by value')
    parser_group.set_defaults(func=groupDbFile)
    
    # Subparser to sort files
    parser_sort = subparsers.add_parser('sort', parents=[parser_parent],
                                        formatter_class=CommonHelpFormatter,
                                        help='Sorts tab-delimited database files by annotation')
    parser_sort.add_argument('-f', action='store', dest='field', type=str, required=True,
                             help='The annotation field by which to sort records')
    parser_sort.add_argument('--num', action='store_true', dest='numeric', default=False,
                             help='Specify to define the sort column as numeric rather than textual')
    parser_sort.add_argument('--max', action='store', dest='max_count', type=int,
                             default=None, help='Maximum number of records in each new file')
    parser_sort.set_defaults(func=sortDbFile)
    
    return parser


if __name__ == '__main__':
    """
    Parses command line arguments and calls appropriate sequencing parsing function
    """
    # Parse arguments
    parser = getArgParser()
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)
    
    # Call appropriate function for each sample file
    del args_dict['command']
    del args_dict['func']    
    del args_dict['db_files']
    for f in args.__dict__['db_files']:
        args_dict['db_file'] = f
        args.func(**args_dict)
