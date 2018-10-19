#!/usr/bin/env python3
"""
Parses tab delimited database files
"""

# Info
__author__ = 'Jason Anthony Vander Heiden'
from changeo import __version__, __date__

# Imports
import csv
import os
import re
from argparse import ArgumentParser
from collections import OrderedDict
from itertools import chain
from textwrap import dedent
from time import time

# Presto and changeo imports
from presto.IO import printLog, printProgress, printMessage
from changeo.Defaults import default_id_field, default_seq_field, default_germ_field, \
                             default_csv_size, default_out_args
from changeo.Commandline import CommonHelpFormatter, checkArgs, getCommonArgParser, parseCommonArgs
from changeo.IO import countDbFile, getOutputHandle, splitName, TSVReader, TSVWriter

# System settings
csv.field_size_limit(default_csv_size)

# Defaults
default_index_field = 'INDEX'


# TODO:  convert SQL-ish operations to modify_func() as per ParseHeaders
def splitDbFile(db_file, field, num_split=None, out_args=default_out_args):
    """
    Divides a tab-delimited database file into segments by description tags

    Arguments:
      db_file : filename of the tab-delimited database file to split
      field : the field name by which to split db_file
      num_split : the numerical threshold by which to group sequences;
                  if None treat field as textual
      out_args : common output argument dictionary from parseCommonArgs

    Returns:
      list : a list of output file names.
    """
    log = OrderedDict()
    log['START'] = 'ParseDb'
    log['COMMAND'] = 'split'
    log['FILE'] = os.path.basename(db_file)
    log['FIELD'] = field
    log['NUM_SPLIT'] = num_split
    printLog(log)

    # Open input
    db_handle = open(db_file, 'rt')
    db_iter = TSVReader(db_handle)
    out_fields = db_iter.fields
    __, __, out_args['out_type'] = splitName(db_file)

    # Determine total numbers of records
    rec_count = countDbFile(db_file)

    start_time = time()
    count = 0
    # Sort records into files based on textual field
    if num_split is None:
        # Create set of unique field tags
        with open(db_file, 'rt') as tmp_handle:
            tmp_iter = TSVReader(tmp_handle)
            tag_list = list(set([row[field] for row in tmp_iter]))

        # Forbidden characters in filename and replacements
        no_good = {'\/':'f','\\':'b','?':'q','\%':'p','*':'s',':':'c',
                   '\|':'pi','\"':'dq','\'':'sq','<':'gt','>':'lt',' ':'_'}
        # Replace forbidden characters in tag_list
        tag_dict = {}
        for tag in tag_list:
            for c,r in no_good.items():
                tag_dict[tag] = (tag_dict.get(tag, tag).replace(c,r) \
                                 if c in tag else tag_dict.get(tag, tag))

        # Create output handles
        handles_dict = {tag: getOutputHandle(db_file,
                                             out_label='%s-%s' % (field, label),
                                             out_name=out_args['out_name'],
                                             out_dir=out_args['out_dir'],
                                             out_type=out_args['out_type'])
                        for tag, label in tag_dict.items()}

        # Create Db writer instances
        writers_dict = {tag: TSVWriter(handles_dict[tag], fields=out_fields)
                        for tag in tag_dict}

        # Iterate over records
        for row in db_iter:
            printProgress(count, rec_count, 0.05, start_time=start_time)
            count += 1
            # Write row to appropriate file
            tag = row[field]
            writers_dict[tag].writeDict(row)

    # Sort records into files based on numeric num_split
    else:
        num_split = float(num_split)

        # Create output handles
        handles_dict = {'under': getOutputHandle(db_file,
                                                 out_label='under-%.1f' % num_split,
                                                 out_name=out_args['out_name'],
                                                 out_dir=out_args['out_dir'],
                                                 out_type=out_args['out_type']),
                        'atleast': getOutputHandle(db_file,
                                                   out_label='atleast-%.1f' % num_split,
                                                   out_name=out_args['out_name'],
                                                   out_dir=out_args['out_dir'],
                                                   out_type=out_args['out_type'])}

        # Create Db writer instances
        writers_dict = {'under': TSVWriter(handles_dict['under'], fields=out_fields),
                        'atleast': TSVWriter(handles_dict['atleast'], fields=out_fields)}

        # Iterate over records
        for row in db_iter:
            printProgress(count, rec_count, 0.05, start_time=start_time)
            count += 1
            tag = row[field]
            tag = 'under' if float(tag) < num_split else 'atleast'
            writers_dict[tag].writeDict(row)

    # Write log
    printProgress(count, rec_count, 0.05, start_time=start_time)
    log = OrderedDict()
    for i, k in enumerate(handles_dict):
        log['OUTPUT%i' % (i + 1)] = os.path.basename(handles_dict[k].name)
    log['RECORDS'] = rec_count
    log['PARTS'] = len(handles_dict)
    log['END'] = 'ParseDb'
    printLog(log)

    # Close output file handles
    db_handle.close()
    for t in handles_dict: handles_dict[t].close()

    return [handles_dict[t].name for t in handles_dict]


def addDbFile(db_file, fields, values, out_file=None, out_args=default_out_args):
    """
    Adds field and value pairs to a database file

    Arguments:
      db_file : the database file name.
      fields : a list of fields to add.
      values : a list of values to assign to all rows of each field.
      out_file : output file name. Automatically generated from the input file if None.
      out_args : common output argument dictionary from parseCommonArgs.

    Returns:
      str : output file name.
    """
    log = OrderedDict()
    log['START'] = 'ParseDb'
    log['COMMAND'] = 'add'
    log['FILE'] = os.path.basename(db_file)
    log['FIELDS'] = ','.join(fields)
    log['VALUES'] = ','.join(values)
    printLog(log)

    # Open inut
    db_handle = open(db_file, 'rt')
    db_iter = TSVReader(db_handle)
    __, __, out_args['out_type'] = splitName(db_file)

    # Add fields
    out_fields = list(db_iter.fields)
    out_fields.extend(fields)

    # Open output
    if out_file is not None:
        pass_handle = open(out_file, 'w')
    else:
        pass_handle = getOutputHandle(db_file, out_label='parse-add', out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'], out_type=out_args['out_type'])
    pass_writer = TSVWriter(pass_handle, out_fields)

    # Count records
    result_count = countDbFile(db_file)

    # Define fields and values to append
    add_dict = {k:v for k,v in zip(fields, values) if k not in db_iter.fields}

    # Iterate over records
    start_time = time()
    rec_count = 0
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time=start_time)
        rec_count += 1
        # Write updated row
        rec.update(add_dict)
        pass_writer.writeDict(rec)

    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time=start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['END'] = 'ParseDb'
    printLog(log)

    # Close file handles
    pass_handle.close()
    db_handle.close()

    return pass_handle.name


def indexDbFile(db_file, field=default_index_field, out_file=None, out_args=default_out_args):
    """
    Adds an index column to a database file

    Arguments:
      db_file : the database file name.
      field : the name of the index field to add.
      out_file : output file name. Automatically generated from the input file if None.
      out_args : common output argument dictionary from parseCommonArgs.

    Returns:
      str : output file name.
    """
    log = OrderedDict()
    log['START'] = 'ParseDb'
    log['COMMAND'] = 'index'
    log['FILE'] = os.path.basename(db_file)
    log['FIELD'] = field
    printLog(log)

    # Open input
    db_handle = open(db_file, 'rt')
    db_iter = TSVReader(db_handle)
    __, __, out_args['out_type'] = splitName(db_file)

    # Append index field
    out_fields = list(db_iter.fields)
    out_fields.append(field)

    # Open output
    if out_file is not None:
        pass_handle = open(out_file, 'w')
    else:
        pass_handle = getOutputHandle(db_file, out_label='parse-index', out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'], out_type=out_args['out_type'])
    pass_writer = TSVWriter(pass_handle, out_fields)

    # Count records
    result_count = countDbFile(db_file)

    # Iterate over records
    start_time = time()
    rec_count = 0
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time=start_time)
        rec_count += 1

        # Add count and write updated row
        rec.update({field:rec_count})
        pass_writer.writeDict(rec)

    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time=start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['END'] = 'ParseDb'
    printLog(log)

    # Close file handles
    pass_handle.close()
    db_handle.close()

    return pass_handle.name


def dropDbFile(db_file, fields, out_file=None, out_args=default_out_args):
    """
    Deletes entire fields from a database file

    Arguments:
      db_file : the database file name.
      fields : a list of fields to drop.
      out_file : output file name. Automatically generated from the input file if None.
      out_args : common output argument dictionary from parseCommonArgs

    Returns:
     str : output file name.
    """
    log = OrderedDict()
    log['START'] = 'ParseDb'
    log['COMMAND'] = 'add'
    log['FILE'] = os.path.basename(db_file)
    log['FIELDS'] = ','.join(fields)
    printLog(log)

    # Open input
    db_handle = open(db_file, 'rt')
    db_iter = TSVReader(db_handle)
    __, __, out_args['out_type'] = splitName(db_file)

    # Exclude dropped field from output
    out_fields = [f for f in db_iter.fields if f not in fields]

    # Open output
    if out_file is not None:
        pass_handle = open(out_file, 'w')
    else:
        pass_handle = getOutputHandle(db_file, out_label='parse-drop', out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'], out_type=out_args['out_type'])
    pass_writer = TSVWriter(pass_handle, out_fields)

    # Count records
    result_count = countDbFile(db_file)

    # Iterate over records
    start_time = time()
    rec_count = 0
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time=start_time)
        rec_count += 1
        # Write row
        pass_writer.writeDict(rec)

    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time=start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['END'] = 'ParseDb'
    printLog(log)

    # Close file handles
    pass_handle.close()

    return pass_handle.name


def deleteDbFile(db_file, fields, values, logic='any', regex=False,
                 out_file=None, out_args=default_out_args):
    """
    Deletes records from a database file

    Arguments: 
      db_file : the database file name.
      fields : a list of fields to check for deletion criteria.
      values : a list of values defining deletion targets.
      logic : one of 'any' or 'all' defining whether one or all fields must have a match.
      regex : if False do exact full string matches; if True allow partial regex matches.
      out_file : output file name. Automatically generated from the input file if None.
      out_args : common output argument dictionary from parseCommonArgs.
                    
    Returns: 
      str : output file name.
    """
    # Define string match function
    if regex:
        def _match_func(x, patterns):  return any([re.search(p, x) for p in patterns])
    else:
        def _match_func(x, patterns):  return x in patterns

    # Define logic function
    if logic == 'any':
        _logic_func = any
    elif logic == 'all':
        _logic_func = all

    log = OrderedDict()
    log['START'] = 'ParseDb'
    log['COMMAND'] = 'delete'
    log['FILE'] = os.path.basename(db_file)
    log['FIELDS'] = ','.join(fields)
    log['VALUES'] = ','.join(values)
    printLog(log)

    # Open input
    db_handle = open(db_file, 'rt')
    db_iter = TSVReader(db_handle)
    out_fields = db_iter.fields
    __, __, out_args['out_type'] = splitName(db_file)

    # Open output
    if out_file is not None:
        pass_handle = open(out_file, 'w')
    else:
        pass_handle = getOutputHandle(db_file, out_label='parse-delete', out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'], out_type=out_args['out_type'])
    pass_writer = TSVWriter(pass_handle, out_fields)

    # Count records
    result_count = countDbFile(db_file)

    # Iterate over records
    start_time = time()
    rec_count, pass_count, fail_count = 0, 0, 0
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time=start_time)
        rec_count += 1
        # Check for deletion values in all fields
        delete = _logic_func([_match_func(rec.get(f, False), values) for f in fields])
        
        # Write sequences
        if not delete:
            pass_count += 1
            pass_writer.writeDict(rec)
        else:
            fail_count += 1
        
    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time=start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['KEPT'] = pass_count
    log['DELETED'] = fail_count
    log['END'] = 'ParseDb'
    printLog(log)

    # Close file handles
    pass_handle.close()
    db_handle.close()
 
    return pass_handle.name


def renameDbFile(db_file, fields, names, out_file=None, out_args=default_out_args):
    """
    Renames fields in a database file

    Arguments:
      db_file : the database file name.
      fields : a list of fields to rename.
      values : a list of new names for fields.
      out_file : output file name. Automatically generated from the input file if None.
      out_args : common output argument dictionary from parseCommonArgs.

    Returns:
      str : output file name.
    """
    log = OrderedDict()
    log['START'] = 'ParseDb'
    log['COMMAND'] = 'rename'
    log['FILE'] = os.path.basename(db_file)
    log['FIELDS'] = ','.join(fields)
    log['NAMES'] = ','.join(names)
    printLog(log)

    # Open file handles
    db_handle = open(db_file, 'rt')
    db_iter = TSVReader(db_handle)
    __, __, out_args['out_type'] = splitName(db_file)

    # Get header and rename fields
    out_fields = list(db_iter.fields)
    for f, n in zip(fields, names):
        i = out_fields.index(f)
        out_fields[i] = n

    # Open writer
    if out_file is not None:
        pass_handle = open(out_file, 'w')
    else:
        pass_handle = getOutputHandle(db_file, out_label='parse-rename', out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'], out_type=out_args['out_type'])
    pass_writer = TSVWriter(pass_handle, out_fields)

    # Count records
    result_count = countDbFile(db_file)

    # Iterate over records
    start_time = time()
    rec_count = 0
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time=start_time)
        rec_count += 1
        # TODO:  repeating renaming is unnecessary.
        # Rename fields
        for f, n in zip(fields, names):
            rec[n] = rec.pop(f)
        # Write
        pass_writer.writeDict(rec)

    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time=start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['END'] = 'ParseDb'
    printLog(log)

    # Close file handles
    pass_handle.close()
    db_handle.close()

    return pass_handle.name


def selectDbFile(db_file, fields, values, logic='any', regex=False,
                 out_file=None, out_args=default_out_args):
    """
    Selects records from a database file

    Arguments:
      db_file : the database file name
      fields : a list of fields to check for selection criteria
      values : a list of values defining selection targets
      logic : one of 'any' or 'all' defining whether one or all fields must have a match.
      regex : if False do exact full string matches; if True allow partial regex matches.
      out_file : output file name. Automatically generated from the input file if None.
      out_args : common output argument dictionary from parseCommonArgs

    Returns:
      str : output file name.
    """
    # Define string match function
    if regex:
        def _match_func(x, patterns):  return any([re.search(p, x) for p in patterns])
    else:
        def _match_func(x, patterns):  return x in patterns

    # Define logic function
    if logic == 'any':
        _logic_func = any
    elif logic == 'all':
        _logic_func = all

    # Print console log
    log = OrderedDict()
    log['START'] = 'ParseDb'
    log['COMMAND'] = 'select'
    log['FILE'] = os.path.basename(db_file)
    log['FIELDS'] = ','.join(fields)
    log['VALUES'] = ','.join(values)
    log['REGEX'] =regex
    printLog(log)

    # Open input
    db_handle = open(db_file, 'rt')
    db_iter = TSVReader(db_handle)
    out_fields = db_iter.fields
    __, __, out_args['out_type'] = splitName(db_file)

    # Open output
    if out_file is not None:
        pass_handle = open(out_file, 'w')
    else:
        pass_handle = getOutputHandle(db_file, out_label='parse-select', out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'], out_type=out_args['out_type'])
    pass_writer = TSVWriter(pass_handle, out_fields)

    # Count records
    result_count = countDbFile(db_file)

    # Iterate over records
    start_time = time()
    rec_count, pass_count, fail_count = 0, 0, 0
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time=start_time)
        rec_count += 1

        # Check for selection values in all fields
        select = _logic_func([_match_func(rec.get(f, False), values) for f in fields])

        # Write sequences
        if select:
            pass_count += 1
            pass_writer.writeDict(rec)
        else:
            fail_count += 1

    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time=start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['SELECTED'] = pass_count
    log['DISCARDED'] = fail_count
    log['END'] = 'ParseDb'
    printLog(log)

    # Close file handles
    pass_handle.close()
    db_handle.close()

    return pass_handle.name


def sortDbFile(db_file, field, numeric=False, descend=False,
               out_file=None, out_args=default_out_args):
    """
    Sorts records by values in an annotation field

    Arguments:
      db_file : the database filename
      field : the field name to sort by
      numeric : if True sort field numerically;
                if False sort field alphabetically
      descend : if True sort in descending order;
                if False sort in ascending order
      out_file : output file name. Automatically generated from the input file if None.
      out_args : common output argument dictionary from parseCommonArgs

    Returns:
      str : output file name
    """
    log = OrderedDict()
    log['START'] = 'ParseDb'
    log['COMMAND'] = 'sort'
    log['FILE'] = os.path.basename(db_file)
    log['FIELD'] = field
    log['NUMERIC'] = numeric
    printLog(log)

    # Open input
    db_handle = open(db_file, 'rt')
    db_iter = TSVReader(db_handle)
    out_fields = db_iter.fields
    __, __, out_args['out_type'] = splitName(db_file)

    # Open output
    if out_file is not None:
        pass_handle = open(out_file, 'w')
    else:
        pass_handle = getOutputHandle(db_file, out_label='parse-sort', out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'], out_type=out_args['out_type'])
    pass_writer = TSVWriter(pass_handle, out_fields)

    # Store all records in a dictionary
    start_time = time()
    printMessage("Indexing: Running", start_time=start_time)
    db_dict = {i:r for i, r in enumerate(db_iter)}
    result_count = len(db_dict)

    # Sort db_dict by field values
    tag_dict = {k:v[field] for k, v in db_dict.items()}
    if numeric:  tag_dict = {k:float(v or 0) for k, v in tag_dict.items()}
    sorted_keys = sorted(tag_dict, key=tag_dict.get, reverse=descend)
    printMessage("Indexing: Done", start_time=start_time, end=True)

    # Iterate over records
    start_time = time()
    rec_count = 0
    for key in sorted_keys:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time=start_time)
        rec_count += 1

        # Write records
        pass_writer.writeDict(db_dict[key])

    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time=start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['END'] = 'ParseDb'
    printLog(log)

    # Close file handles
    pass_handle.close()
    db_handle.close()

    return pass_handle.name


def updateDbFile(db_file, field, values, updates, out_file=None, out_args=default_out_args):
    """
    Updates field and value pairs to a database file

    Arguments:
      db_file : the database file name.
      field : the field to update.
      values : a list of values to specifying which rows to update.
      updates : a list of values to update each value with.
      out_file : output file name. Automatically generated from the input file if None.
      out_args : common output argument dictionary from parseCommonArgs.

    Returns:
      str : output file name
    """
    log = OrderedDict()
    log['START'] = 'ParseDb'
    log['COMMAND'] = 'update'
    log['FILE'] = os.path.basename(db_file)
    log['FIELD'] = field
    log['VALUES'] = ','.join(values)
    log['UPDATES'] = ','.join(updates)
    printLog(log)

    # Open input
    db_handle = open(db_file, 'rt')
    db_iter = TSVReader(db_handle)
    out_fields = db_iter.fields
    __, __, out_args['out_type'] = splitName(db_file)

    # Open output
    if out_file is not None:
        pass_handle = open(out_file, 'w')
    else:
        pass_handle = getOutputHandle(db_file, out_label='parse-update', out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'], out_type=out_args['out_type'])
    pass_writer = TSVWriter(pass_handle, out_fields)

    # Count records
    result_count = countDbFile(db_file)

    # Iterate over records
    start_time = time()
    rec_count, pass_count = 0, 0
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time=start_time)
        rec_count += 1

        # Updated values if found
        for x, y in zip(values, updates):
            if rec[field] == x:
                rec[field] = y
                pass_count += 1

        # Write records
        pass_writer.writeDict(rec)

    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time=start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['UPDATED'] = pass_count
    log['END'] = 'ParseDb'
    printLog(log)

    # Close file handles
    pass_handle.close()
    db_handle.close()

    return pass_handle.name


def mergeDbFiles(db_files, drop=False, out_file=None, out_args=default_out_args):
    """
    Updates field and value pairs to a database file

    Arguments:
      db_files : list of database file names.
      drop : if True drop columns not present in all files.
      out_file : output file name. Automatically generated from the input file if None.
      out_args : common output argument dictionary from parseCommonArgs.

    Returns:
      str : output file name.
    """
    log = OrderedDict()
    log['START'] = 'ParseDb'
    log['COMMAND'] = 'merge'
    log['FILES'] = ','.join([os.path.basename(f) for f in db_files])
    log['DROP'] = drop
    printLog(log)

    # Open input
    db_handles = [open(f, 'rt') for f in db_files]
    db_iters = [TSVReader(x) for x in db_handles]
    result_count = sum([countDbFile(f) for f in db_files])

    # Define output fields
    field_list = [x.fields for x in db_iters]
    if drop:
        field_set = set.intersection(*map(set, field_list))
    else:
        field_set = set.union(*map(set, field_list))
    field_order = OrderedDict([(f, None) for f in chain(*field_list)])
    out_fields = [f for f in field_order if f in field_set]

    # Open output file
    if out_file is not None:
        pass_handle = open(out_file, 'w')
    else:
        __, __, out_args['out_type'] = splitName(db_files[0])
        pass_handle = getOutputHandle(db_files[0], out_label='parse-merge', out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'], out_type=out_args['out_type'])
    pass_writer = TSVWriter(pass_handle, out_fields)

    # Iterate over records
    start_time = time()
    rec_count = 0
    for db in db_iters:
        for rec in db:
            # Print progress for previous iteration
            printProgress(rec_count, result_count, 0.05, start_time=start_time)
            rec_count += 1

            # Write records
            pass_writer.writeDict(rec)

    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time=start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['END'] = 'ParseDb'
    printLog(log)

    # Close file handles
    pass_handle.close()
    for x in db_handles: x.close()

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
    fields = dedent(
             '''
             output files:
                 sequences
                     FASTA formatted sequences output from the subcommands fasta and clip.
                 <field>-<value>
                     database files partitioned by annotation <field> and <value>.
                 parse-<command>
                     output of the database modification functions where <command> is one of
                     the subcommands add, index, drop, delete, rename, select, sort or update.

             required fields:
                 SEQUENCE_ID
                 
             optional fields:
                 JUNCTION, SEQUENCE_IMGT, SEQUENCE_VDJ, GERMLINE_IMGT, GERMLINE_VDJ,
                 GERMLINE_IMGT_D_MASK, GERMLINE_VDJ_D_MASK,
                 GERMLINE_IMGT_V_REGION, GERMLINE_VDJ_V_REGION
             ''')
    
    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            formatter_class=CommonHelpFormatter, add_help=False)
    group_help = parser.add_argument_group('help')
    group_help.add_argument('--version', action='version',
                            version='%(prog)s:' + ' %s-%s' %(__version__, __date__))
    group_help.add_argument('-h', '--help', action='help', help='show this help message and exit')
    subparsers = parser.add_subparsers(title='subcommands', dest='command', metavar='',
                                       help='Database operation')
    # TODO:  This is a temporary fix for Python issue 9253
    subparsers.required = True

    # Define parent parsers
    default_parent = getCommonArgParser(failed=False, log=False, format=False)
    multi_parent = getCommonArgParser(out_file=False, failed=False, log=False, format=False)

    # Subparser to add records
    parser_add = subparsers.add_parser('add', parents=[default_parent],
                                       formatter_class=CommonHelpFormatter, add_help=False,
                                       help='Adds field and value pairs.',
                                       description='Adds field and value pairs.')
    group_add = parser_add.add_argument_group('parsing arguments')
    group_add.add_argument('-f', nargs='+', action='store', dest='fields', required=True,
                           help='The name of the fields to add.')
    group_add.add_argument('-u', nargs='+', action='store', dest='values', required=True,
                           help='The value to assign to all rows for each field.')
    parser_add.set_defaults(func=addDbFile)

    # Subparser to delete records
    parser_delete = subparsers.add_parser('delete', parents=[default_parent],
                                          formatter_class=CommonHelpFormatter, add_help=False,
                                          help='Deletes specific records.',
                                          description='Deletes specific records.')
    group_delete = parser_delete.add_argument_group('parsing arguments')
    group_delete.add_argument('-f', nargs='+', action='store', dest='fields', required=True,
                               help='The name of the fields to check for deletion criteria.')
    group_delete.add_argument('-u', nargs='+', action='store', dest='values', default=['', 'NA'],
                               help='''The values defining which records to delete. A value
                                    may appear in any of the fields specified with -f.''')
    group_delete.add_argument('--logic', action='store', dest='logic',
                               choices=('any', 'all'), default='any',
                               help='''Defines whether a value may appear in any field (any)
                                    or whether it must appear in all fields (all).''')
    group_delete.add_argument('--regex', action='store_true', dest='regex',
                               help='''If specified, treat values as regular expressions
                                    and allow partial string matches.''')
    parser_delete.set_defaults(func=deleteDbFile)

    # Subparser to drop fields
    parser_drop = subparsers.add_parser('drop', parents=[default_parent],
                                        formatter_class=CommonHelpFormatter, add_help=False,
                                        help='Deletes entire fields.',
                                        description='Deletes entire fields.')
    group_drop = parser_drop.add_argument_group('parsing arguments')
    group_drop.add_argument('-f', nargs='+', action='store', dest='fields', required=True,
                               help='The name of the fields to delete from the database.')
    parser_drop.set_defaults(func=dropDbFile)

    # Subparser to index fields
    parser_index = subparsers.add_parser('index', parents=[default_parent],
                                         formatter_class=CommonHelpFormatter, add_help=False,
                                         help='Adds a numeric index field.',
                                         description='Adds a numeric index field.')
    group_index = parser_index.add_argument_group('parsing arguments')
    group_index.add_argument('-f', action='store', dest='field',
                              default=default_index_field,
                              help='The name of the index field to add to the database.')
    parser_index.set_defaults(func=indexDbFile)

    # Subparser to rename fields
    parser_rename = subparsers.add_parser('rename', parents=[default_parent],
                                          formatter_class=CommonHelpFormatter, add_help=False,
                                          help='Renames fields.',
                                          description='Renames fields.')
    group_rename = parser_rename.add_argument_group('parsing arguments')
    group_rename.add_argument('-f', nargs='+', action='store', dest='fields', required=True,
                               help='List of fields to rename.')
    group_rename.add_argument('-k', nargs='+', action='store', dest='names', required=True,
                               help='List of new names for each field.')
    parser_rename.set_defaults(func=renameDbFile)

    # Subparser to select records
    parser_select = subparsers.add_parser('select', parents=[default_parent],
                                          formatter_class=CommonHelpFormatter, add_help=False,
                                          help='Selects specific records.',
                                          description='Selects specific records.')
    group_select = parser_select.add_argument_group('parsing arguments')
    group_select.add_argument('-f', nargs='+', action='store', dest='fields', required=True,
                               help='The name of the fields to check for selection criteria.')
    group_select.add_argument('-u', nargs='+', action='store', dest='values', required=True,
                               help='''The values defining with records to select. A value
                                    may appear in any of the fields specified with -f.''')
    group_select.add_argument('--logic', action='store', dest='logic',
                               choices=('any', 'all'), default='any',
                               help='''Defines whether a value may appear in any field (any)
                                    or whether it must appear in all fields (all).''')
    group_select.add_argument('--regex', action='store_true', dest='regex',
                               help='''If specified, treat values as regular expressions
                                    and allow partial string matches.''')
    parser_select.set_defaults(func=selectDbFile)

    # Subparser to sort file by records
    parser_sort = subparsers.add_parser('sort', parents=[default_parent],
                                        formatter_class=CommonHelpFormatter, add_help=False,
                                        help='Sorts records by field values.',
                                        description='Sorts records by field values.')
    group_sort = parser_sort.add_argument_group('parsing arguments')
    group_sort.add_argument('-f', action='store', dest='field', type=str, required=True,
                             help='The annotation field by which to sort records.')
    group_sort.add_argument('--num', action='store_true', dest='numeric', default=False,
                             help='''Specify to define the sort column as numeric rather
                                  than textual.''')
    group_sort.add_argument('--descend', action='store_true', dest='descend',
                             help='''If specified, sort records in descending, rather
                             than ascending, order by values in the target field.''')
    parser_sort.set_defaults(func=sortDbFile)

    # Subparser to update records
    parser_update = subparsers.add_parser('update', parents=[default_parent],
                                          formatter_class=CommonHelpFormatter, add_help=False,
                                          help='Updates field and value pairs.',
                                          description='Updates field and value pairs.')
    group_update = parser_update.add_argument_group('parsing arguments')
    group_update.add_argument('-f', action='store', dest='field', required=True,
                               help='The name of the field to update.')
    group_update.add_argument('-u', nargs='+', action='store', dest='values', required=True,
                               help='The values that will be replaced.')
    group_update.add_argument('-t', nargs='+', action='store', dest='updates', required=True,
                               help='''The new value to assign to each selected row.''')
    parser_update.set_defaults(func=updateDbFile)

    # Subparser to merge files
    parser_merge = subparsers.add_parser('merge', parents=[multi_parent],
                                         formatter_class=CommonHelpFormatter, add_help=False,
                                         help='Merges files.',
                                         description='Merges files.')
    group_merge = parser_merge.add_argument_group('parsing arguments')
    group_merge.add_argument('-o', action='store', dest='out_file', default=None,
                              help='''Explicit output file name. Note, this argument cannot be used with 
                                   the --failed, --outdir or --outname arguments.''')
    group_merge.add_argument('--drop', action='store_true', dest='drop',
                              help='''If specified, drop fields that do not exist in all input files.
                                   Otherwise, include all columns in all files and fill missing data 
                                   with empty strings.''')
    parser_merge.set_defaults(func=mergeDbFiles)

    # Subparser to partition files by annotation values
    parser_split = subparsers.add_parser('split', parents=[multi_parent],
                                         formatter_class=CommonHelpFormatter, add_help=False,
                                         help='Splits database files by field values.',
                                         description='Splits database files by field values')
    group_split = parser_split.add_argument_group('parsing arguments')
    group_split.add_argument('-f', action='store', dest='field', type=str, required=True,
                              help='Annotation field by which to split database files.')
    group_split.add_argument('--num', action='store', dest='num_split', type=float, default=None,
                              help='''Specify to define the field as numeric and group
                                   records by whether they are less than or at least
                                   (greater than or equal to) the specified value.''')
    parser_split.set_defaults(func=splitDbFile)

    return parser


if __name__ == '__main__':
    """
    Parses command line arguments and calls main function
    """
    # Parse arguments
    parser = getArgParser()
    checkArgs(parser)
    args = parser.parse_args()
    if args.command == 'merge':
        args_dict = parseCommonArgs(args, in_list=True)
    else:
        args_dict = parseCommonArgs(args)
    # Delete command declaration from argument dictionary
    del args_dict['command']
    del args_dict['func']

    # Check argument pairs
    if args.command == 'add' and len(args_dict['fields']) != len(args_dict['values']):
        parser.error('You must specify exactly one value (-u) per field (-f)')
    elif args.command == 'rename' and len(args_dict['fields']) != len(args_dict['names']):
        parser.error('You must specify exactly one new name (-k) per field (-f)')
    elif args.command == 'update' and len(args_dict['values']) != len(args_dict['updates']):
        parser.error('You must specify exactly one value (-u) per replacement (-t)')

    # Call parser function for each database file
    if args.command == 'merge':
        args.func(**args_dict)
    elif args.command == 'split':
        del args_dict['db_files']
        for f in args.__dict__['db_files']:
            args_dict['db_file'] = f
            args.func(**args_dict)
    else:
        del args_dict['db_files']
        del args_dict['out_files']
        for i, f in enumerate(args.__dict__['db_files']):
            args_dict['db_file'] = f
            args_dict['out_file'] = args.__dict__['out_files'][i] \
                if args.__dict__['out_files'] else None
            args.func(**args_dict)
 
