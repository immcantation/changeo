#!/usr/bin/env python3
"""
Reconstructs germline sequences from alignment data
"""

# Info
__author__ = 'Namita Gupta, Jason Anthony Vander Heiden'
from changeo import __version__, __date__

# Imports
import os
from argparse import ArgumentParser
from collections import OrderedDict
from itertools import groupby
from textwrap import dedent
from time import time

# Presto and change imports
from presto.Defaults import default_out_args
from presto.IO import printLog, printMessage, printProgress, printError, printWarning
from changeo.Defaults import default_v_field, default_d_field, default_j_field, default_clone_field, \
                             default_seq_field, default_format
from changeo.Commandline import CommonHelpFormatter, checkArgs, getCommonArgParser, parseCommonArgs
from changeo.Gene import buildGermline, buildClonalGermline
from changeo.IO import countDbFile, getDbFields, getFormatOperators, getOutputHandle, readGermlines, \
                       checkFields

# Defaults
default_germ_types = 'dmask'


def createGermlines(db_file, references, seq_field=default_seq_field, v_field=default_v_field,
                    d_field=default_d_field, j_field=default_j_field,
                    cloned=False, clone_field=default_clone_field, germ_types=default_germ_types,
                    format=default_format, out_file=None, out_args=default_out_args):
    """
    Write germline sequences to tab-delimited database file

    Arguments:
      db_file : input tab-delimited database file.
      references : folders and/or files containing germline repertoire data in FASTA format.
      seq_field : field in which to look for sequence.
      v_field : field in which to look for V call.
      d_field : field in which to look for D call.
      j_field : field in which to look for J call.
      cloned : if True build germlines by clone, otherwise build individual germlines.
      clone_field : field containing clone identifiers; ignored if cloned=False.
      germ_types : list of germline sequence types to be output from the set of 'full', 'dmask', 'vonly', 'regions'
      format : input and output format.
      out_file : output file name. Automatically generated from the input file if None.
      out_args : arguments for output preferences.

    Returns:
      dict: names of the 'pass' and 'fail' output files.
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'CreateGermlines'
    log['FILE'] = os.path.basename(db_file)
    log['GERM_TYPES'] = ','.join(germ_types)
    log['SEQ_FIELD'] = seq_field
    log['V_FIELD'] = v_field
    log['D_FIELD'] = d_field
    log['J_FIELD'] = j_field
    log['CLONED'] = cloned
    if cloned:  log['CLONE_FIELD'] = clone_field
    printLog(log)

    # Define format operators
    try:
        reader, writer, schema = getFormatOperators(format)
    except ValueError:
        printError('Invalid format %s' % format)
    out_args['out_type'] = schema.out_type

    # TODO: this won't work for AIRR necessarily
    # Define output germline fields
    germline_fields = OrderedDict()
    seq_type = seq_field.split('_')[-1]
    if 'full' in germ_types:  germline_fields['full'] = 'germline_' + seq_type
    if 'dmask' in germ_types:  germline_fields['dmask'] = 'germline_' + seq_type + '_d_mask'
    if 'vonly' in germ_types:  germline_fields['vonly'] = 'germline_' + seq_type + '_v_region'
    if 'regions' in germ_types:  germline_fields['regions'] = 'germline_regions'
    if cloned:
        germline_fields['v'] = 'germline_v_call'
        germline_fields['d'] = 'germline_d_call'
        germline_fields['j'] = 'germline_j_call'
    out_fields = getDbFields(db_file,
                             add=[schema.fromReceptor(f) for f in germline_fields.values()],
                             reader=reader)

    # Get repertoire and open Db reader
    reference_dict = readGermlines(references)
    db_handle = open(db_file, 'rt')
    db_iter = reader(db_handle)

    # Check for required columns
    try:
        required = ['v_germ_start_imgt', 'd_germ_start', 'j_germ_start',
                    'np1_length', 'np2_length']
        checkFields(required, db_iter.fields, schema=schema)
    except LookupError as e:
        printError(e)

    # Check for IMGT-gaps in germlines
    if all('...' not in x for x in reference_dict.values()):
        printWarning('Germline reference sequences do not appear to contain IMGT-numbering spacers. Results may be incorrect.')

    # Count input
    total_count = countDbFile(db_file)

    # Check for existence of fields
    for f in [v_field, d_field, j_field, seq_field]:
        if f not in db_iter.fields:
            printError('%s field does not exist in input database file.' % f)

    # Translate to Receptor attribute names
    v_field = schema.toReceptor(v_field)
    d_field = schema.toReceptor(d_field)
    j_field = schema.toReceptor(j_field)
    seq_field = schema.toReceptor(seq_field)
    clone_field = schema.toReceptor(clone_field)

    # Define Receptor iterator
    if cloned:
        start_time = time()
        printMessage('Sorting by clone', start_time=start_time, width=20)
        sorted_records = sorted(db_iter, key=lambda x: x.getField(clone_field))
        printMessage('Done', start_time=start_time, end=True, width=20)
        receptor_iter = groupby(sorted_records, lambda x: x.getField(clone_field))
    else:
        receptor_iter = ((x.sequence_id, [x]) for x in db_iter)

    # Define log handle
    if out_args['log_file'] is None:
        log_handle = None
    else:
        log_handle = open(out_args['log_file'], 'w')

    # Initialize handles, writers and counters
    pass_handle, pass_writer = None, None
    fail_handle, fail_writer = None, None
    rec_count, pass_count, fail_count = 0, 0, 0
    start_time = time()

    # Iterate over rows
    for key, records in receptor_iter:
        # Print progress
        printProgress(rec_count, total_count, 0.05, start_time=start_time)

        # Define iteration variables
        records = list(records)
        rec_log = OrderedDict([('ID', key)])
        rec_count += len(records)

        # Build germline for records
        if len(records) == 1:
            germ_log, germlines, genes = buildGermline(records[0], reference_dict, seq_field=seq_field, v_field=v_field,
                                                       d_field=d_field, j_field=j_field)
        else:
            germ_log, germlines, genes = buildClonalGermline(records, reference_dict, seq_field=seq_field, v_field=v_field,
                                                             d_field=d_field, j_field=j_field)
        rec_log.update(germ_log)

        # Write row to pass or fail file
        if germlines is not None:
            pass_count += len(records)

            # Add germlines to Receptor record
            annotations = {}
            if 'full' in germ_types:  annotations[germline_fields['full']] = germlines['full']
            if 'dmask' in germ_types:  annotations[germline_fields['dmask']] = germlines['dmask']
            if 'vonly' in germ_types:  annotations[germline_fields['vonly']] = germlines['vonly']
            if 'regions' in germ_types:  annotations[germline_fields['regions']] = germlines['regions']
            if cloned:
                annotations[germline_fields['v']] = genes['v']
                annotations[germline_fields['d']] = genes['d']
                annotations[germline_fields['j']] = genes['j']

            # Write records
            try:
                for r in records:
                    r.setDict(annotations)
                    pass_writer.writeReceptor(r)
            except AttributeError:
                # Create output file handle and writer
                if out_file is not None:
                    pass_handle = open(out_file, 'w')
                else:
                    pass_handle = getOutputHandle(db_file,
                                                  out_label='germ-pass',
                                                  out_dir=out_args['out_dir'],
                                                  out_name=out_args['out_name'],
                                                  out_type=out_args['out_type'])
                pass_writer = writer(pass_handle, fields=out_fields)
                for r in records:
                    r.setDict(annotations)
                    pass_writer.writeReceptor(r)
        else:
            fail_count += len(records)
            if out_args['failed']:
                try:
                    fail_writer.writeReceptor(records)
                except AttributeError:
                    fail_handle = getOutputHandle(db_file,
                                                  out_label='germ-fail',
                                                  out_dir=out_args['out_dir'],
                                                  out_name=out_args['out_name'],
                                                  out_type=out_args['out_type'])
                    fail_writer = writer(fail_handle, fields=out_fields)
                    fail_writer.writeReceptor(records)

        # Write log
        printLog(rec_log, handle=log_handle)

    # Print log
    printProgress(rec_count, total_count, 0.05, start_time=start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name) if pass_handle is not None else None
    log['RECORDS'] = rec_count
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    log['END'] = 'CreateGermlines'
    printLog(log)

    # Close file handles
    db_handle.close()
    output = {'pass': None, 'fail': None}
    if pass_handle is not None:
        output['pass'] = pass_handle.name
        pass_handle.close()
    if fail_handle is not None:
        output['fail'] = fail_handle.name
        fail_handle.close()
    if log_handle is not None:
        log_handle.close()

    return output


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
                 germ-pass
                    database with assigned germline sequences.
                 germ-fail
                    database with records failing germline assignment.

             required fields:
                 SEQUENCE_ID, SEQUENCE_IMGT, V_CALL, D_CALL, J_CALL,
                 V_SEQ_START, V_SEQ_LENGTH, V_GERM_START_IMGT, V_GERM_LENGTH_IMGT,
                 D_SEQ_START, D_SEQ_LENGTH, D_GERM_START, D_GERM_LENGTH,
                 J_SEQ_START, J_SEQ_LENGTH, J_GERM_START, J_GERM_LENGTH,
                 NP1_LENGTH, NP2_LENGTH

             optional fields:
                 N1_LENGTH, N2_LENGTH, P3V_LENGTH, P5D_LENGTH, P3D_LENGTH, P5J_LENGTH,
                 CLONE


             output fields:
                 GERMLINE_IMGT, GERMLINE_IMGT_D_MASK, GERMLINE_IMGT_V_REGION,
                 GERMLINE_V_CALL, GERMLINE_D_CALL, GERMLINE_J_CALL,
                 GERMLINE_REGIONS
              ''')
    # Define argument parser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            parents=[getCommonArgParser(format=False)],
                            formatter_class=CommonHelpFormatter, add_help=False)

    # Germlines arguments
    group = parser.add_argument_group('germline construction arguments')
    group.add_argument('-r', nargs='+', action='store', dest='references', required=True,
                        help='''List of folders and/or fasta files (with .fasta, .fna or .fa
                         extension) with germline sequences. When using the default
                         Change-O sequence and coordinate fields, these reference sequences 
                         must contain IMGT-numbering spacers (gaps) in the V segment. 
                         Alternative numbering schemes, or no numbering, may work for alternative 
                         sequence and coordinate definitions that define a valid alignment, but 
                         a warning will be issued.''')
    group.add_argument('-g', action='store', dest='germ_types', default=default_germ_types,
                        nargs='+', choices=('full', 'dmask', 'vonly', 'regions'),
                        help='''Specify type(s) of germlines to include full germline,
                             germline with D segment masked, or germline for V segment only.''')
    group.add_argument('--cloned', action='store_true', dest='cloned',
                        help='''Specify to create only one germline per clone. Note, if allele 
                             calls are ambiguous within a clonal group, this will place the germline call 
                             used for the entire clone within the
                             GERMLINE_V_CALL, GERMLINE_D_CALL and GERMLINE_J_CALL fields.''')
    group.add_argument('--sf', action='store', dest='seq_field', default=default_seq_field,
                        help='Field containing the aligned sequence.')
    group.add_argument('--vf', action='store', dest='v_field', default=default_v_field,
                        help='Field containing the germline V segment call.')
    group.add_argument('--df', action='store', dest='d_field', default=default_d_field,
                        help='Field containing the germline D segment call.')
    group.add_argument('--jf', action='store', dest='j_field', default=default_j_field,
                        help='Field containing the germline J segment call.')

    return parser


if __name__ == '__main__':
    """
    Parses command line arguments and calls main
    """

    # Parse command line arguments
    parser = getArgParser()
    checkArgs(parser)
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)

    # # Set default fields if not specified.
    # default_fields = {'seq_field': default_seq_field,
    #                   'v_field': default_v_field,
    #                   'd_field': default_d_field,
    #                   'j_field': default_j_field}
    #
    # # Default Change-O fields
    # if args_dict['format'] == 'changeo':
    #     for f in default_fields:
    #         if args_dict[f] is None:  args_dict[f] = default_fields[f]
    #         else: args_dict[f] = args_dict[f].upper()
    #
    # # Default AIRR fields
    # if args_dict['format'] == 'airr':
    #     for f in default_fields:
    #         if args_dict[f] is None:  args_dict[f] = ChangeoSchema.fromReceptor(default_fields[f])
    #         else: args_dict[f] = args_dict[f].lower()

    # Check that reference files exist
    for f in args_dict['references']:
        if not os.path.exists(f):
            parser.error('Germline reference file or folder %s does not exist.' % f)

    # Clean arguments dictionary
    del args_dict['db_files']
    if 'out_files' in args_dict: del args_dict['out_files']

    # Call main function for each input file
    for i, f in enumerate(args.__dict__['db_files']):
        args_dict['db_file'] = f
        args_dict['out_file'] = args.__dict__['out_files'][i] \
            if args.__dict__['out_files'] else None
        createGermlines(**args_dict)
