#!/usr/bin/env python3
"""
Multiple aligns sequence fields
"""

# Info
__author__ = 'Jason Anthony Vander Heiden'
from changeo import __version__, __date__

# Imports
import os
import shutil
from argparse import ArgumentParser
from collections import OrderedDict
from itertools import chain
from textwrap import dedent
from Bio.SeqRecord import SeqRecord

# Presto and changeo import
from presto.Defaults import default_out_args, default_muscle_exec
from presto.Applications import runMuscle
from presto.IO import printLog, printError, printWarning
from presto.Multiprocessing import manageProcesses
from changeo.Commandline import CommonHelpFormatter, checkArgs, getCommonArgParser, parseCommonArgs
from changeo.IO import getDbFields, getFormatOperators
from changeo.Multiprocessing import DbResult, feedDbQueue, processDbQueue, collectDbQueue


# TODO:  maybe not bothering with 'set' is best. can just work off field identity
def groupRecords(records, fields=None, calls=['v', 'j'], mode='gene', action='first'):
    """
    Groups Receptor objects based on gene or annotation

    Arguments:
      records : an iterator of Receptor objects to group.
      fields : gene field to group by.
      calls : allele calls to use for grouping.
              one or more of ('v', 'd', 'j').
      mode : specificity of alignment call to use for allele call fields.
             one of ('allele', 'gene').
      action : only 'first' is currently supported.

    Returns:
    dictionary of grouped records
    """
    # Define functions for grouping keys
    if mode == 'allele' and fields is None:
        def _get_key(rec, calls, action):
            return tuple(rec.getAlleleCalls(calls, action))
    elif mode == 'gene' and fields is None:
        def _get_key(rec, calls, action):
            return tuple(rec.getGeneCalls(calls, action))
    elif mode == 'allele' and fields is not None:
        def _get_key(rec, calls, action):
            vdj = rec.getAlleleCalls(calls, action)
            ann = [rec.getChangeo(k) for k in fields]
            return tuple(chain(vdj, ann))
    elif mode == 'gene' and fields is not None:
        def _get_key(rec, calls, action):
            vdj = rec.getGeneCalls(calls, action)
            ann = [rec.getChangeo(k) for k in fields]
            return tuple(chain(vdj, ann))

    rec_index = {}
    for rec in records:
        key = _get_key(rec, calls, action)
        # Assigned grouped records to individual keys and all failed to a single key
        if all([k is not None for k in key]):
            rec_index.setdefault(key, []).append(rec)
        else:
            rec_index.setdefault(None, []).append(rec)

    return rec_index


def alignBlocks(data, field_map, muscle_exec=default_muscle_exec):
    """
    Multiple aligns blocks of sequence fields together

    Arguments:
      data : DbData object with Receptor objects to process.
      field_map : a dictionary of {input sequence : output sequence) field names to multiple align.
      muscle_exec : the MUSCLE executable.

    Returns:
      changeo.Multiprocessing.DbResult : object containing Receptor objects with multiple aligned sequence fields.
    """
    # Define return object
    result = DbResult(data.id, data.data)
    result.results = data.data
    result.valid = True

    # Fail invalid groups
    if result.id is None:
        result.log = None
        result.valid = False
        return result

    seq_fields = list(field_map.keys())
    seq_list = [SeqRecord(r.getSeq(f), id='%s_%s' % (r.sequence_id.replace(' ', '_'), f)) for f in seq_fields \
                for r in data.data]
    seq_aln = runMuscle(seq_list, aligner_exec=muscle_exec)
    if seq_aln is not None:
        aln_map = {x.id: i for i, x in enumerate(seq_aln)}
        for i, r in enumerate(result.results, start=1):
            for f in seq_fields:
                idx = aln_map['%s_%s' % (r.sequence_id.replace(' ', '_'), f)]
                seq = str(seq_aln[idx].seq)
                r.annotations[field_map[f]] = seq
                result.log['%s-%s' % (f, r.sequence_id)] = seq

    else:
        result.valid = False

    #for r in result.results:  print r.annotations
    return result


def alignAcross(data, field_map, muscle_exec=default_muscle_exec):
    """
    Multiple aligns sequence fields column wise

    Arguments:
      data : DbData object with Receptor objects to process.
      field_map : a dictionary of {input sequence : output sequence) field names to multiple align.
      muscle_exec : the MUSCLE executable.

    Returns:
      changeo.Multiprocessing.DbResult : object containing Receptor objects with multiple aligned sequence fields.
    """
    # Define return object
    result = DbResult(data.id, data.data)
    result.results = data.data
    result.valid = True

    # Fail invalid groups
    if result.id is None:
        result.log = None
        result.valid = False
        return result

    seq_fields = list(field_map.keys())
    for f in seq_fields:
        seq_list = [SeqRecord(r.getSeq(f), id=r.sequence_id.replace(' ', '_')) for r in data.data]
        seq_aln = runMuscle(seq_list, aligner_exec=muscle_exec)
        if seq_aln is not None:
            aln_map = {x.id: i for i, x in enumerate(seq_aln)}
            for i, r in enumerate(result.results, start=1):
                idx = aln_map[r.sequence_id.replace(' ', '_')]
                seq = str(seq_aln[idx].seq)
                r.annotations[field_map[f]] = seq
                result.log['%s-%s' % (f, r.sequence_id)] = seq
        else:
            result.valid = False

    #for r in result.results:  print r.annotations
    return result


def alignWithin(data, field_map, muscle_exec=default_muscle_exec):
    """
    Multiple aligns sequence fields within a row

    Arguments:
      data : DbData object with Receptor objects to process.
      field_map : a dictionary of {input sequence : output sequence) field names to multiple align.
      muscle_exec : the MUSCLE executable.

    Returns:
      changeo.Multiprocessing.DbResult : object containing Receptor objects with multiple aligned sequence fields.
    """
    # Define return object
    result = DbResult(data.id, data.data)
    result.results = data.data
    result.valid = True

    # Fail invalid groups
    if result.id is None:
        result.log = None
        result.valid = False
        return result

    record = data.data
    seq_fields = list(field_map.keys())
    seq_list = [SeqRecord(record.getSeq(f), id=f) for f in seq_fields]
    seq_aln = runMuscle(seq_list, aligner_exec=muscle_exec)
    if seq_aln is not None:
        aln_map = {x.id: i for i, x in enumerate(seq_aln)}
        for f in seq_fields:
            idx = aln_map[f]
            seq = str(seq_aln[idx].seq)
            record.annotations[field_map[f]] = seq
            result.log[f] = seq
    else:
        result.valid = False

    return result


def alignRecords(db_file, seq_fields, group_func, align_func, group_args={}, align_args={},
                 format='changeo', out_file=None, out_args=default_out_args, nproc=None, queue_size=None):
    """
    Performs a multiple alignment on sets of sequences

    Arguments: 
      db_file : filename of the input database.
      seq_fields : the sequence fields to multiple align.
      group_func : function to use to group records.
      align_func : function to use to multiple align sequence groups.
      group_args : dictionary of arguments to pass to group_func.
      align_args : dictionary of arguments to pass to align_func.
      format : output format. One of 'changeo' or 'airr'.
      out_file : output file name. Automatically generated from the input file if None.
      out_args : common output argument dictionary from parseCommonArgs.
      nproc : the number of processQueue processes.
              if None defaults to the number of CPUs.
      queue_size : maximum size of the argument queue.
                   if None defaults to 2*nproc.
                      
    Returns: 
      dict : names of the 'pass' and 'fail' output files.
    """
    # Define subcommand label dictionary
    cmd_dict = {alignAcross: 'across', alignWithin: 'within', alignBlocks: 'block'}
    
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'AlignRecords'
    log['COMMAND'] = cmd_dict.get(align_func, align_func.__name__)
    log['FILE'] = os.path.basename(db_file)
    log['SEQ_FIELDS'] = ','.join(seq_fields)
    if 'group_fields' in group_args: log['GROUP_FIELDS'] = ','.join(group_args['group_fields'])
    if 'mode' in group_args: log['MODE'] = group_args['mode']
    if 'action' in group_args: log['ACTION'] = group_args['action']
    log['NPROC'] = nproc
    printLog(log)

    # Define format operators
    try:
        reader, writer, schema = getFormatOperators(format)
    except ValueError:
        printError('Invalid format %s.' % format)

    # Define feeder function and arguments
    if 'group_fields' in group_args and group_args['group_fields'] is not None:
        group_args['group_fields'] = [schema.toReceptor(f) for f in group_args['group_fields']]
    feed_func = feedDbQueue
    feed_args = {'db_file': db_file,
                 'reader': reader,
                 'group_func': group_func,
                 'group_args': group_args}
    # Define worker function and arguments
    field_map = OrderedDict([(schema.toReceptor(f), '%s_align' % f) for f in seq_fields])
    align_args['field_map'] = field_map
    work_func = processDbQueue
    work_args = {'process_func': align_func,
                 'process_args': align_args}
    # Define collector function and arguments
    out_fields = getDbFields(db_file, add=list(field_map.values()), reader=reader)
    out_args['out_type'] = schema.out_type
    collect_func = collectDbQueue
    collect_args = {'db_file': db_file,
                    'label': 'align',
                    'fields': out_fields,
                    'writer': writer,
                    'out_file': out_file,
                    'out_args': out_args}
    
    # Call process manager
    result = manageProcesses(feed_func, work_func, collect_func, 
                             feed_args, work_args, collect_args, 
                             nproc, queue_size)
        
    # Print log
    result['log']['END'] = 'AlignRecords'
    printLog(result['log'])
    output = {k: v for k, v in result.items() if k in ('pass', 'fail')}

    return output


def getArgParser():
    """
    Defines the ArgumentParser

    Arguments: 
    None
                      
    Returns: 
    an ArgumentParser object
    """
    # Define output file names and header fields
    fields = dedent(
             '''
             output files:
                 align-pass
                     database with multiple aligned sequences.
                 align-fail
                     database with records failing alignment.

             required fields:
                 SEQUENCE_ID, V_CALL, J_CALL
                 <field>
                     user specified sequence fields to align.

             output fields:
                 <field>_ALIGN
             ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            formatter_class=CommonHelpFormatter, add_help=False)
    group_help = parser.add_argument_group('help')
    group_help.add_argument('--version', action='version',
                            version='%(prog)s:' + ' %s-%s' %(__version__, __date__))
    group_help.add_argument('-h', '--help', action='help', help='show this help message and exit')
    subparsers = parser.add_subparsers(title='subcommands', dest='command', metavar='',
                                       help='alignment method')
    # TODO:  This is a temporary fix for Python issue 9253
    subparsers.required = True

    # Parent parser
    parser_parent = getCommonArgParser(format=False, multiproc=True)

    # Argument parser for column-wise alignment across records
    parser_across = subparsers.add_parser('across', parents=[parser_parent],
                                           formatter_class=CommonHelpFormatter, add_help=False,
                                           help='''Multiple aligns sequence columns within groups 
                                                 and across rows using MUSCLE.''')
    group_across = parser_across.add_argument_group('alignment arguments')
    group_across.add_argument('--sf', nargs='+', action='store', dest='seq_fields', required=True,
                               help='The sequence fields to multiple align within each group.')
    group_across.add_argument('--gf', nargs='+', action='store', dest='group_fields', default=None,
                               help='Additional (not allele call) fields to use for grouping.')
    group_across.add_argument('--calls', nargs='+', action='store', dest='calls',
                               choices=('v', 'd', 'j'), default=['v', 'j'],
                               help='Segment calls (allele assignments) to use for grouping.')
    group_across.add_argument('--mode', action='store', dest='mode',
                              choices=('allele', 'gene'), default='gene',
                              help='''Specifies whether to use the V(D)J allele or gene when
                                   an allele call field (--calls) is specified.''')
    group_across.add_argument('--act', action='store', dest='action', default='first',
                               choices=('first', ),
                               help='''Specifies how to handle multiple values within default
                                     allele call fields. Currently, only "first" is supported.''')
    group_across.add_argument('--exec', action='store', dest='muscle_exec',
                               default=default_muscle_exec,
                               help='The location of the MUSCLE executable')
    parser_across.set_defaults(group_func=groupRecords, align_func=alignAcross)


    # Argument parser for alignment of fields within records
    parser_within = subparsers.add_parser('within', parents=[parser_parent],
                                          formatter_class=CommonHelpFormatter, add_help=False,
                                          help='Multiple aligns sequence fields within rows using MUSCLE')
    group_within = parser_within.add_argument_group('alignment arguments')
    group_within.add_argument('--sf', nargs='+', action='store', dest='seq_fields', required=True,
                               help='The sequence fields to multiple align within each record.')
    group_within.add_argument('--exec', action='store', dest='muscle_exec',
                              default=default_muscle_exec,
                              help='The location of the MUSCLE executable')
    parser_within.set_defaults(group_func=None, align_func=alignWithin)

    # Argument parser for column-wise alignment across records
    parser_block = subparsers.add_parser('block', parents=[parser_parent],
                                        formatter_class=CommonHelpFormatter, add_help=False,
                                        help='''Multiple aligns sequence groups across both 
                                             columns and rows using MUSCLE.''')
    group_block = parser_block.add_argument_group('alignment arguments')
    group_block.add_argument('--sf', nargs='+', action='store', dest='seq_fields', required=True,
                               help='The sequence fields to multiple align within each group.')
    group_block.add_argument('--gf', nargs='+', action='store', dest='group_fields', default=None,
                               help='Additional (not allele call) fields to use for grouping.')
    group_block.add_argument('--calls', nargs='+', action='store', dest='calls',
                               choices=('v', 'd', 'j'), default=['v', 'j'],
                               help='Segment calls (allele assignments) to use for grouping.')
    group_block.add_argument('--mode', action='store', dest='mode',
                              choices=('allele', 'gene'), default='gene',
                              help='''Specifies whether to use the V(D)J allele or gene when
                                   an allele call field (--calls) is specified.''')
    group_block.add_argument('--act', action='store', dest='action', default='first',
                               choices=('first', ),
                               help='''Specifies how to handle multiple values within default
                                     allele call fields. Currently, only "first" is supported.''')
    group_block.add_argument('--exec', action='store', dest='muscle_exec',
                               default=default_muscle_exec,
                               help='The location of the MUSCLE executable')
    parser_block.set_defaults(group_func=groupRecords, align_func=alignBlocks)

    return parser


if __name__ == '__main__':
    """
    Parses command line arguments and calls main function
    """
    # Parse arguments
    parser = getArgParser()
    checkArgs(parser)
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)

    # Convert case of fields
    # if 'seq_fields' in args_dict:
    #     args_dict['seq_fields'] = [f.upper() for f in args_dict['seq_fields']]
    # if 'group_fields' in args_dict and args_dict['group_fields'] is not None:
    #     args_dict['group_fields'] = [f.upper() for f in args_dict['group_fields']]

    # Check if a valid MUSCLE executable was specified for muscle mode
    if not shutil.which(args.muscle_exec):
        parser.error('%s does not exist or is not executable.' % args.muscle_exec)
    # Define align_args
    args_dict['align_args'] = {'muscle_exec': args_dict['muscle_exec']}
    del args_dict['muscle_exec']

    # Define group_args
    if args_dict['group_func'] is groupRecords:
        args_dict['group_args'] = {'fields':args_dict['group_fields'],
                                   'calls':args_dict['calls'],
                                   'mode':args_dict['mode'],
                                   'action':args_dict['action']}
        del args_dict['group_fields']
        del args_dict['calls']
        del args_dict['mode']
        del args_dict['action']

    # Clean arguments dictionary
    del args_dict['command']
    del args_dict['db_files']
    if 'out_files' in args_dict: del args_dict['out_files']

    # Call main function for each input file
    for i, f in enumerate(args.__dict__['db_files']):
        args_dict['db_file'] = f
        args_dict['out_file'] = args.__dict__['out_files'][i] \
            if args.__dict__['out_files'] else None
        alignRecords(**args_dict)

