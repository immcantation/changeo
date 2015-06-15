#!/usr/bin/env python
"""
Multiple aligns sequence fields
"""

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2013 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.2.0'
__date__      = '2015.05.30'

# Imports
import csv, os, sys, textwrap
from argparse import ArgumentParser
from collections import deque, OrderedDict
from cStringIO import StringIO
from itertools import chain, izip
from subprocess import PIPE, Popen
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MuscleCommandline
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# IgCore imports
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from IgCore import default_out_args, default_separator
from IgCore import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from IgCore import printLog
from IgCore import manageProcesses
from DbCore import feedDbQueue, processDbQueue, collectDbQueue
from DbCore import DbData, DbResult

# Globals
# TODO:  not convinced this is the best way to deal with fails in grouping. a class may be better.
FAIL_GROUP_KEY = 'FAIL_GROUP_KEY'

# Defaults
default_muscle_exec = r'/usr/local/bin/muscle'

# TODO:  maybe not bothering with 'set' is best. can just work off field identity
def groupRecords(records, fields=None, calls=['v', 'j'], mode='gene', action='first',
                 separator=default_separator):
    """
    Groups IgRecords based on gene or annotation

    Arguments:
    records = an iterator of IgRecord objects to group
    fields = gene field to group by
    calls = allele calls to use for grouping;
            one or more of ('v', 'd', 'j')
    mode = specificity of alignment call to use for allele call fields;
           one of ('allele', 'gene')
    action = how to handle multiple calls in alleles fields;
             one of ('first', 'set')
    separator = the delimiter separating values within a field

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
            ann = [rec.toDict().get(k, None) for k in fields]
            return tuple(chain(vdj, ann))
    elif mode == 'gene' and fields is not None:
        def _get_key(rec, call, action):
            vdj = rec.getGeneCalls(calls, action)
            ann = [rec.toDict().get(k, None) for k in fields]
            return tuple(chain(vdj, ann))

    rec_index = {}
    for rec in records:
        key = _get_key(rec, calls, action)
        # Assigned grouped records to individual keys and all failed to a single key
        if all([k is not None for k in key]):
            rec_index.setdefault(key, []).append(rec)
        else:
            rec_index.setdefault(FAIL_GROUP_KEY, []).append(rec)

    return rec_index


def alignRecords(data, seq_fields, muscle_exec=default_muscle_exec):
    """
    Multiple aligns sequence fields

    Arguments:
    data = a DbData object with IgRecords to process
    seq_fields = the sequence fields to multiple align
    muscle_exec = the MUSCLE executable

    Returns:
    a list of modified IgRecords with multiple aligned sequence fields
    """
    # Define return object
    result = DbResult(data.id, data.data)
    result.results = data.data
    result.valid = True

    for f in seq_fields:
        seq_list = [SeqRecord(r.getField(f), id=r.id) for r in data.data]
        seq_aln = alignSeqSet(seq_list, muscle_exec=muscle_exec)
        if seq_aln is not None:
            for i, r in enumerate(result.results):
                r.annotations['%s_ALIGN' % f] = str(seq_aln[i].seq)
        else:
            result.valid = False

    #for r in result.results:  print r.annotations
    return result


def alignSeqSet(seq_list, muscle_exec=default_muscle_exec):
    """
    Multiple aligns a set of sequences

    Arguments: 
    seq_list = a list of SeqRecord objects to align
    muscle_exec = the MUSCLE executable
    
    Returns: 
    a MultipleSeqAlignment object containing the alignment
    """
    # Return sequence if only one sequence in seq_list
    if len(seq_list) < 2:
        align = MultipleSeqAlignment(seq_list)
        return align
    
    # Set MUSCLE command
    cmd = MuscleCommandline(muscle_exec, maxiters=2, diags=True)

    # Convert sequences to FASTA and write to string
    stdin_handle = StringIO()
    SeqIO.write(seq_list, stdin_handle, 'fasta')
    stdin_str = stdin_handle.getvalue()
    stdin_handle.close()
    
    # Open MUSCLE process
    child = Popen(str(cmd), stdin=PIPE, stdout=PIPE, stderr=PIPE,
                  shell=(sys.platform != 'win32'))

    # Send sequences to MUSCLE stdin and retrieve stdout, stderr
    stdout_str, __ = child.communicate(stdin_str)

    # Capture sequences from MUSCLE stdout
    stdout_handle = StringIO(stdout_str)
    align = AlignIO.read(stdout_handle, 'fasta')
    stdout_handle.close()

    return align


def gapSeq(db_file, seq_fields, group_func, align_func, group_args={}, align_args={},
           out_args=default_out_args, nproc=None, queue_size=None):
    """
    Performs a multiple alignment on sets of sequences

    Arguments: 
    db_file = filename of the input database
    seq_fields = the sequence fields to multiple align
    group_func = function to use to group records
    align_func = function to use to multiple align sequence groups
    group_args = dictionary of arguments to pass to group_func
    align_args = dictionary of arguments to pass to gap_func
    out_args = common output argument dictionary from parseCommonArgs
    nproc = the number of processQueue processes;
            if None defaults to the number of CPUs
    queue_size = maximum size of the argument queue;
                 if None defaults to 2*nproc
                      
    Returns: 
    a tuple of (align-pass, align-fail) filenames
    """
    # Define subcommand label dictionary
    cmd_dict = {alignRecords:'align'}
    
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'GapSeq'
    log['COMMAND'] = cmd_dict.get(align_func, align_func.__name__)
    log['FILE'] = os.path.basename(db_file)
    log['SEQ_FIELDS'] = ','.join(seq_fields)
    if 'group_fields' in group_args: log['GROUP_FIELDS'] = ','.join(group_args['group_fields'])
    if 'mode' in group_args: log['MODE'] = group_args['mode']
    if 'action' in group_args: log['ACTION'] = group_args['action']
    log['NPROC'] = nproc
    printLog(log)
 
    # Define feeder function and arguments
    feed_func = feedDbQueue
    feed_args = {'db_file': db_file,
                 'group_func': group_func,
                 'group_args': group_args}
    # Define worker function and arguments
    align_args['seq_fields'] = seq_fields
    work_func = processDbQueue
    work_args = {'process_func': align_func,
                 'process_args': align_args}
    # Define collector function and arguments
    collect_func = collectDbQueue
    collect_args = {'db_file': db_file,
                    'task_label': 'gap',
                    'out_args': out_args,
                    'add_fields': ['%s_ALIGN' % f for f in seq_fields]}
    
    # Call process manager
    result = manageProcesses(feed_func, work_func, collect_func, 
                             feed_args, work_args, collect_args, 
                             nproc, queue_size)
        
    # Print log
    result['log']['END'] = 'GapSeq'
    printLog(result['log'])
        
    return result['out_files']


def getArgParser():
    """
    Defines the ArgumentParser

    Arguments: 
    None
                      
    Returns: 
    an ArgumentParser object
    """
    # Define output file names and header fields
    fields = textwrap.dedent(
             '''
             output files:
               gap-pass       database with multiple aligned sequences.
               gap-fail       database with records failing alignment.

             required fields:
               SEQUENCE_ID
               SEQUENCE_VDJ
               V_CALL
               J_CALL

             output fields:
               <sequence field>_ALIGN
             ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            version='%(prog)s:' + ' v%s-%s' %(__version__, __date__),
                            formatter_class=CommonHelpFormatter)
    subparsers = parser.add_subparsers(title='subcommands', dest='command', metavar='',
                                       help='Gapping method')
    
    # Parent parser
    parser_parent = getCommonArgParser(seq_in=False, seq_out=False, db_in=True,
                                       multiproc=True)
    # MUSCLE mode argument parser
    parser_align = subparsers.add_parser('align', parents=[parser_parent],
                                         formatter_class=CommonHelpFormatter,
                                         help='Multiple aligns sequence groups using MUSCLE')

    parser_align.add_argument('--sf', nargs='+', action='store', dest='seq_fields',
                              default=['SEQUENCE_VDJ', 'GERMLINE_VDJ'],
                              help='The sequence field to multiple align within each group.')
    parser_align.add_argument('--gf', nargs='+', action='store', dest='group_fields',
                              default=None,
                              help='Additional (not allele call) fields to use for grouping.')
    parser_align.add_argument('-a', nargs='+', action='store', dest='calls',
                              choices=('v', 'd', 'j'), default=['v', 'j'],
                              help='Segment calls (allele assignments) to use for grouping.')
    parser_align.add_argument('--mode', action='store', dest='mode',
                              choices=('allele', 'gene'), default='gene',
                              help='''Specifies whether to use the V(D)J allele or gene when
                                   an allele call field (--cf) is specified.''')
    parser_align.add_argument('--act', action='store', dest='action', default='first',
                              choices=('first'),
                              help='''Specifies how to handle multiple values within default
                                   allele call fields. Currently, only "first" is supported.''')
    parser_align.add_argument('--exec', action='store', dest='muscle_exec',
                              default=default_muscle_exec,
                              help='The location of the MUSCLE executable')
    parser_align.set_defaults(group_func=groupRecords, align_func=alignRecords)
    
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
    if 'seq_fields' in args_dict and args_dict['seq_fields'] is not None:
        args_dict['seq_fields'] = [f.upper() for f in args_dict['seq_fields']]
    if 'group_fields' in args_dict and args_dict['group_fields'] is not None:
        args_dict['group_fields'] = [f.upper() for f in args_dict['group_fields']]

    # Check if a valid MUSCLE executable was specified for muscle mode
    if args.command == 'align' and not os.path.isfile(args.muscle_exec):
        parser.error('%s does not exist' % args.muscle_exec)

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

    # Define align_args
    if args_dict['align_func'] is alignRecords:
        args_dict['align_args'] = {'muscle_exec':args_dict['muscle_exec']}
        del args_dict['muscle_exec']


    # Call gapSeq for each input file
    del args_dict['command']
    del args_dict['db_files']
    for f in args.__dict__['db_files']:
        args_dict['db_file'] = f
        gapSeq(**args_dict)

