#!/usr/bin/env python
"""
Multiple aligns sequence fields
"""

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2013 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.0'
__date__      = '2015.03.30'

# Imports
import csv, os, sys, textwrap
from argparse import ArgumentParser
from collections import deque, OrderedDict
from cStringIO import StringIO
from itertools import izip
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
from IgCore import printLog, printMessage, printProgress, getFileType, getOutputHandle
from IgCore import manageProcesses
from DbCore import DbData, DbResult

# TODO: remove when feed, process and collect functions are ported
from IgCore import countSeqFile, readSeqFile, countSeqSets, SeqData, SeqResult

# Defaults
default_muscle_exec = r'/usr/local/bin/muscle'


def feedDbQueue(alive, data_queue, seq_file, index_func=None, index_args={}):
    """
    Feeds the data queue with SeqRecord objects

    Arguments:
    alive = a multiprocessing.Value boolean controlling whether processing
            continues; when False function returns
    data_queue = a multiprocessing.Queue to hold data for processing
    seq_file = the sequence file to read input from
    index_func = the function to use to define sequence sets
                 if None do not index sets and feed individual records
    index_args = a dictionary of arguments to pass to index_func

    Returns:
    None
    """
    try:
        # Read input file and index sequence sets if required
        if index_func is None:
            seq_iter = readSeqFile(seq_file)
            data_iter = ((s.id, s) for s in seq_iter)
        else:
            seq_dict = readSeqFile(seq_file, index=True)
            index_dict = index_func(seq_dict, **index_args)
            data_iter = ((k, [seq_dict[i] for i in v]) \
                         for k, v in index_dict.iteritems())
    except:
        alive.value = False
        raise

    try:
        # Iterate over data_iter and feed data queue
        while alive.value:
            # Get data from queue
            if data_queue.full():  continue
            else:  data = next(data_iter, None)
            # Exit upon reaching end of iterator
            if data is None:  break

            # Feed queue
            data_queue.put(SeqData(*data))
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
    except:
        alive.value = False
        raise

    return None

def processDbQueue(alive, data_queue, result_queue, process_func, process_args={}):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments:
    alive = a multiprocessing.Value boolean controlling whether processing
            continues; when False function returns
    data_queue = a multiprocessing.Queue holding data to process
    result_queue = a multiprocessing.Queue to hold processed results
    process_func = the function to use for filtering sequences
    process_args = a dictionary of arguments to pass to process_func

    Returns:
    None
    """
    try:
        # Iterator over data queue until sentinel object reached
        while alive.value:
            # Get data from queue
            if data_queue.empty():  continue
            else:  data = data_queue.get()
            # Exit upon reaching sentinel
            if data is None:  break

            # Perform work
            result = process_func(data, **process_args)

            #import cProfile
            #prof = cProfile.Profile()
            #result = prof.runcall(process_func, data, **process_args)
            #prof.dump_stats('worker-%d.prof' % os.getpid())

            # Feed results to result queue
            result_queue.put(result)
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
    except:
        alive.value = False
        sys.stderr.write('Error processing sequence with ID: %s.\n' % data.id)
        raise

    return None


def collectDbQueue(alive, result_queue, collect_queue, seq_file,
                    task_label, out_args, index_field=None):
    """
    Pulls from results queue, assembles results and manages log and file IO

    Arguments:
    alive = a multiprocessing.Value boolean controlling whether processing
            continues; when False function returns
    result_queue = a multiprocessing.Queue holding worker results
    collect_queue = a multiprocessing.Queue to store collector return values
    seq_file = the sample sequence file name
    task_label = the task label used to tag the output files
    out_args = common output argument dictionary from parseCommonArgs
    index_field = the field defining set membership for sequence sets
                  if None data queue contained individual records

    Returns:
    None
    (adds a dictionary of {log: log object, out_files: output file names} to collect_queue)
    """
    try:
        # Count records
        if index_field is None:
            result_count = countSeqFile(seq_file)
        else:
            result_count = countSeqSets(seq_file, index_field, out_args['delimiter'])

        # Define output format
        out_type = getFileType(seq_file) if out_args['out_type'] is None \
                   else out_args['out_type']

        # Defined valid alignment output handle
        pass_handle = getOutputHandle(seq_file,
                                      '%s-pass' % task_label,
                                      out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'],
                                      out_type=out_type)
        # Defined failed alignment output handle
        if out_args['failed']:
            fail_handle = getOutputHandle(seq_file,
                                          '%s-fail'  % task_label,
                                          out_dir=out_args['out_dir'],
                                          out_name=out_args['out_name'],
                                          out_type=out_type)
        else:
            fail_handle = None

        # Define log handle
        if out_args['log_file'] is None:
            log_handle = None
        else:
            log_handle = open(out_args['log_file'], 'w')
    except:
        alive.value = False
        raise

    try:
        # Iterator over results queue until sentinel object reached
        start_time = time()
        set_count = seq_count = pass_count = fail_count = 0
        while alive.value:
            # Get result from queue
            if result_queue.empty():  continue
            else:  result = result_queue.get()
            # Exit upon reaching sentinel
            if result is None:  break

            # Print progress for previous iteration
            printProgress(set_count, result_count, 0.05, start_time)

            # Update counts for current iteration
            set_count += 1
            seq_count += result.data_count

            # Write log
            printLog(result.log, handle=log_handle)

            # Write alignments
            if result:
                pass_count += 1
                SeqIO.write(result.results, pass_handle, out_type)
            else:
                fail_count += 1
                if fail_handle is not None:
                    SeqIO.write(result.data, fail_handle, out_type)
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None

        # Print total counts
        printProgress(set_count, result_count, 0.05, start_time)

        # Update return values
        log = OrderedDict()
        log['OUTPUT'] = os.path.basename(pass_handle.name)
        log['SEQUENCES'] = seq_count
        if index_field is not None:
            log['SETS'] = set_count
        log['PASS'] = pass_count
        log['FAIL'] = fail_count
        collect_dict = {'log':log, 'out_files': [pass_handle.name]}
        collect_queue.put(collect_dict)

        # Close file handles
        pass_handle.close()
        if fail_handle is not None:  fail_handle.close()
        if log_handle is not None:  log_handle.close()
    except:
        alive.value = False
        raise

    return None


def groupRecords(records, gene_fields=None, ann_fields=None, action='first',
                 gene_mode='gene', separator=default_separator):
    """
    Groups IgRecords bases on gene or annotation

    Arguments:
    records = an iterator of IgRecord objects to group
    gene_fields = gene fields to group by
    ann_fields = annotation fields to group by
    gene_mode = specificity of alignment call to use;
                one of ('allele', 'gene')
    action = how to handle multiple value fields;
             one of ('first', 'set')
    separator = the delimiter separating values within a field

    Returns:
    dictionary of grouped records
    """
    pass


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


def gapSeq(db_file, group_func, gap_func, group_args={}, gap_args={},
           out_args=default_out_args, nproc=None, queue_size=None):
    """
    Performs a multiple alignment on sets of sequences

    Arguments: 
    db_file = filename of the input database
    group_func = function to use to group records
    gap_func = function to use to multiple align sequence groups
    group_args = dictionary of arguments to pass to group_func
    gap_args = dictionary of arguments to pass to gap_func
    out_args = common output argument dictionary from parseCommonArgs
    nproc = the number of processQueue processes;
            if None defaults to the number of CPUs
    queue_size = maximum size of the argument queue;
                 if None defaults to 2*nproc
                      
    Returns: 
    a tuple of (gap-pass, gap-fail) filenames
    """
    # Define subcommand label dictionary
    cmd_dict = {alignSeqSet:'align'}
    
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'AlignSets'
    log['COMMAND'] = cmd_dict.get(gap_func, gap_func.__name__)
    log['FILE'] = os.path.basename(seq_file)
    if 'mode' in gap_args: log['MODE'] = gap_args['mode']
    if 'field' in gap_args: log['OFFSET_FIELD'] = gap_args['field']
    log['NPROC'] = nproc
    printLog(log)
 
    # Define feeder function and arguments
    feed_func = feedDbQueue
    feed_args = {'db_file': db_file,
                 'group_func': groupRecords,
                 'group_args': group_args}
    # Define worker function and arguments
    work_func = processDbQueue
    work_args = {'gap_func': gap_func,
                 'gap_args': gap_args}
    # Define collector function and arguments
    collect_func = collectDbQueue
    collect_args = {'db_file': db_file,
                    'task_label': 'gap',
                    'out_args': out_args}
    
    # Call process manager
    result = manageProcesses(feed_func, work_func, collect_func, 
                             feed_args, work_args, collect_args, 
                             nproc, queue_size)
        
    # Print log
    result['log']['END'] = 'AlignSets'
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
             gap-pass     database with gapped sequences.
             gap-fail     database with records failing gapping.

         required fields:
             SEQUENCE_ID
             SEQUENCE
             V_START
             V_LENGTH
             J_START
             J_LENGTH

         output fields:
             SEQUENCE_ALN
             SEQUENCE_IMGT
         ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            version='%(prog)s:' + ' v%s-%s' %(__version__, __date__),
                            formatter_class=CommonHelpFormatter)
    subparsers = parser.add_subparsers(title='subcommands', dest='command', metavar='',
                                       help='Gapping method')
    
    # Parent parser    
    parser_parent = getCommonArgParser(multiproc=True)

    # MUSCLE mode argument parser
    parser_align = subparsers.add_parser('align', parents=[parser_parent],
                                         formatter_class=CommonHelpFormatter,
                                         help='Multiple align groups of sequence using MUSCLE')
    parser_align.add_argument('-f', nargs='+', action='store', dest='fields', default=None,
                             help='Additional fields to use for grouping')
    parser_align.add_argument('--mode', action='store', dest='mode',
                             choices=('allele', 'gene'), default='gene',
                             help='Specifies whether to use the V(D)J allele or gene for initial grouping')
    parser_align.add_argument('--act', action='store', dest='action', default='set',
                             choices=('first', 'set'),
                             help='Specifies how to handle multiple V(D)J assignments for initial grouping')

    parser_align.add_argument('--exec', action='store', dest='muscle_exec', default=default_muscle_exec,
                               help='The location of the MUSCLE executable')
    parser_align.set_defaults(gap_func=alignSeqSet)
    
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
    if 'barcode_field' in args_dict and args_dict['barcode_field']:
        args_dict['barcode_field'] = args_dict['barcode_field'].upper()
    if 'primer_field' in args_dict and args_dict['primer_field']:
        args_dict['primer_field'] = args_dict['primer_field'].upper()
    
    # Check if a valid MUSCLE executable was specified for muscle mode
    if args.command == 'align' and not os.path.isfile(args.muscle_exec):
        parser.error('%s does not exist' % args.muscle_exec)
    
    # Define gap_args
    if args_dict['gap_func'] is alignSeqSet:
        args_dict['gap_args'] = {'muscle_exec':args_dict['muscle_exec']}
        del args_dict['muscle_exec']
        
    # Call gapSeq for each input file
    del args_dict['command']
    del args_dict['db_files']
    for f in args.__dict__['db_files']:
        args_dict['db_file'] = f
        gapSeq(**args_dict)

