#!/usr/bin/env python
"""
Assign Ig sequences into clones
"""

__author__    = 'Namita Gupta, Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2013 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.0'
__date__      = '2013.10.10'

# Imports
import csv, os, re, sys
import multiprocessing as mp
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import OrderedDict
from time import time

# IgCore imports
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from DbCore import default_delimiter, default_out_args
from DbCore import getCommonParser, parseCommonArgs
from DbCore import getOutputHandle, printLog, printProgress

# Defaults

# >>> into DbCore when finished
class IgRecord:
    """
    A class defining a V(D)J germline sequence alignment
    """
    def __init__(self, i, s, v, d, j):
        self.id = i
        self.seq = s
        self.v_call = v
        self.d_call = d
        self.j_call = j
    
    def getVAllele(self):
        return self.v_call
    
    def getDAllele(self):
        return self.d_call
    
    def getJAllele(self):
        return self.j_call


# >>> into DbCore when finished
def readDbFile(db_file):
    """
    Reads database files

    Arguments: 
    db_file = a tab delimited database file
    
    Returns: 
    a tuple of (input file type, database record iterator)
    """
    # Read and check file
    try:
        db_type = os.path.splitext(db_file)[1].lower().lstrip('.')
        db_records = ['test_record']
    except IOError:
        sys.exit('ERROR:  File %s cannot be read' % db_file)
    except:
        sys.exit('ERROR:  File %s is invalid' % db_file)
    
    return db_type, db_records


# >>> into DbCore when finished
def countDbFile(db_file):
    """
    Counts the records in database files

    Arguments: 
    db_file = a tab delimited database file

    Returns: 
    the count of records in the database file
    """
    # Count records and check file
    try:
        db_type = os.path.splitext(db_file)[1].lower().lstrip('.')
        db_count = 1
    except IOError:
        sys.exit('ERROR:  File %s cannot be read' % db_file)
    except:
        sys.exit('ERROR:  File %s is invalid' % db_file)
    else:
        if db_count == 0:  sys.exit('ERROR:  File %s is empty' % db_file)
        
    return db_count


# >>> into DbCore when finished
def feedQueue(data_queue, nproc, db_iter):
    """
    Feeds the data queue with Ig records for processQueue processes

    Arguments: 
    data_queue = a multiprocessing.Queue to hold data for processing
    nproc = the number of processQueue processes
    db_iter = an iterator of IgRecords returned by readDbFile
    
    Returns: 
    None
    """
    # Iterate over db_ter and define processQueue arguments 
    count = 0
    for rec in db_iter:
        # Feed queue
        #data_queue.put({'id':rec.id, 'rec':rec})
        data_queue.put({'id':count, 'rec':rec})
        count += 1
    
    # Add sentinel object for each processQueue process
    for __ in range(nproc):
        data_queue.put(None)

    return None


def distanceClones(ig_set):
    """
    Separates a set of IgRecords into clones

    Arguments: 
    ig_set = a iterator of IgRecords
    
    Returns: 
    a list of IgRecords with a clone annotation
    """

    return None


def processQueue(data_queue, result_queue, clone_func, clone_args, delimiter=default_delimiter):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments: 
    data_queue = a multiprocessing.Queue holding data to process
    result_queue = a multiprocessing.Queue to hold processed results
    clone_func = the function to call for clonal assignment
    clone_args = a dictionary of arguments to pass to clone_func
    delimiter = a tuple of delimiters for (fields, values, value lists) 

    Returns: 
    None
    """
    # Iterator over data queue until sentinel object reached
    for args in iter(data_queue.get, None):
        in_rec = args['rec']
        
        # Define result dictionary for iteration
        results = {'id':args['id'],
                   'in_rec':in_rec,
                   'out_rec':None,
                   'pass':False,
                   'log':OrderedDict([('ID', args['id'])])}
 
        # Assign clones
        preclones = clone_func(in_rec)
        #clones = clone_func(in_rec, **clone_args)

        # Feed results to result queue
        result_queue.put(results)

    return None

# >>> into DbCore when finished
def collectQueue(result_queue, result_count, collect_dict, task_label, db_file, out_args):
    """
    Assembles results from a queue of individual sequence results and manages log/file I/O

    Arguments: 
    result_queue = a multiprocessing.Queue holding processQueue results
    result_count = the total number of results expected
    collect_dict = a multiprocessing.Manager.dict to store return values
    task_label = the task label used to tag the output files
    db_file = the input database file name
    out_args = common output argument dictionary from parseCommonArgs
    
    Returns: 
    None
    (adds 'log' and 'out_files' to collect_dict)
    """
    # Defined successful output handle
    pass_handle = getOutputHandle(db_file, 
                                  '%s-pass' % task_label, 
                                  out_dir=out_args['out_dir'], 
                                  out_name=out_args['out_name'], 
                                  out_type=out_args['out_type'])
    # Defined failed alignment output handle
    if out_args['clean']:   
        fail_handle = None
    else:  
        fail_handle = getOutputHandle(db_file, 
                                      '%s-fail' % task_label, 
                                      out_dir=out_args['out_dir'], 
                                      out_name=out_args['out_name'], 
                                      out_type=out_args['out_type'])
    # Define log handle
    if out_args['log_file'] is None:  
        log_handle = None
    else:  
        log_handle = open(out_args['log_file'], 'w')
        
    # Iterator over results queue until sentinel object reached
    start_time = time()
    rec_count = pass_count = fail_count = 0
    for result in iter(result_queue.get, None): 
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time) 
        
        # Update counts
        rec_count += 1

        # Write log
        printLog(result['log'], handle=log_handle)
 
    # Print total counts
    printProgress(rec_count, result_count, 0.05, start_time) 

    # Update return list
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    collect_dict['log'] = log
    collect_dict['out_files'] = [pass_handle.name]

    # Close file handles
    pass_handle.close()
    if fail_handle is not None:  fail_handle.close()
    if log_handle is not None:  log_handle.close()
        
    return None


def defineClones(db_file, mode, action, fields=None, out_args=default_out_args, 
                 nproc=None, queue_size=None):
    """
    Define clonally related sequences
    
    Arguments:
    db_file = filename of input database
    mode = specificity of alignment call to use for assigning preclones;
           one of ['allele', 'gene']
    action = how to handle multiple value fields when assigning preclones;
             one of ['first', 'intersect', 'union']
    fields = additional fields to use to group preclones
             if None use only V_CALL, J_CALL, J_LENGTH
    out_args = common output argument dictionary from parseCommonArgs
    nproc = the number of processQueue processes;
            if None defaults to the number of CPUs
    queue_size = maximum size of the argument queue;
                 if None defaults to 2*nproc    
    
    Returns:
    a list of successful output file names
    """
    # Define number of processes and queue size
    if nproc is None:  nproc = mp.cpu_count()
    if queue_size is None:  queue_size = nproc * 2
    
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'DefineClones'
    log['DB_FILE'] = os.path.basename(db_file)
    log['MODE'] = mode
    log['ACTION'] = action
    log['FIELDS'] = ','.join([str(x) for x in fields]) \
                    if fields is not None else None
    log['NPROC'] = nproc
    printLog(log)
    
    # Define database iterator
    in_type, db_iter = readDbFile(db_file)
    if out_args['out_type'] is None:  out_args['out_type'] = in_type
    
    # Count records
    result_count = countDbFile(db_file)
    
    clone_func = distanceClones
    clone_args = {'mode':mode, 
                  'action':action,
                  'fields':fields}
    
    # Define shared data objects
    manager = mp.Manager()
    collect_dict = manager.dict()
    data_queue = mp.Queue(queue_size)
    result_queue = mp.Queue(queue_size)
    
    # Initiate feeder process
    feeder = mp.Process(target=feedQueue, args=(data_queue, nproc, db_iter))
    feeder.start()

    # Initiate processQueue processes
    workers = []
    for __ in range(nproc):
        w = mp.Process(target=processQueue, args=(data_queue, result_queue, clone_func, clone_args, 
                                                  out_args['delimiter']))
        w.start()
        workers.append(w)

    # Initiate collector process
    collector = mp.Process(target=collectQueue, args=(result_queue, result_count, collect_dict, 
                                                      'clone', db_file, out_args))
    collector.start()

    # Wait for feeder and worker processes to finish, add sentinel to result_queue
    feeder.join()
    for w in workers:  w.join()
    result_queue.put(None)
    
    # Wait for collector process to finish and shutdown manager
    collector.join()
    log = collect_dict['log']
    out_files = collect_dict['out_files']
    manager.shutdown()
    
    # Print log
    log['END'] = 'DefineClones'
    printLog(log)
        
    return out_files


def getArgParser():
    """
    Defines the ArgumentParser

    Arguments: 
    None
                      
    Returns: 
    an ArgumentParser object
    """
    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, 
                            version='%(prog)s:' + ' v%s-%s' %(__version__, __date__), 
                            parents=[getCommonParser(seq_in=False, seq_out=False, db_in=True, multiproc=True)], 
                            formatter_class=ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-f', nargs='+', action='store', dest='fields', default=None,
                        help='Additional fields to use for grouping clones (non VDJ)')
    parser.add_argument('--mode', action='store', dest='mode', 
                        choices=('allele', 'gene'), default='gene', 
                        help='Specifies whether to use the V(D)J allele or gene for preclone assignment')
    parser.add_argument('--act', action='store', dest='action', default='first',
                        choices=['first', 'intersect', 'union'],
                        help='Specifies how to handle multiple V(D)J assignments for preclone assignment') 
        
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
    if args_dict['fields']:  args_dict['fields'] = map(str.upper, args_dict['fields'])

    # Call defineClones
    del args_dict['db_files']
    for f in args.__dict__['db_files']:
        args_dict['db_file'] = f
        defineClones(**args_dict)