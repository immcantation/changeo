#!/usr/bin/env python
"""
Assign Ig sequences into clones
"""

__author__    = 'Jason Anthony Vander Heiden, Namita Gupta, Gur Yaari, Mohamed Uduman'
__copyright__ = 'Copyright 2013 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.0'
__date__      = '2013.10.14'

# Imports
import csv, linecache, os, re, sys
import multiprocessing as mp
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import OrderedDict
from time import time

# IgCore imports
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from DbCore import default_separator, default_out_args
from DbCore import getCommonParser, parseCommonArgs
from DbCore import getOutputHandle, printLog, printProgress
from DbCore import countDbRecords, readDbFile, getFileType

# Defaults


def indexPreclones(db_iter, fields=None, action='first', separator=default_separator):
    """
    Identifies preclonal groups by V, J and junction length

    Arguments: 
    db_iter = an iterator of IgRecords defined by readDbFile
    fields = additional fields to use to group preclones;
             if None use only V_CALL, J_CALL, JUNCTION_GAP_LENGTH
    action = how to handle multiple value fields when assigning preclones;
             one of ('first', 'set')
    separator = the delimiter separating values within a field 
    
    Returns: 
    an iterator of 
    """
    clone_index = {}
    for i, rec in enumerate(db_iter):
        key = (rec.getVAllele('first'), rec.getJAllele('first'), len(rec.junction))
        #key = (rec.getVAllele('set'), rec.getJAllele('set'), len(rec.junction))
        #key = (rec.getVGene('set'), rec.getJGene('set'), len(rec.junction))
        print key
        #print rec.__dict__
        if all([k is not None for k in key]):
            clone_index.setdefault(key, []).append(rec)
            #clone_index.setdefault(key, []).append(i)
        
    return clone_index


def distanceClones(rec_list):
    """
    Separates a set of IgRecords into clones

    Arguments: 
    rec_list = an iterator of IgRecords
    
    Returns: 
    a list of IgRecords with a clone annotation
    """
    #print ''
    for rec in rec_list:
        pass
        #print rec.id, rec.junction
        #print rec.__dict__
        
    return None


# >>> into DbCore when finished?
def feedQueue(data_queue, nproc, db_file):
    """
    Feeds the data queue with Ig records

    Arguments: 
    data_queue = a multiprocessing.Queue to hold data for processing
    nproc = the number of processQueue processes
    db_file = the Ig record database file
    
    Returns: 
    None
    """
    # Iterate over Ig records and feed data queue
    db_iter = readDbFile(db_file)
    clone_dict = indexPreclones(db_iter)
    
    #db_iter = readDbFile(db_file)
    for k, v in clone_dict.iteritems():
        # Feed queue
        data_queue.put({'id':k, 'rec_list':v})
    
    # Add sentinel object for each processQueue process
    for __ in range(nproc):
        data_queue.put(None)

    return None


def processQueue(data_queue, result_queue, clone_func, clone_args, separator=default_separator):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments: 
    data_queue = a multiprocessing.Queue holding data to process
    result_queue = a multiprocessing.Queue to hold processed results
    clone_func = the function to call for clonal assignment
    clone_args = a dictionary of arguments to pass to clone_func
    separator = the delimiter separating values within a field 

    Returns: 
    None
    """
    # Iterator over data queue until sentinel object reached
    for args in iter(data_queue.get, None):
        rec_list = args['rec_list']
        
        # Define result dictionary for iteration
        results = {'id':args['id'],
                   'rec_list':rec_list,
                   'out_rec':None,
                   'pass':False,
                   'log':OrderedDict([('ID', args['id'])])}
 
        # Assign clones
        clones = clone_func(rec_list, **clone_args)

        # Feed results to result queue
        result_queue.put(results)

    return None

# >>> into DbCore when finished?
def collectQueue(result_queue, collect_dict, db_file, out_args):
    """
    Assembles results from a queue of individual sequence results and manages log/file I/O

    Arguments: 
    result_queue = a multiprocessing.Queue holding processQueue results
    result_count = the total number of results expected
    collect_dict = a multiprocessing.Manager.dict to store return values
    db_file = the input database file name
    out_args = common output argument dictionary from parseCommonArgs
    
    Returns: 
    None
    (adds 'log' and 'out_files' to collect_dict)
    """
    # Count records and define output format 
    out_type = getFileType(db_file) if out_args['out_type'] is None \
               else out_args['out_type']
    result_count = countDbRecords(db_file)

    # Defined successful output handle
    pass_handle = getOutputHandle(db_file, 
                                  'clone-pass', 
                                  out_dir=out_args['out_dir'], 
                                  out_name=out_args['out_name'], 
                                  out_type=out_type)
    # Defined failed alignment output handle
    if out_args['clean']:   
        fail_handle = None
    else:  
        fail_handle = getOutputHandle(db_file, 
                                      'clone-fail', 
                                      out_dir=out_args['out_dir'], 
                                      out_name=out_args['out_name'], 
                                      out_type=out_type)
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
        
        if result['pass']:
            pass_count += 1
        else:
            fail_count += 1

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
           one of ('allele', 'gene')
    action = how to handle multiple value fields when assigning preclones;
             one of ('first', 'set')
    fields = additional fields to use to group preclones;
             if None use only V_CALL, J_CALL, JUNCTION_GAP_LENGTH
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
    
    clone_func = distanceClones
    #clone_args = {'mode':mode, 
    #              'action':action,
    #              'fields':fields}
    clone_args = {}
    
    # Define shared data objects
    manager = mp.Manager()
    collect_dict = manager.dict()
    data_queue = mp.Queue(queue_size)
    result_queue = mp.Queue(queue_size)
    
    # Initiate feeder process
    feeder = mp.Process(target=feedQueue, args=(data_queue, nproc, db_file))
    feeder.start()

    # Initiate processQueue processes
    workers = []
    for __ in range(nproc):
        w = mp.Process(target=processQueue, args=(data_queue, result_queue, clone_func, clone_args, 
                                                  out_args['delimiter']))
        w.start()
        workers.append(w)

    # Initiate collector process
    collector = mp.Process(target=collectQueue, args=(result_queue, collect_dict, db_file, 
                                                      out_args))
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
                        choices=('first', 'set'),
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