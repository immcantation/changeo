#!/usr/bin/env python
"""
Assign Ig sequences into clones
"""

__author__    = 'Jason Anthony Vander Heiden, Namita Gupta, Gur Yaari, Mohamed Uduman'
__copyright__ = 'Copyright 2013 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.0'
__date__      = '2013.11.10'

# Imports
import csv, linecache, os, re, sys
import multiprocessing as mp
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import OrderedDict
from itertools import chain
from time import time

# IgCore imports
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from IgCore import default_separator, default_out_args
from IgCore import getCommonArgParser, parseCommonArgs
from IgCore import getFileType, getOutputHandle, printLog, printProgress
from DbCore import countDbFile, readDbFile, getDbWriter

# R imports
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage as STAP

# >>> THIS IS AN AWFUL HACK. NEEDS FIXING.
# Source ClonesByDist
py_wd = os.getcwd()
r_lib = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'rlib')
os.chdir(r_lib)
with open('ClonesByDist.R', 'r') as f:
    r_script = ''.join(f.readlines())
ClonesByDist = STAP(r_script, 'ClonesByDist')
os.chdir(py_wd)

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
    for rec in db_iter:
        key = (rec.getVAllele('first'), rec.getJAllele('first'), len(rec.junction))
        #key = (rec.getVAllele('set'), rec.getJAllele('set'), len(rec.junction))
        #key = (rec.getVGene('set'), rec.getJGene('set'), len(rec.junction))
        #print key
        #print rec.__dict__
        #print rec.to_dict()
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
    a list of lists defining [IgRecords clonal groups]
    """
    # Define unique junction mapping
    junc_map = {}
    for r in rec_list:
        junc_map.setdefault(str(r.junction.upper()), []).append(r)
    
    # Process records    
    if len(junc_map) > 1:
        # Call R function
        junc_str = '|'.join(junc_map.keys())
        #junc_str = '|'.join([str(r.junction.upper()) for r in rec_list])
        clone_str = ClonesByDist.getClones(junc_str, 5)
        # Parse R output
        clone_list = []
        for clone in str(clone_str[0]).split('|'):
            junc_r = clone.split(',')
            #print list(chain(*[junc_map[j.upper()] for j in junc_r]))
            clone_list.append(list(chain(*[junc_map[j.upper()] for j in junc_r])))
    else:
        clone_list = [rec_list]

    return clone_list


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
                   'clone_list':None,
                   'pass':False,
                   'log':OrderedDict([('ID', args['id'])])}
 
        # Assign clones
        clone_list = clone_func(rec_list, **clone_args)
        if clone_list is not None:
            results['clone_list'] = clone_list
            results['pass'] = True
            
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
    result_count = countDbFile(db_file)

    # Defined successful output handle
    pass_handle = getOutputHandle(db_file, 
                                  out_label='clone-pass', 
                                  out_dir=out_args['out_dir'], 
                                  out_name=out_args['out_name'], 
                                  out_type=out_type)
    pass_writer = getDbWriter(pass_handle, db_file, add_fields='CLONE_ID')
    
    # Defined failed alignment output handle
    if out_args['clean']:   
        fail_handle = None
        fail_writer = None
    else:  
        fail_handle = getOutputHandle(db_file, 
                                      out_label='clone-fail', 
                                      out_dir=out_args['out_dir'], 
                                      out_name=out_args['out_name'], 
                                      out_type=out_type)
        fail_writer = getDbWriter(fail_handle, db_file)

    # Define log handle
    if out_args['log_file'] is None:  
        log_handle = None
    else:  
        log_handle = open(out_args['log_file'], 'w')
        
    # Iterator over results queue until sentinel object reached
    start_time = time()
    clone_count = pass_count = fail_count = 0
    for result in iter(result_queue.get, None): 
        # Print progress for previous iteration
        printProgress(pass_count + fail_count, result_count, 0.05, start_time) 
        
        if result['pass']:
            for clone in result['clone_list']:
                clone_count += 1
                for c in clone:
                    row = c.to_dict()
                    row['CLONE_ID'] = clone_count
                    pass_writer.writerow(row)
                    pass_count += 1
        else:
            fail_count += len(result['rec_list'])

        # Write log
        printLog(result['log'], handle=log_handle)
 
    # Print total counts
    printProgress(pass_count + fail_count, result_count, 0.05, start_time) 

    # Update return list
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = pass_count + fail_count
    log['CLONES'] = clone_count
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
                            parents=[getCommonArgParser(seq_in=False, seq_out=False, db_in=True, multiproc=True)], 
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