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
from rpy2.robjects.vectors import StrVector

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
default_distance = 5


def indexJunctions(db_iter, fields=None, mode='allele', action='first', 
                   separator=default_separator):
    """
    Identifies preclonal groups by V, J and junction length

    Arguments: 
    db_iter = an iterator of IgRecords defined by readDbFile
    fields = additional annotation fields to use to group preclones;
             if None use only V, J and junction length
    mode = specificity of alignment call to use for assigning preclones;
           one of ('allele', 'gene')
    action = how to handle multiple value fields when assigning preclones;
             one of ('first', 'set')
    separator = the delimiter separating values within a field 
    
    Returns: 
    a dictionary of {(V, J, junction length):[IgRecords]}
    """
    if mode == 'allele' and fields is None:
        def _get_key(rec, act):
            return (rec.getVAllele(act), rec.getJAllele(act), len(rec.junction))
    elif mode == 'gene' and fields is None:
        def _get_key(rec, act):  
            return (rec.getVGene(act), rec.getJGene(act), len(rec.junction))
    elif mode == 'allele' and fields is not None:
        def _get_key(rec, act):
            vdj = [rec.getVAllele(act), rec.getJAllele(act), len(rec.junction)]
            ann = [rec.toDict().get(k, None) for k in fields]
            return tuple(chain(vdj, ann))
    elif mode == 'gene' and fields is not None:
        def _get_key(rec, act):
            vdj = [rec.getVGene(act), rec.getJGene(act), len(rec.junction)]
            ann = [rec.toDict().get(k, None) for k in fields]
            return tuple(chain(vdj, ann))
 
    clone_index = {}
    for rec in db_iter:
        key = _get_key(rec, action)
        #print key
        # Assigned passed preclone records to key and failed to index None
        if all([k is not None for k in key]):
            clone_index.setdefault(key, []).append(rec)
        else:
            clone_index.setdefault(None, []).append(rec)
    
    if action == 'set':
        print 'WARNING: set mode is unimplemented\n'
    #    for k in clone_index:
    #        pass
        
    return clone_index


def distanceClones(records, distance=default_distance):
    """
    Separates a set of IgRecords into clones

    Arguments: 
    records = an iterator of IgRecords
    distance = the distance threshold to assign clonal groups
    
    Returns: 
    a dictionary of lists defining {clone number: [IgRecords clonal group]}
    """
    # Define unique junction mapping
    junc_map = {}
    for r in records:
        junc_map.setdefault(str(r.junction.upper()), []).append(r)
    
    # Process records    
    if any([len(j) == 0 for j in junc_map]):
        clone_dict = None
    elif len(junc_map) == 1:
        clone_dict = {1:records}
    else:
        # Call R function
        #junctions = '|'.join(junc_map.keys())
        junctions = StrVector(junc_map.keys())
        #print junctions
        clone_list = ClonesByDist.getClones(junctions, distance)
        # Parse R output        
        clone_dict = {}
        for i, c in enumerate(clone_list):
            clone_dict[i + 1] = list(chain(*[junc_map[str(j).upper()] for j in c]))
        # Parse R output
        #clone_dict = {}        
        #clone_split = str(clone_str[0]).split('|')
        #for i, c in enumerate(clone_split):
        #    junc = c.split(',')
            #print list(chain(*[junc_map[j.upper()] for j in junc_r]))
            #clone_list.append(list(chain(*[junc_map[j.upper()] for j in junc_r])))
        #    clone_dict[i + 1] = list(chain(*[junc_map[j.upper()] for j in junc]))
    
    #print clone_dict
    return clone_dict


def feedQueue(data_queue, nproc, db_file, group_func, group_args={}):
    """
    Feeds the data queue with Ig records

    Arguments: 
    data_queue = a multiprocessing.Queue to hold data for processing
    nproc = the number of processQueue processes
    db_file = the Ig record database file
    group_func = the function to use for assigning preclones
    group_args = a dictionary of arguments to pass to group_func
    
    Returns: 
    None
    """
    # Iterate over Ig records and feed data queue
    db_iter = readDbFile(db_file)
    clone_dict = group_func(db_iter, **group_args)
    
    #db_iter = readDbFile(db_file)
    for k, v in clone_dict.iteritems():
        # Feed queue
        group_pass = False if k is None else True
        data_queue.put({'id':k, 'records':v, 'pass':group_pass})
    
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
        records = args['records']
        
        # Define result dictionary for iteration
        results = {'id':args['id'],
                   'records':records,
                   'clones':None,
                   'pass':False,
                   'log':OrderedDict([('ID', args['id'])])}
 
        # Add V(D)J to log
        results['log']['VALLELE'] = ','.join(set([(r.getVAllele() or '') for r in records]))
        results['log']['DALLELE'] = ','.join(set([(r.getDAllele() or '') for r in records]))
        results['log']['JALLELE'] = ','.join(set([(r.getJAllele() or '') for r in records]))
        results['log']['JUNCLEN'] = ','.join(set([(str(len(r.junction)) or '0') for r in records]))
        results['log']['READS'] = len(records)
        
        # Checking for preclone failure and assign clones
        clones = clone_func(records, **clone_args) if args['pass'] else None
        if clones is not None:
            results['clones'] = clones
            results['pass'] = True
            results['log']['CLONES'] = len(clones)
        else:
            results['log']['CLONES'] = 0

        # Feed results to result queue
        result_queue.put(results)

    return None


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
    pass_writer = getDbWriter(pass_handle, db_file, add_fields='CLONE')
    
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
    rec_count = clone_count = pass_count = fail_count = 0
    for result in iter(result_queue.get, None): 
        # Print progress for previous iteration and update record count
        printProgress(rec_count, result_count, 0.05, start_time) 
        rec_count += len(result['records'])
        
        # Write passed and failed records
        if result['pass']:
            for clone in result['clones'].itervalues():
                clone_count += 1
                for i, rec in enumerate(clone):
                    rec.annotations['CLONE'] = clone_count
                    pass_writer.writerow(rec.toDict())
                    pass_count += 1
                    result['log']['CLONE%i-%i' % (clone_count, i + 1)] = str(rec.junction)

        else:
            for i, rec in enumerate(result['records']):
                fail_writer.writerow(rec.toDict())
                fail_count += 1
                result['log']['CLONE0-%i' % (i + 1)] = str(rec.junction)
                
        # Write log
        printLog(result['log'], handle=log_handle)
 
    # Print total counts
    printProgress(rec_count, result_count, 0.05, start_time) 

    # Update return list
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['CLONES'] = clone_count
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


def defineClones(db_file, group_func, clone_func, group_args={}, clone_args={}, 
                 out_args=default_out_args, nproc=None, queue_size=None):
    """
    Define clonally related sequences
    
    Arguments:
    db_file = filename of input database
    group_func = the function to use for assigning preclones
    clone_func = the function to use for determining clones within preclonal groups
    group_args = a dictionary of arguments to pass to group_func
    clone_args = a dictionary of arguments to pass to clone_func
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
    log['GROUP_FUNC'] = group_func.__name__
    log['GROUP_ARGS'] = group_args
    log['CLONE_FUNC'] = clone_func.__name__
    log['CLONE_ARGS'] = clone_args
    log['NPROC'] = nproc
    printLog(log)
    
    # Define shared data objects
    manager = mp.Manager()
    collect_dict = manager.dict()
    data_queue = mp.Queue(queue_size)
    result_queue = mp.Queue(queue_size)

    try:  
        # Initiate feeder process
        feeder = mp.Process(target=feedQueue, args=(data_queue, nproc, db_file, 
                                                    group_func, group_args))
        feeder.start()
    
        # Initiate processQueue processes
        workers = []
        for __ in range(nproc):
            w = mp.Process(target=processQueue, args=(data_queue, result_queue, 
                                                      clone_func, clone_args, 
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
    
    except Exception as e:
        print "EXCEPTION: %s" % e
        feeder.terminate()
        feeder.join()
        for w in workers:
            w.terminate()
            w.join()
        collector.terminate()
        collector.join
        manager.shutdown()
        
        return None


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
                            formatter_class=ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(title='subcommands', dest='command', 
                                       help='Cloning method', metavar='')
    
    # Parent parser    
    parser_parent = getCommonArgParser(seq_in=False, seq_out=False, db_in=True, multiproc=True)
    
    # S5F distance cloning method
    parser_dist = subparsers.add_parser('dist', parents=[parser_parent],
                                        formatter_class=ArgumentDefaultsHelpFormatter,
                                        help='Defines clones V assignment, J assignment and junction \
                                              length with S5F distance model')
    parser_dist.add_argument('-f', nargs='+', action='store', dest='fields', default=None,
                             help='Additional fields to use for grouping clones (non VDJ)')
    parser_dist.add_argument('--mode', action='store', dest='mode', 
                             choices=('allele', 'gene'), default='gene', 
                             help='Specifies whether to use the V(D)J allele or gene for preclone assignment')
    parser_dist.add_argument('--act', action='store', dest='action', default='first',
                             choices=('first', 'set'),
                             help='Specifies how to handle multiple V(D)J assignments for preclone assignment') 
    parser_dist.add_argument('--dist', action='store', dest='distance', type=float, 
                             default=default_distance,
                             help='The junction distance threshold for clonal grouping')
    parser_dist.set_defaults(group_func=indexJunctions)    
    parser_dist.set_defaults(clone_func=distanceClones)
        
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
    if 'fields' in args_dict and args_dict['fields'] is not None:  
        args_dict['fields'] = map(str.upper, args_dict['fields'])
    
    # Define clone_args
    if args.command == 'dist':
        args_dict['group_args'] = {'fields': args_dict['fields'],
                                   'action': args_dict['action'], 
                                   'mode':args_dict['mode']}
        args_dict['clone_args'] = {'distance': args_dict['distance']}
        del args_dict['fields']
        del args_dict['action']
        del args_dict['mode']
        del args_dict['distance']
    
    # Call defineClones
    del args_dict['command']
    del args_dict['db_files']
    for f in args.__dict__['db_files']:
        args_dict['db_file'] = f
        defineClones(**args_dict)