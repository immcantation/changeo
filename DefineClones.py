#!/usr/bin/env python
"""
Assign Ig sequences into clones
"""

__author__    = 'Jason Anthony Vander Heiden, Namita Gupta, Gur Yaari, Mohamed Uduman'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.0'
__date__      = '2014.4.13'

# Imports
import os, signal, sys
import multiprocessing as mp
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import OrderedDict
from ctypes import c_bool
from itertools import chain
from time import time

# IgCore imports
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from IgCore import default_separator, default_out_args
from IgCore import getCommonArgParser, parseCommonArgs
from IgCore import getFileType, getOutputHandle, printLog, printProgress
from DbCore import countDbFile, readDbFile, getDbWriter, IgRecord

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


class DbData:
    """
    A class defining IgRecord data objects for worker processes
    """
    # Instantiation
    def __init__(self, key, records):
        self.id = key
        self.data = records
        self.valid = (key is not None and records is not None)
    
    # Boolean evaluation
    def __nonzero__(self): 
        return self.valid
    
    # Length evaluation
    def __len__(self):
        if isinstance(self.data, IgRecord):
            return 1
        elif self.data is None:
            return 0
        else:
            return len(self.data)


class DbResult:
    """
    A class defining IgRecord result objects for collector processes
    """
    # Instantiation
    def __init__(self, key, records):
        self.id = key
        self.data = records
        self.results = None
        self.valid = False
        self.log = OrderedDict([('ID', key)])
        #if isinstance(values, list):
        #    for v in values:  setattr(self, v, None)
        #else:
        #    setattr(self, values, None)
            
    # Boolean evaluation
    def __nonzero__(self): 
        return self.valid
    
    # Length evaluation
    def __len__(self):
        if isinstance(self.results, IgRecord):
            return 1
        elif self.data is None:
            return 0
        else:
            return len(self.results)


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
    # Define functions for grouping keys
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


def feedQueue(alive, data_queue, db_file, group_func, group_args={}):
    """
    Feeds the data queue with Ig records

    Arguments: 
    alive = a multiprocessing.Value boolean controlling whether processing continues
            if False exit process
    data_queue = a multiprocessing.Queue to hold data for processing
    db_file = the Ig record database file
    group_func = the function to use for assigning preclones
    group_args = a dictionary of arguments to pass to group_func
    
    Returns: 
    None
    """
    # Open input file and perform grouping
    try:
        # Iterate over Ig records and assign groups
        db_iter = readDbFile(db_file)
        clone_dict = group_func(db_iter, **group_args)
    except:
        #sys.stderr.write('Exception in feeder grouping step\n')
        alive.value = False
        raise
    
    # Add groups to data queue
    try:
        #print 'START FEED', alive.value
        # Iterate over groups and feed data queue
        clone_iter = clone_dict.iteritems()
        while alive.value:
            # Get data from queue
            if data_queue.full():  continue
            else:  data = next(clone_iter, None)
            # Exit upon reaching end of iterator
            if data is None:  break
            #print "FEED", alive.value, k
            
            # Feed queue
            data_queue.put(DbData(*data))
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
    except:
        #sys.stderr.write('Exception in feeder queue feeding step\n')
        alive.value = False
        raise

    return None


def processQueue(alive, data_queue, result_queue, clone_func, clone_args):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments: 
    alive = a multiprocessing.Value boolean controlling whether processing continues
            if False exit process
    data_queue = a multiprocessing.Queue holding data to process
    result_queue = a multiprocessing.Queue to hold processed results
    clone_func = the function to call for clonal assignment
    clone_args = a dictionary of arguments to pass to clone_func

    Returns: 
    None
    """
    try:
        #print 'START WORK', alive.value
        # Iterator over data queue until sentinel object reached
        while alive.value:
            # Get data from queue
            if data_queue.empty():  continue
            else:  data = data_queue.get()
            # Exit upon reaching sentinel
            if data is None:  break
            #print "WORK", alive.value, data['id']

            # Define result object for iteration and get data records
            records = data.data
            result = DbResult(data.id, records)
             
            # Add V(D)J to log
            result.log['VALLELE'] = ','.join(set([(r.getVAllele() or '') for r in records]))
            result.log['DALLELE'] = ','.join(set([(r.getDAllele() or '') for r in records]))
            result.log['JALLELE'] = ','.join(set([(r.getJAllele() or '') for r in records]))
            result.log['JUNCLEN'] = ','.join(set([(str(len(r.junction)) or '0') for r in records]))
            result.log['READS'] = len(records)
             
            # Checking for preclone failure and assign clones
            clones = clone_func(records, **clone_args) if data else None
            if clones is not None:
                result.results = clones
                result.valid = True
                result.log['CLONES'] = len(clones)
            else:
                result.log['CLONES'] = 0
  
            # Feed results to result queue
            result_queue.put(result)
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
    except:
        #sys.stderr.write('Exception in worker\n')
        alive.value = False
        raise
    
    return None


def collectQueue(alive, result_queue, collect_dict, db_file, out_args):
    """
    Assembles results from a queue of individual sequence results and manages log/file I/O

    Arguments: 
    alive = a multiprocessing.Value boolean controlling whether processing continues
            if False exit process
    result_queue = a multiprocessing.Queue holding processQueue results
    result_count = the total number of results expected
    collect_dict = a multiprocessing.Manager.dict to store return values
    db_file = the input database file name
    out_args = common output argument dictionary from parseCommonArgs
    
    Returns: 
    None
    (adds 'log' and 'out_files' to collect_dict)
    """
    # Open output files
    try:
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
    except:
        #sys.stderr.write('Exception in collector file opening step\n')
        alive.value = False
        raise

    # Get results from queue and write to files
    try:
        #print 'START COLLECT', alive.value
        # Iterator over results queue until sentinel object reached
        start_time = time()
        rec_count = clone_count = pass_count = fail_count = 0
        while alive.value:
            # Get result from queue
            if result_queue.empty():  continue
            else:  result = result_queue.get()
            # Exit upon reaching sentinel
            if result is None:  break
            #print "COLLECT", alive.value, result['id']
            
            # Print progress for previous iteration and update record count
            printProgress(rec_count, result_count, 0.05, start_time) 
            rec_count += len(result.data)
            
            # Write passed and failed records
            if result:
                for clone in result.results.itervalues():
                    clone_count += 1
                    for i, rec in enumerate(clone):
                        rec.annotations['CLONE'] = clone_count
                        pass_writer.writerow(rec.toDict())
                        pass_count += 1
                        result.log['CLONE%i-%i' % (clone_count, i + 1)] = str(rec.junction)
    
            else:
                for i, rec in enumerate(result.data):
                    fail_writer.writerow(rec.toDict())
                    fail_count += 1
                    result.log['CLONE0-%i' % (i + 1)] = str(rec.junction)
                    
            # Write log
            printLog(result.log, handle=log_handle)
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
        
        # Print total counts
        printProgress(rec_count, result_count, 0.05, start_time)

        # Close file handles
        pass_handle.close()
        if fail_handle is not None:  fail_handle.close()
        if log_handle is not None:  log_handle.close()
                
        # Update return list
        log = OrderedDict()
        log['OUTPUT'] = os.path.basename(pass_handle.name)
        log['CLONES'] = clone_count
        log['RECORDS'] = rec_count
        log['PASS'] = pass_count
        log['FAIL'] = fail_count
        collect_dict['log'] = log
        collect_dict['out_files'] = [pass_handle.name]
    except:
        #sys.stderr.write('Exception in collector result processing step\n')
        alive.value = False
        raise

    return None


def manageProcesses(feed_func, work_func, collect_func, 
                    feed_args={}, work_args={}, collect_args={}, 
                    nproc=None, queue_size=None):
    """
    Manages feeder, worker and collector processes
    
    Arguments:
    feed_func = the data Queue feeder function
    work_func = the worker function
    collect_func = the result Queue collector function
    feed_args = a dictionary of arguments to pass to feed_func
    work_args = a dictionary of arguments to pass to work_func
    collect_args = a dictionary of arguments to pass to collect_func
    nproc = the number of processQueue processes;
            if None defaults to the number of CPUs
    queue_size = maximum size of the argument queue;
                 if None defaults to 2*nproc    
    
    Returns:
    a dictionary of collector results
    """
    # Raise KeyboardInterrupt
    def _signalHandler(s, f):
        raise SystemExit

    # Terminate processes
    def _terminate():
        sys.stderr.write('Terminating child processes...')
        # Terminate feeders
        feeder.terminate()
        feeder.join()
        # Terminate workers
        for w in workers:
            w.terminate()
            w.join()
        # Terminate collector
        collector.terminate()
        collector.join
        # Shutdown manager
        manager.shutdown()
        sys.stderr.write('  Done.\n')
    
    # Raise SystemExit upon termination signal
    signal.signal(signal.SIGTERM, _signalHandler)
        
    # Define number of processes and queue size
    if nproc is None:  nproc = mp.cpu_count()
    if queue_size is None:  queue_size = nproc * 2
    
    # Define shared child process keep alive flag
    alive = mp.Value(c_bool, True)
    
    # Initiate manager and define shared data objects
    manager = mp.Manager()
    data_queue = manager.Queue(queue_size)
    result_queue = manager.Queue(queue_size)
    collect_dict = manager.dict()

    try:  
        # Initiate feeder process
        feeder = mp.Process(target=feed_func, 
                            args=(alive, data_queue), 
                            kwargs=feed_args)
        #feeder.daemon = True
        feeder.start()
    
        # Initiate processQueue processes
        workers = []
        for __ in range(nproc):
            w = mp.Process(target=work_func, 
                           args=(alive, data_queue, result_queue), 
                           kwargs=work_args) 
            #w.daemon = True
            w.start()
            workers.append(w)
    
        # Initiate collector process
        collector = mp.Process(target=collect_func, 
                               args=(alive, result_queue, collect_dict), 
                               kwargs=collect_args)
        #collector.daemon = True
        collector.start()
    
        # Wait for feeder to finish
        feeder.join()
        # Add sentinel object to data queue for each worker process
        for __ in range(nproc):  data_queue.put(None)
        
        # Wait for worker processes to finish
        for w in workers:  w.join()
        # Add sentinel to result queue
        result_queue.put(None)
        
        # Wait for collector process to finish
        collector.join()

        # Copy collector results and shutdown manager
        result = dict(collect_dict)
        manager.shutdown()
    except (KeyboardInterrupt, SystemExit):
        sys.stderr.write('Exit signal received\n')
        _terminate()
        sys.exit()
    except Exception as e:
        sys.stderr.write('Error:  %s\n' % e)
        _terminate()
        sys.exit()
    except:
        sys.stderr.write('Error:  Exiting with unknown exception\n')
        _terminate()
        sys.exit()
    else:
        if not alive.value:
            sys.stderr.write('Error:  Exiting due to child process error\n')
            _terminate()
            sys.exit()
    
    return result


    
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
    
    # Define feeder function and arguments
    feed_func = feedQueue
    feed_args = {'db_file': db_file,
                 'group_func': group_func, 
                 'group_args': group_args}
    # Define worker function and arguments
    work_func = processQueue
    work_args = {'clone_func': clone_func, 
                 'clone_args': clone_args}
    # Define collector function and arguments
    collect_func = collectQueue
    collect_args = {'db_file': db_file,
                    'out_args': out_args}
    
    # Call process manager
    result = manageProcesses(feed_func, work_func, collect_func, 
                             feed_args, work_args, collect_args, 
                             nproc, queue_size)
        
    # Print log
    result['log']['END'] = 'DefineClones'
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
                                              length with old substitution distance model')
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