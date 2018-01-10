#!/usr/bin/env python3
"""
Assign Ig sequences into clones
"""
# Info
__author__ = 'Namita Gupta, Jason Anthony Vander Heiden, Gur Yaari, Mohamed Uduman'
from changeo import __version__, __date__

# Imports
import os
import re
import sys
import csv
import numpy as np
from argparse import ArgumentParser
from collections import OrderedDict
from itertools import chain
from textwrap import dedent
from time import time
from Bio import pairwise2
from Bio.Seq import translate

# Presto and changeo imports
from presto.Defaults import default_out_args
from presto.IO import getOutputHandle, printLog, printProgress
from presto.Multiprocessing import manageProcesses
from changeo.Commandline import CommonHelpFormatter, checkArgs, getCommonArgParser, parseCommonArgs
from changeo.Distance import distance_models, calcDistances, formClusters
from changeo.IO import countDbFile, getDbFields
from changeo.Multiprocessing import DbData, DbResult, feedDbQueue
from changeo.Parsers import ChangeoSchema, ChangeoReader, ChangeoWriter

# Defaults
default_translate = False
default_distance = 0.0
default_index_mode = 'gene'
default_index_action = 'set'
default_distance_model = 'ham'
default_seq_field = 'JUNCTION'
default_norm = 'len'
default_sym = 'avg'
default_linkage = 'single'
default_max_missing=0
choices_distance_model = ('ham', 'aa', 'hh_s1f', 'hh_s5f', 'mk_rs1nf', 'mk_rs5nf', 'hs1f_compat', 'm1n_compat')


def filterMissing(data, field=default_seq_field, max_missing=default_max_missing):
    """
    Splits a list of Receptor into passed and failed groups based on the number
    of missing characters in the sequence

    Arguments:
        data : list of Receptor records
        field : sequence field to filter on
        max_missing : maximum number of missing characters (non-ACGT) to permit before failing the record

    Returns:
        dict : a dictionary of {pass : list of passing records, fail : list of failing records}
    """
    # Function to validate the sequence string
    def _pass(seq):
        if len(seq) > 0 and len(re.findall(r'[^ACGT]', seq)) <= max_missing:
            return True
        else:
            return False

    # Define return object
    result = {'pass': None, 'fail': None}

    pass_list = []
    fail_list = []
    for rec in data:
        seq = str(rec.getSeq(field))
        if _pass(seq):  pass_list.append(rec)
        else:  fail_list.append(rec)

    return {'pass': pass_list, 'fail': fail_list}


def indexByIdentity(index, key, rec, fields=None):
    """
    Updates a preclone index with a simple key

    Arguments:
      index : preclone index from groupByGene
      key : index key
      rec : Receptor to add to the index
      fields : additional annotation fields to use to group preclones;
               if None use only V, J and junction length

    Returns:
      None : Updates index with new key and records.
    """
    index.setdefault(tuple(key), []).append(rec)


def indexByUnion(index, key, rec, fields=None):
    """
    Updates a preclone index with the union of nested keys

    Arguments:
      index : preclone index from groupByGene
      key : index key
      rec : Receptor to add to the index
      fields : additional annotation fields to use to group preclones;
               if None use only V, J and junction length

    Returns:
      None : Updates index with new key and records.
    """
    # List of values for this/new key
    val = [rec]
    f_range = list(range(2, 3 + (len(fields) if fields else 0)))

    # See if field/junction length combination exists in index
    outer_dict = index
    for field in f_range:
        try:
            outer_dict = outer_dict[key[field]]
        except (KeyError):
            outer_dict = None
            break
    # If field combination exists, look through Js
    j_matches = []
    if outer_dict is not None:
        for j in outer_dict.keys():
            if not set(key[1]).isdisjoint(set(j)):
                key[1] = tuple(set(key[1]).union(set(j)))
                j_matches += [j]
    # If J overlap exists, look through Vs for each J
    for j in j_matches:
        v_matches = []
        # Collect V matches for this J
        for v in outer_dict[j].keys():
            if not set(key[0]).isdisjoint(set(v)):
                key[0] = tuple(set(key[0]).union(set(v)))
                v_matches += [v]
        # If there are V overlaps for this J, pop them out
        if v_matches:
            val += list(chain(*(outer_dict[j].pop(v) for v in v_matches)))
            # If the J dict is now empty, remove it
            if not outer_dict[j]:
                outer_dict.pop(j, None)

    # Add value(s) into index nested dictionary
    # OMG Python pointers are the best!
    # Add field dictionaries into index
    outer_dict = index
    for field in f_range:
        outer_dict.setdefault(key[field], {})
        outer_dict = outer_dict[key[field]]
    # Add J, then V into index
    if key[1] in outer_dict:
        outer_dict[key[1]].update({key[0]: val})
    else:
        outer_dict[key[1]] = {key[0]: val}


def groupByGene(db_iter, fields=None, mode=default_index_mode,
                   action=default_index_action):
    """
    Identifies preclonal groups by V, J and junction length

    Arguments: 
      db_iter : an iterator of Receptor objects defined by ChangeoReader
      fields : additional annotation fields to use to group preclones;
               if None use only V, J and junction length
      mode : specificity of alignment call to use for assigning preclones;
             one of ('allele', 'gene')
      action : how to handle multiple value fields when assigning preclones;
               one of ('first', 'set')
    
    Returns: 
      dict: dictionary of {(V, J, junction length):[Receptor]}
    """
    # print(fields)
    # Define functions for grouping keys
    if mode == 'allele' and fields is None:
        def _get_key(rec, act):
            return [rec.getVAllele(act), rec.getJAllele(act),
                    None if rec.junction is None else len(rec.junction)]
    elif mode == 'gene' and fields is None:
        def _get_key(rec, act):  
            return [rec.getVGene(act), rec.getJGene(act),
                    None if rec.junction is None else len(rec.junction)]
    elif mode == 'allele' and fields is not None:
        def _get_key(rec, act):
            vdj = [rec.getVAllele(act), rec.getJAllele(act),
                    None if rec.junction is None else len(rec.junction)]
            ann = [rec.getChangeo(k) for k in fields]
            return list(chain(vdj, ann))
    elif mode == 'gene' and fields is not None:
        def _get_key(rec, act):
            vdj = [rec.getVGene(act), rec.getJGene(act),
                    None if rec.junction is None else len(rec.junction)]
            ann = [rec.getChangeo(k) for k in fields]
            return list(chain(vdj, ann))

    # Function to flatten nested dictionary
    def _flatten_dict(d, parent_key=''):
        items = []
        for k, v in d.items():
            new_key = parent_key + [k] if parent_key else [k]
            if isinstance(v, dict):
                items.extend(_flatten_dict(v, new_key).items())
            else:
                items.append((new_key, v))
        flat_dict = {None if None in i[0] else tuple(i[0]): i[1] for i in items}
        return flat_dict

    if action == 'first':
        index_func = indexByIdentity
    elif action == 'set':
        index_func = indexByUnion
    else:
        sys.stderr.write('Unrecognized action: %s.\n' % action)

    start_time = time()
    clone_index = {}
    rec_count = 0
    for rec in db_iter:
        key = _get_key(rec, action)

        # Print progress
        if rec_count == 0:
            print('PROGRESS> Grouping sequences')

        printProgress(rec_count, step=1000, start_time=start_time)
        rec_count += 1

        # Assigned passed preclone records to key and failed to index None
        if all([k is not None and k != '' for k in key]):
            # Update index dictionary
            index_func(clone_index, key, rec, fields)
        else:
            clone_index.setdefault(None, []).append(rec)

    printProgress(rec_count, step=1000, start_time=start_time, end=True)

    if action == 'set':
        clone_index = _flatten_dict(clone_index)

    return clone_index


def distanceClones(records, model=default_distance_model, distance=default_distance,
                   dist_mat=None, norm=default_norm, sym=default_sym,
                   linkage=default_linkage, seq_field=default_seq_field):
    """
    Separates a set of Receptor objects into clones

    Arguments: 
      records : an iterator of Receptor objects
      model : substitution model used to calculate distance
      distance : the distance threshold to assign clonal groups
      dist_mat : pandas DataFrame of pairwise nucleotide or amino acid distances
      norm : normalization method
      sym : symmetry method
      linkage : type of linkage
      seq_field : sequence field used to calculate distance between records

    Returns: 
      dict : dictionary of lists defining {clone number: [Receptor objects clonal group]}
    """
    # Get distance matrix if not provided
    if dist_mat is None:
        try:
            dist_mat = distance_models[model]
        except KeyError:
            sys.exit('Unrecognized distance model: %s' % args_dict['model'])

    # TODO:  can be cleaned up with abstract model class
    # Determine length of n-mers
    if model in ['hs1f_compat', 'm1n_compat', 'aa', 'ham', 'hh_s1f', 'mk_rs1nf']:
        nmer_len = 1
    elif model in ['hh_s5f', 'mk_rs5nf']:
        nmer_len = 5
    else:
        sys.exit('Unrecognized distance model: %s.\n' % model)

    # Define unique junction mapping
    seq_map = {}
    for rec in records:
        seq = rec.getChangeo(seq_field, seq=True)
        # Check if sequence length is 0
        if len(seq) == 0:
            return None

        seq = re.sub('[\.-]', 'N', str(seq))
        if model == 'aa':  seq = translate(seq)

        seq_map.setdefault(seq, []).append(rec)

    # Process records
    if len(seq_map) == 1:
        return {1:records}

    # Define sequences
    seqs = list(seq_map.keys())

    # Calculate pairwise distance matrix
    dists = calcDistances(seqs, nmer_len, dist_mat, sym=sym, norm=norm)

    # Perform hierarchical clustering
    clusters = formClusters(dists, linkage, distance)

    # Turn clusters into clone dictionary
    clone_dict = {}
    for i, c in enumerate(clusters):
        clone_dict.setdefault(c, []).extend(seq_map[seqs[i]])

    return clone_dict


def processQueue(alive, data_queue, result_queue, max_missing=default_max_missing,
                 clone_func=distanceClones, clone_args={}):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments: 
      alive : a multiprocessing.Value boolean controlling whether processing continues
              if False exit process
      data_queue : a multiprocessing.Queue holding data to process
      result_queue : a multiprocessing.Queue to hold processed results
      max_missing : maximum number of non-ACGT characters to allow in the junction sequence.
      clone_func : the function to call for clonal assignment
      clone_args : a dictionary of arguments to pass to clone_func

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

            # Define result object for iteration and get data records
            result = DbResult(data.id, data.data)

            # Check for invalid data (due to failed indexing) and add failed result
            if not data:
                result_queue.put(result)
                continue

            # Filter records based on missing content
            seq_field = clone_args['seq_field'] if 'seq_field' in clone_args else 'JUNCTION'
            filtered = filterMissing(data.data, field=seq_field, max_missing=max_missing)
            records = filtered['pass']
            result.failed = filtered['fail']

            # Add V(D)J to log
            result.log['ID'] = ','.join([str(x) for x in data.id])
            result.log['VALLELE'] = ','.join(set([(r.getVAllele() or '') for r in data.data]))
            result.log['DALLELE'] = ','.join(set([(r.getDAllele() or '') for r in data.data]))
            result.log['JALLELE'] = ','.join(set([(r.getJAllele() or '') for r in data.data]))
            result.log['JUNCLEN'] = ','.join(set([(str(len(r.junction)) or '0') for r in data.data]))
            result.log['PASSCOUNT'] = len(records)
            result.log['FAILCOUNT'] = len(result.failed)
             
            # Checking for preclone failure and assign clones
            clones = clone_func(records, **clone_args) if records else None

            # import cProfile
            # prof = cProfile.Profile()
            # clones = prof.runcall(clone_func, records, **clone_args)
            # prof.dump_stats('worker-%d.prof' % os.getpid())

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


def collectQueue(alive, result_queue, collect_queue, db_file, out_args=default_out_args):
    """
    Assembles results from a queue of individual sequence results and manages log/file I/O

    Arguments: 
      alive = a multiprocessing.Value boolean controlling whether processing continues
              if False exit process
      result_queue : a multiprocessing.Queue holding processQueue results
      collect_queue : a multiprocessing.Queue to store collector return values
      db_file : the input database file name
      out_args : common output argument dictionary from parseCommonArgs
    
    Returns:
      None
    """
    # Open output files
    try:
        # Count records and define output format
        result_count = countDbFile(db_file)
        out_fields = getDbFields(db_file, add='CLONE')

        # Defined successful output handle
        pass_handle = getOutputHandle(db_file, 
                                      out_label='clone-pass', 
                                      out_dir=out_args['out_dir'], 
                                      out_name=out_args['out_name'], 
                                      out_type='tsv')
        pass_writer = ChangeoWriter(pass_handle, fields=out_fields)

        # Defined failed alignment output handle
        if out_args['failed']:
            fail_handle = getOutputHandle(db_file,
                                          out_label='clone-fail', 
                                          out_dir=out_args['out_dir'], 
                                          out_name=out_args['out_name'], 
                                          out_type='tsv')
            fail_writer = ChangeoWriter(fail_handle, fields=out_fields)
        else:
            fail_handle = None
            fail_writer = None

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
            if rec_count == 0:
                print('PROGRESS> Assigning clones')
            printProgress(rec_count, result_count, 0.05, start_time) 
            rec_count += len(result.data)
            
            # Write passed and failed records
            if result:
                # Writing passing sequences
                for clone in result.results.values():
                    clone_count += 1
                    for i, rec in enumerate(clone, start=1):
                        rec.annotations['clone'] = clone_count
                        pass_writer.writeReceptor(rec)
                        pass_count += 1
                        result.log['CLONE%i-%i' % (clone_count, i)] = str(rec.junction)
                # Write failed sequences from passing sets
                if result.failed:
                    for i, rec in enumerate(result.failed, start=1):
                        fail_count += 1
                        if fail_writer is not None: fail_writer.writeReceptor(rec)
                        result.log['FAIL%i-%i' % (clone_count, i)] = str(rec.junction)
            else:
                for i, rec in enumerate(result.data, start=1):
                    if fail_writer is not None: fail_writer.writeReceptor(rec)
                    fail_count += 1
                    result.log['CLONE0-%i' % (i)] = str(rec.junction)
                    
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
        collect_dict = {'log':log, 'out_files': [pass_handle.name]}
        collect_queue.put(collect_dict)
    except:
        #sys.stderr.write('Exception in collector result processing step\n')
        alive.value = False
        raise

    return None


def defineClones(db_file, group_func=groupByGene, group_args={},
                 clone_func=distanceClones, clone_args={},
                 max_missing=default_max_missing,
                 out_args=default_out_args, nproc=None, queue_size=None):
    """
    Define clonally related sequences
    
    Arguments:
      db_file : filename of input database
      group_func : the function to use for assigning preclones
      group_args : a dictionary of arguments to pass to group_func
      clone_func : the function to use for determining clones within preclonal groups
      clone_args : a dictionary of arguments to pass to clone_func
      max_missing : maximum number of non-ACGT characters to allow in the junction sequence.
      out_args : common output argument dictionary from parseCommonArgs
      nproc : the number of processQueue processes;
              if None defaults to the number of CPUs
      queue_size : maximum size of the argument queue;
                   if None defaults to 2*nproc
    
    Returns:
      list : successful output file names
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'DefineClones'
    log['DB_FILE'] = os.path.basename(db_file)
    log['MAX_MISSING'] = max_missing
    log['GROUP_FUNC'] = group_func.__name__
    log['GROUP_ARGS'] = group_args
    log['CLONE_FUNC'] = clone_func.__name__

    # TODO:  this is yucky, but can be fixed by using a model class
    clone_log = clone_args.copy()
    if 'dist_mat' in clone_log:  del clone_log['dist_mat']
    log['CLONE_ARGS'] = clone_log

    log['NPROC'] = nproc
    printLog(log)
    
    # Define feeder function and arguments
    feed_args = {'db_file': db_file,
                 'group_func': group_func, 
                 'group_args': group_args}
    # Define worker function and arguments
    work_args = {'max_missing': max_missing,
                 'clone_func': clone_func,
                 'clone_args': clone_args}
    # Define collector function and arguments
    collect_args = {'db_file': db_file,
                    'out_args': out_args}

    # Call process manager
    result = manageProcesses(feed_func=feedDbQueue, work_func=processQueue, collect_func=collectQueue,
                             feed_args=feed_args, work_args=work_args, collect_args=collect_args,
                             nproc=nproc, queue_size=queue_size)
        
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
    # Define input and output fields
    fields = dedent(
             '''
             output files:
                 clone-pass
                     database with assigned clonal group numbers.
                 clone-fail
                     database with records failing clonal grouping.

             required fields:
                 SEQUENCE_ID, V_CALL or V_CALL_GENOTYPED, D_CALL, J_CALL, JUNCTION

                 <field>
                     sequence field specified by the --sf parameter
                
             output fields:
                 CLONE
              ''')

    # Parent parser
    parser_parent = getCommonArgParser(seq_in=False, seq_out=False, db_in=True,
                                       multiproc=True)
    # Define argument parser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            parents=[parser_parent],
                            formatter_class=CommonHelpFormatter)
    parser.add_argument('--version', action='version',
                        version='%(prog)s:' + ' %s-%s' %(__version__, __date__))

    # Distance cloning method
    parser.add_argument('-f', nargs='+', action='store', dest='fields', default=None,
                             help='Additional fields to use for grouping clones (non VDJ)')
    parser.add_argument('--mode', action='store', dest='mode',
                             choices=('allele', 'gene'), default=default_index_mode,
                             help='''Specifies whether to use the V(D)J allele or gene for
                                  initial grouping.''')
    parser.add_argument('--act', action='store', dest='action',
                             choices=('first', 'set'), default=default_index_action,
                             help='''Specifies how to handle multiple V(D)J assignments
                                  for initial grouping.''')
    parser.add_argument('--model', action='store', dest='model',
                             choices=choices_distance_model,
                             default=default_distance_model,
                             help='''Specifies which substitution model to use for calculating distance
                                  between sequences. The "ham" model is nucleotide Hamming distance and
                                  "aa" is amino acid Hamming distance. The "hh_s1f" and "hh_s5f" models are
                                  human specific single nucleotide and 5-mer content models, respectively,
                                  from Yaari et al, 2013. The "mk_rs1nf" and "mk_rs5nf" models are
                                  mouse specific single nucleotide and 5-mer content models, respectively,
                                  from Cui et al, 2016. The "m1n_compat" and "hs1f_compat" models are
                                  deprecated models provided backwards compatibility with the "m1n" and
                                  "hs1f" models in Change-O v0.3.3 and SHazaM v0.1.4. Both
                                  5-mer models should be considered experimental.''')
    parser.add_argument('--dist', action='store', dest='distance', type=float,
                             default=default_distance,
                             help='The distance threshold for clonal grouping')
    parser.add_argument('--norm', action='store', dest='norm',
                             choices=('len', 'mut', 'none'), default=default_norm,
                             help='''Specifies how to normalize distances. One of none
                                  (do not normalize), len (normalize by length),
                                  or mut (normalize by number of mutations between sequences).''')
    parser.add_argument('--sym', action='store', dest='sym',
                             choices=('avg', 'min'), default=default_sym,
                             help='''Specifies how to combine asymmetric distances. One of avg
                                  (average of A->B and B->A) or min (minimum of A->B and B->A).''')
    parser.add_argument('--link', action='store', dest='linkage',
                             choices=('single', 'average', 'complete'), default=default_linkage,
                             help='''Type of linkage to use for hierarchical clustering.''')
    parser.add_argument('--maxmiss', action='store', dest='max_missing', type=int,
                                default=default_max_missing,
                                help='''The maximum number of non-ACGT characters (gaps or Ns) to 
                                     permit in the junction sequence before excluding the record 
                                     from clonal assignment. Warning, under single linkage 
                                     non-informative positions can create artifactual links 
                                     between unrelated sequences. Use with caution.''')
    parser.add_argument('--sf', action='store', dest='seq_field',
                                default=default_seq_field,
                                help='''The name of the field to be used to calculate
                                     distance between records''')
    parser.set_defaults(group_func=groupByGene)
    parser.set_defaults(clone_func=distanceClones)
        
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
    if 'seq_field' in args_dict:
        args_dict['seq_field'] = args_dict['seq_field'].upper()
    if 'fields' in args_dict and args_dict['fields'] is not None:  
        args_dict['fields'] = [f.upper() for f in args_dict['fields']]
    
    # Define grouping and cloning function arguments
    args_dict['group_args'] = {'fields': args_dict['fields'],
                               'action': args_dict['action'],
                               'mode':args_dict['mode']}
    args_dict['clone_args'] = {'model':  args_dict['model'],
                               'distance':  args_dict['distance'],
                               'norm': args_dict['norm'],
                               'sym': args_dict['sym'],
                               'linkage': args_dict['linkage'],
                               'seq_field': args_dict['seq_field']}

    # Get distance matrix
    try:
        args_dict['clone_args']['dist_mat'] = distance_models[args_dict['model']]
    except KeyError:
        sys.exit('Unrecognized distance model: %s' % args_dict['model'])

    # Clean argument dictionary
    del args_dict['fields']
    del args_dict['action']
    del args_dict['mode']
    del args_dict['model']
    del args_dict['distance']
    del args_dict['norm']
    del args_dict['sym']
    del args_dict['linkage']
    del args_dict['seq_field']

    # Call defineClones
    del args_dict['db_files']
    for f in args.__dict__['db_files']:
        args_dict['db_file'] = f
        defineClones(**args_dict)
