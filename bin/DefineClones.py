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
from presto.IO import getFileType, getOutputHandle, printLog, printProgress
from presto.Multiprocessing import manageProcesses
from presto.Sequence import getDNAScoreDict
from changeo.Commandline import CommonHelpFormatter, checkArgs, getCommonArgParser, parseCommonArgs
from changeo.Distance import distance_models, calcDistances, formClusters
from changeo.IO import getDbWriter, readDbFile, countDbFile
from changeo.Multiprocessing import DbData, DbResult

# Defaults
default_translate = False
default_distance = 0.0
default_index_mode = 'gene'
default_index_action = 'set'
default_bygroup_model = 'ham'
default_hclust_model = 'chen2010'
default_seq_field = 'JUNCTION'
default_norm = 'len'
default_sym = 'avg'
default_linkage = 'single'
choices_bygroup_model = ('ham', 'aa', 'hh_s1f', 'hh_s5f', 'mk_rs1nf', 'mk_rs5nf', 'hs1f_compat', 'm1n_compat')


def indexByIdentity(index, key, rec, fields=None):
    """
    Updates a preclone index with a simple key

    Arguments:
    index = preclone index from indexJunctions
    key = index key
    rec = IgRecord to add to the index
    fields = additional annotation fields to use to group preclones;
             if None use only V, J and junction length

    Returns:
    None. Updates index with new key and records.
    """
    index.setdefault(tuple(key), []).append(rec)


def indexByUnion(index, key, rec, fields=None):
    """
    Updates a preclone index with the union of nested keys

    Arguments:
    index = preclone index from indexJunctions
    key = index key
    rec = IgRecord to add to the index
    fields = additional annotation fields to use to group preclones;
             if None use only V, J and junction length

    Returns:
    None. Updates index with new key and records.
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


def indexJunctions(db_iter, fields=None, mode=default_index_mode,
                   action=default_index_action):
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
    
    Returns: 
    a dictionary of {(V, J, junction length):[IgRecords]}
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
            ann = [rec.toDict().get(k, None) for k in fields]
            return list(chain(vdj, ann))
    elif mode == 'gene' and fields is not None:
        def _get_key(rec, act):
            vdj = [rec.getVGene(act), rec.getJGene(act),
                    None if rec.junction is None else len(rec.junction)]
            ann = [rec.toDict().get(k, None) for k in fields]
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


def distanceClones(records, model=default_bygroup_model, distance=default_distance,
                   dist_mat=None, norm=default_norm, sym=default_sym,
                   linkage=default_linkage, seq_field=default_seq_field):
    """
    Separates a set of IgRecords into clones

    Arguments: 
    records = an iterator of IgRecords
    model = substitution model used to calculate distance
    distance = the distance threshold to assign clonal groups
    dist_mat = pandas DataFrame of pairwise nucleotide or amino acid distances
    norm = normalization method
    sym = symmetry method
    linkage = type of linkage
    seq_field = sequence field used to calculate distance between records

    Returns: 
    a dictionary of lists defining {clone number: [IgRecords clonal group]}
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
    for ig in records:
        seq = ig.getSeqField(seq_field)
        # Check if sequence length is 0
        if len(seq) == 0:
            return None

        seq = re.sub('[\.-]', 'N', str(seq))
        if model == 'aa':  seq = translate(seq)

        seq_map.setdefault(seq, []).append(ig)

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


def distChen2010(records):
    """
    Calculate pairwise distances as defined in Chen 2010
    
    Arguments:
    records = list of IgRecords where first is query to be compared to others in list
    
    Returns:
    list of distances
    """
    # Pull out query sequence and V/J information
    query = records.popitem(last=False)
    query_cdr3 = query.junction[3:-3]
    query_v_allele = query.getVAllele()
    query_v_gene = query.getVGene()
    query_v_family = query.getVFamily()
    query_j_allele = query.getJAllele()
    query_j_gene = query.getJGene()
    # Create alignment scoring dictionary
    score_dict = getDNAScoreDict()
    
    scores = [0]*len(records)    
    for i in range(len(records)):
        ld = pairwise2.align.globalds(query_cdr3, records[i].junction[3:-3],
                                      score_dict, -1, -1, one_alignment_only=True)
        # Check V similarity
        if records[i].getVAllele() == query_v_allele: ld += 0
        elif records[i].getVGene() == query_v_gene: ld += 1
        elif records[i].getVFamily() == query_v_family: ld += 3
        else: ld += 5
        # Check J similarity
        if records[i].getJAllele() == query_j_allele: ld += 0
        elif records[i].getJGene() == query_j_gene: ld += 1
        else: ld += 3
        # Divide by length
        scores[i] = ld/max(len(records[i].junction[3:-3]), query_cdr3)
        
    return scores


def distAdemokun2011(records):
    """
    Calculate pairwise distances as defined in Ademokun 2011
    
    Arguments:
    records = list of IgRecords where first is query to be compared to others in list
    
    Returns:
    list of distances
    """
    # Pull out query sequence and V family information
    query = records.popitem(last=False)
    query_cdr3 = query.junction[3:-3]
    query_v_family = query.getVFamily()
    # Create alignment scoring dictionary
    score_dict = getDNAScoreDict()
    
    scores = [0]*len(records)    
    for i in range(len(records)):
        
        if abs(len(query_cdr3) - len(records[i].junction[3:-3])) > 10:
            scores[i] = 1
        elif query_v_family != records[i].getVFamily(): 
            scores[i] = 1
        else: 
            ld = pairwise2.align.globalds(query_cdr3, records[i].junction[3:-3], 
                                          score_dict, -1, -1, one_alignment_only=True)
            scores[i] = ld/min(len(records[i].junction[3:-3]), query_cdr3)
    
    return scores


def hierClust(dist_mat, method='chen2010'):
    """
    Calculate hierarchical clustering
    
    Arguments:
    dist_mat = square-formed distance matrix of pairwise CDR3 comparisons
    
    Returns:
    list of cluster ids
    """
    if method == 'chen2010':
        clusters = formClusters(dist_mat, 'average', 0.32)
    elif method == 'ademokun2011':
        clusters = formClusters(dist_mat, 'complete', 0.25)
    else: clusters = np.ones(dist_mat.shape[0])
        
    return clusters

# TODO:  Merge duplicate feed, process and collect functions.
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
        clone_iter = iter(clone_dict.items())
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


def feedQueueClust(alive, data_queue, db_file, group_func=None, group_args={}):
    """
    Feeds the data queue with Ig records

    Arguments: 
    alive = a multiprocessing.Value boolean controlling whether processing continues
            if False exit process
    data_queue = a multiprocessing.Queue to hold data for processing
    db_file = the Ig record database file
    
    Returns: 
    None
    """
    # Open input file and perform grouping
    try:
        # Iterate over Ig records and order by junction length
        records = {}
        db_iter = readDbFile(db_file)
        for rec in db_iter:
            records[rec.id] = rec
        records = OrderedDict(sorted(list(records.items()), key=lambda i: i[1].junction_length))
        dist_dict = {}
        for __ in range(len(records)):
            k,v = records.popitem(last=False)
            dist_dict[k] = [v].append(list(records.values()))
    except:
        #sys.stderr.write('Exception in feeder grouping step\n')
        alive.value = False
        raise
    
    # Add groups to data queue
    try:
        # print 'START FEED', alive.value
        # Iterate over groups and feed data queue
        dist_iter = iter(dist_dict.items())
        while alive.value:
            # Get data from queue
            if data_queue.full():  continue
            else:  data = next(dist_iter, None)
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
        # Iterator over data queue until sentinel object reached
        while alive.value:
            # Get data from queue
            if data_queue.empty():  continue
            else:  data = data_queue.get()
            # Exit upon reaching sentinel
            if data is None:  break

            # Define result object for iteration and get data records
            records = data.data
            # print(data.id)
            result = DbResult(data.id, records)

            # Check for invalid data (due to failed indexing) and add failed result
            if not data:
                result_queue.put(result)
                continue

            # Add V(D)J to log
            result.log['ID'] = ','.join([str(x) for x in data.id])
            result.log['VALLELE'] = ','.join(set([(r.getVAllele() or '') for r in records]))
            result.log['DALLELE'] = ','.join(set([(r.getDAllele() or '') for r in records]))
            result.log['JALLELE'] = ','.join(set([(r.getJAllele() or '') for r in records]))
            result.log['JUNCLEN'] = ','.join(set([(str(len(r.junction)) or '0') for r in records]))
            result.log['SEQUENCES'] = len(records)
             
            # Checking for preclone failure and assign clones
            clones = clone_func(records, **clone_args) if data else None

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


def processQueueClust(alive, data_queue, result_queue, clone_func, clone_args):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments: 
    alive = a multiprocessing.Value boolean controlling whether processing continues
            if False exit process
    data_queue = a multiprocessing.Queue holding data to process
    result_queue = a multiprocessing.Queue to hold processed results
    clone_func = the function to call for calculating pairwise distances between sequences
    clone_args = a dictionary of arguments to pass to clone_func

    Returns: 
    None
    """
    
    try:
        # print 'START WORK', alive.value
        # Iterator over data queue until sentinel object reached
        while alive.value:
            # Get data from queue
            if data_queue.empty():  continue
            else:  data = data_queue.get()
            # Exit upon reaching sentinel
            if data is None:  break
            # print "WORK", alive.value, data['id']

            # Define result object for iteration and get data records
            records = data.data
            result = DbResult(data.id, records)
             
            # Create row of distance matrix and check for error
            dist_row = clone_func(records, **clone_args) if data else None
            if dist_row is not None:
                result.results = dist_row
                result.valid = True
  
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


def collectQueue(alive, result_queue, collect_queue, db_file, out_args, cluster_func=None, cluster_args={}):
    """
    Assembles results from a queue of individual sequence results and manages log/file I/O

    Arguments: 
    alive = a multiprocessing.Value boolean controlling whether processing continues
            if False exit process
    result_queue = a multiprocessing.Queue holding processQueue results
    collect_queue = a multiprocessing.Queue to store collector return values
    db_file = the input database file name
    out_args = common output argument dictionary from parseCommonArgs
    cluster_func = the function to call for carrying out clustering on distance matrix
    cluster_args = a dictionary of arguments to pass to cluster_func
    
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
        if out_args['failed']:
            fail_handle = getOutputHandle(db_file,
                                          out_label='clone-fail', 
                                          out_dir=out_args['out_dir'], 
                                          out_name=out_args['out_name'], 
                                          out_type=out_type)
            fail_writer = getDbWriter(fail_handle, db_file)
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
                for clone in result.results.values():
                    clone_count += 1
                    for i, rec in enumerate(clone):
                        rec.annotations['CLONE'] = clone_count
                        pass_writer.writerow(rec.toDict())
                        pass_count += 1
                        result.log['CLONE%i-%i' % (clone_count, i + 1)] = str(rec.junction)
    
            else:
                for i, rec in enumerate(result.data):
                    if fail_writer is not None: fail_writer.writerow(rec.toDict())
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
        collect_dict = {'log':log, 'out_files': [pass_handle.name]}
        collect_queue.put(collect_dict)
    except:
        #sys.stderr.write('Exception in collector result processing step\n')
        alive.value = False
        raise

    return None


def collectQueueClust(alive, result_queue, collect_queue, db_file, out_args, cluster_func, cluster_args):
    """
    Assembles results from a queue of individual sequence results and manages log/file I/O

    Arguments: 
    alive = a multiprocessing.Value boolean controlling whether processing continues
            if False exit process
    result_queue = a multiprocessing.Queue holding processQueue results
    collect_queue = a multiprocessing.Queue to store collector return values
    db_file = the input database file name
    out_args = common output argument dictionary from parseCommonArgs
    cluster_func = the function to call for carrying out clustering on distance matrix
    cluster_args = a dictionary of arguments to pass to cluster_func
    
    Returns: 
    None
    (adds 'log' and 'out_files' to collect_dict)
    """
    # Open output files
    try:
               
        # Iterate over Ig records to count and order by junction length
        result_count = 0
        records = {}
        # print 'Reading file...'
        db_iter = readDbFile(db_file)
        for rec in db_iter:
            records[rec.id] = rec
            result_count += 1
        records = OrderedDict(sorted(list(records.items()), key=lambda i: i[1].junction_length))
                
        # Define empty matrix to store assembled results
        dist_mat = np.zeros((result_count,result_count))
        
        # Count records and define output format 
        out_type = getFileType(db_file) if out_args['out_type'] is None \
                   else out_args['out_type']
                   
        # Defined successful output handle
        pass_handle = getOutputHandle(db_file, 
                                      out_label='clone-pass', 
                                      out_dir=out_args['out_dir'], 
                                      out_name=out_args['out_name'], 
                                      out_type=out_type)
        pass_writer = getDbWriter(pass_handle, db_file, add_fields='CLONE')
        
        # Defined failed cloning output handle
        if out_args['failed']:
            fail_handle = getOutputHandle(db_file,
                                          out_label='clone-fail', 
                                          out_dir=out_args['out_dir'], 
                                          out_name=out_args['out_name'], 
                                          out_type=out_type)
            fail_writer = getDbWriter(fail_handle, db_file)
        else:
            fail_handle = None
            fail_writer = None

        # Open log file
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
        row_count = rec_count = 0
        while alive.value:
            # Get result from queue
            if result_queue.empty():  continue
            else:  result = result_queue.get()
            # Exit upon reaching sentinel
            if result is None:  break

            # Print progress for previous iteration
            if row_count == 0:
                print('PROGRESS> Assigning clones')
            printProgress(row_count, result_count, 0.05, start_time)
            
            # Update counts for iteration
            row_count += 1
            rec_count += len(result)
            
            # Add result row to distance matrix
            if result:
                dist_mat[list(range(result_count-len(result),result_count)),result_count-len(result)] = result.results
                
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None    
        
        # Calculate linkage and carry out clustering
        # print dist_mat
        clusters = cluster_func(dist_mat, **cluster_args) if dist_mat is not None else None
        clones = {}
        # print clusters
        for i, c in enumerate(clusters):
            clones.setdefault(c, []).append(records[list(records.keys())[i]])
        
        # Write passed and failed records
        clone_count = pass_count = fail_count = 0
        if clones:
            for clone in clones.values():
                clone_count += 1
                for i, rec in enumerate(clone):
                    rec.annotations['CLONE'] = clone_count
                    pass_writer.writerow(rec.toDict())
                    pass_count += 1
                    #result.log['CLONE%i-%i' % (clone_count, i + 1)] = str(rec.junction)

        else:
            for i, rec in enumerate(result.data):
                fail_writer.writerow(rec.toDict())
                fail_count += 1
                #result.log['CLONE0-%i' % (i + 1)] = str(rec.junction)
        
        # Print final progress
        printProgress(row_count, result_count, 0.05, start_time)
    
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
        alive.value = False
        raise
    
    return None


def defineClones(db_file, feed_func, work_func, collect_func, clone_func, cluster_func=None,
                 group_func=None, group_args={}, clone_args={}, cluster_args={}, 
                 out_args=default_out_args, nproc=None, queue_size=None):
    """
    Define clonally related sequences
    
    Arguments:
    db_file = filename of input database
    feed_func = the function that feeds the queue
    work_func = the worker function that will run on each CPU
    collect_func = the function that collects results from the workers
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
    if group_func is not None:
        log['GROUP_FUNC'] = group_func.__name__
        log['GROUP_ARGS'] = group_args
    log['CLONE_FUNC'] = clone_func.__name__

    # TODO:  this is yucky, but can be fixed by using a model class
    clone_log = clone_args.copy()
    if 'dist_mat' in clone_log:  del clone_log['dist_mat']
    log['CLONE_ARGS'] = clone_log

    if cluster_func is not None:
        log['CLUSTER_FUNC'] = cluster_func.__name__
        log['CLUSTER_ARGS'] = cluster_args
    log['NPROC'] = nproc
    printLog(log)
    
    # Define feeder function and arguments
    feed_args = {'db_file': db_file,
                 'group_func': group_func, 
                 'group_args': group_args}
    # Define worker function and arguments
    work_args = {'clone_func': clone_func, 
                 'clone_args': clone_args}
    # Define collector function and arguments
    collect_args = {'db_file': db_file,
                    'out_args': out_args,
                    'cluster_func': cluster_func,
                    'cluster_args': cluster_args}
    
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

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            formatter_class=CommonHelpFormatter)
    parser.add_argument('--version', action='version',
                        version='%(prog)s:' + ' %s-%s' %(__version__, __date__))
    subparsers = parser.add_subparsers(title='subcommands', dest='command', metavar='',
                                       help='Cloning method')
    # TODO:  This is a temporary fix for Python issue 9253
    subparsers.required = True
    
    # Parent parser    
    parser_parent = getCommonArgParser(seq_in=False, seq_out=False, db_in=True, 
                                       multiproc=True)
    
    # Distance cloning method
    parser_bygroup = subparsers.add_parser('bygroup', parents=[parser_parent],
                                           formatter_class=CommonHelpFormatter,
                                           help='''Defines clones as having same V assignment,
                                                J assignment, and junction length with
                                                specified substitution distance model.''',
                                           description='''Defines clones as having same V assignment,
                                                       J assignment, and junction length with
                                                       specified substitution distance model.''')
    parser_bygroup.add_argument('-f', nargs='+', action='store', dest='fields', default=None,
                             help='Additional fields to use for grouping clones (non VDJ)')
    parser_bygroup.add_argument('--mode', action='store', dest='mode', 
                             choices=('allele', 'gene'), default=default_index_mode,
                             help='''Specifies whether to use the V(D)J allele or gene for
                                  initial grouping.''')
    parser_bygroup.add_argument('--act', action='store', dest='action',
                             choices=('first', 'set'), default=default_index_action,
                             help='''Specifies how to handle multiple V(D)J assignments
                                  for initial grouping.''')
    parser_bygroup.add_argument('--model', action='store', dest='model', 
                             choices=choices_bygroup_model,
                             default=default_bygroup_model,
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
    parser_bygroup.add_argument('--dist', action='store', dest='distance', type=float, 
                             default=default_distance,
                             help='The distance threshold for clonal grouping')
    parser_bygroup.add_argument('--norm', action='store', dest='norm',
                             choices=('len', 'mut', 'none'), default=default_norm,
                             help='''Specifies how to normalize distances. One of none
                                  (do not normalize), len (normalize by length),
                                  or mut (normalize by number of mutations between sequences).''')
    parser_bygroup.add_argument('--sym', action='store', dest='sym',
                             choices=('avg', 'min'), default=default_sym,
                             help='''Specifies how to combine asymmetric distances. One of avg
                                  (average of A->B and B->A) or min (minimum of A->B and B->A).''')
    parser_bygroup.add_argument('--link', action='store', dest='linkage',
                             choices=('single', 'average', 'complete'), default=default_linkage,
                             help='''Type of linkage to use for hierarchical clustering.''')
    parser_bygroup.add_argument('--sf', action='store', dest='seq_field',
                                default=default_seq_field,
                                help='''The name of the field to be used to calculate
                                     distance between records''')
    parser_bygroup.set_defaults(feed_func=feedQueue)
    parser_bygroup.set_defaults(work_func=processQueue)
    parser_bygroup.set_defaults(collect_func=collectQueue)  
    parser_bygroup.set_defaults(group_func=indexJunctions)  
    parser_bygroup.set_defaults(clone_func=distanceClones)
    
    # Chen2010
    parser_chen = subparsers.add_parser('chen2010', parents=[parser_parent],
                                        formatter_class=CommonHelpFormatter,
                                        help='''Defines clones by method specified in Chen, 2010.''',
                                        description='''Defines clones by method specified in Chen, 2010.''')
    parser_chen.set_defaults(feed_func=feedQueueClust)
    parser_chen.set_defaults(work_func=processQueueClust)
    parser_chen.set_defaults(collect_func=collectQueueClust)
    parser_chen.set_defaults(cluster_func=hierClust)

    # Ademokun2011
    parser_ade = subparsers.add_parser('ademokun2011', parents=[parser_parent],
                                        formatter_class=CommonHelpFormatter,
                                        help='''Defines clones by method specified in Ademokun, 2011.''',
                                        description='''Defines clones by method specified in Ademokun, 2011.''')
    parser_ade.set_defaults(feed_func=feedQueueClust)
    parser_ade.set_defaults(work_func=processQueueClust)
    parser_ade.set_defaults(collect_func=collectQueueClust)
    parser_ade.set_defaults(cluster_func=hierClust)
        
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
    
    # Define clone_args
    if args.command == 'bygroup':
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

        del args_dict['fields']
        del args_dict['action']
        del args_dict['mode']
        del args_dict['model']
        del args_dict['distance']
        del args_dict['norm']
        del args_dict['sym']
        del args_dict['linkage']
        del args_dict['seq_field']

    # Define clone_args
    if args.command == 'chen2010':
        args_dict['clone_func'] = distChen2010
        args_dict['cluster_args'] = {'method': args.command }

    if args.command == 'ademokun2011':
        args_dict['clone_func'] = distAdemokun2011
        args_dict['cluster_args'] = {'method': args.command }
    
    # Call defineClones
    del args_dict['command']
    del args_dict['db_files']
    for f in args.__dict__['db_files']:
        args_dict['db_file'] = f
        defineClones(**args_dict)
