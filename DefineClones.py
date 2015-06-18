#!/usr/bin/env python
"""
Assign Ig sequences into clones
"""

__author__    = 'Namita Gupta, Jason Anthony Vander Heiden, Gur Yaari, Mohamed Uduman'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.2.1'
__date__      = '2015.06.18'

# Imports
import os, signal, sys, textwrap, re
import multiprocessing as mp
import numpy as np
from argparse import ArgumentParser
from collections import OrderedDict
from ctypes import c_bool
from itertools import chain
from pandas import DataFrame
from pandas.io.parsers import read_csv
from time import time
from Bio import pairwise2
from Bio.Seq import translate, Seq

# IgCore imports
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from IgCore import default_out_args
from IgCore import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from IgCore import getFileType, getOutputHandle, printLog, printProgress
from IgCore import getScoreDict
from IgCore import manageProcesses
from DbCore import countDbFile, readDbFile, getDbWriter
from DbCore import DbData, DbResult, getDistMat
from DbCore import calcDistances, formClusters

# Defaults
default_translate = False
default_distance = 0.0
default_bygroup_model = 'hs1f'
default_hclust_model = 'chen2010'
default_norm = 'length'
default_linkage = 'single'


# TODO:  moved models into core, can wait until package conversion
# Path to model files
model_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'models')

# Amino acid Hamming distance
aa_model = getDistMat(n_score=1, gap_score=0, alphabet='aa')

# DNA Hamming distance
ham_model = getDistMat(n_score=0, gap_score=0, alphabet='dna')

# Human 1-mer model
hs1f = DataFrame([[0,2.08,1,1.75],
                  [2.08,0,1.75,1],
                  [1,1.75,0,2.08],
                  [1.75,1,2.08,0]],
                 index=['A','C','G','T'],
                 columns=['A','C','G','T'],
                 dtype=float)
hs1f_model = getDistMat(hs1f)

# Mouse 1-mer model
smith96 = DataFrame([[0,2.86,1,2.14],
                     [2.86,0,2.14,1],
                     [1,2.14,0,2.86],
                     [2.14,1,2.86,0]],
                    index=['A','C','G','T'],
                    columns=['A','C','G','T'],
                    dtype=float)
m1n_model = getDistMat(smith96)

# Human 5-mer DNA model
hs5f_file = os.path.join(model_path, 'HS5F_Distance.tab')
hs5f_model = read_csv(hs5f_file, sep='\t', index_col=0)

# TODO:  Merge duplicate feed, process and collect functions.
# TODO:  Update feed, process and collect functions to current pRESTO implementation.


# TODO:  this function is an abstraction to facilitate later cleanup
def getModelMatrix(model):
    """
    Simple wrapper to get distance matrix from model name

    Arguments:
    model = model name

    Return:
    a pandas.DataFrame containing the character distance matrix
    """
    if model == 'aa':
        return(aa_model)
    elif model == 'ham':
        return(ham_model)
    elif model == 'm1n':
        return(m1n_model)
    elif model == 'hs1f':
        return(hs1f_model)
    elif model == 'hs5f':
        return(hs5f_model)
    else:
        sys.stderr.write('Unrecognized distance model: %s.\n' % model)


def indexJunctions(db_iter, fields=None, mode='gene', action='first'):
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

    start_time = time()
    clone_index = {}
    rec_count = 0
    for rec in db_iter:
        key = _get_key(rec, action)

        # Print progress
        if rec_count == 0:
            print 'PROGRESS> Grouping sequences'

        printProgress(rec_count, step=1000, start_time=start_time)
        rec_count += 1

        # Assigned passed preclone records to key and failed to index None
        if all([k is not None for k in key]):
            #print key
            # TODO:  Has much slow. Should have less slow.
            if action == 'set':
                
                f_range = range(2, 3 + (len(fields) if fields else 0) )
                vdj_range = range(2)
                
                # Check for any keys that have matching columns and junction length and overlapping genes/alleles
                to_remove = []
                if len(clone_index) > 0 and not key in clone_index:
                    key = list(key)
                    for k in clone_index:
                        if all([key[i] == k[i] for i in f_range]):
                            if all([not set(key[i]).isdisjoint(set(k[i])) for i in vdj_range]):
                                for i in vdj_range:  key[i] = tuple(set(key[i]).union(set(k[i])))
                                to_remove.append(k)
                
                # Remove original keys, replace with union of all genes/alleles and append values to new key
                val = [rec]
                val += list(chain(*(clone_index.pop(k) for k in to_remove)))
                clone_index[tuple(key)] = clone_index.get(tuple(key),[]) + val 

            elif action == 'first':
                clone_index.setdefault(key, []).append(rec)
        else:
            # TODO: weird return object for missing data case
            clone_index.setdefault((0,0,0), []).append(rec)

    printProgress(rec_count, step=1000, start_time=start_time, end=True)

    return clone_index


def distanceClones(records, model=default_bygroup_model, distance=default_distance,
                   dist_mat=None, norm=default_norm, linkage=default_linkage):
    """
    Separates a set of IgRecords into clones

    Arguments: 
    records = an iterator of IgRecords
    model = substitution model used to calculate distance
    distance = the distance threshold to assign clonal groups
    dist_mat = pandas DataFrame of pairwise nucleotide or amino acid distances
    norm = normalization method
    linkage = type of linkage

    Returns: 
    a dictionary of lists defining {clone number: [IgRecords clonal group]}
    """
    # Get distance matrix if not provided
    if dist_mat is None:  dist_mat = getModelMatrix(model)

    # Determine length of n-mers
    if model in ['hs1f', 'm1n', 'aa', 'ham']:
        nmer_len = 1
    elif model in ['hs5f']:
        nmer_len = 5
    else:
        sys.stderr.write('Unrecognized distance model: %s.\n' % model)

    # Define unique junction mapping
    junc_map = {}
    for ig in records:
        # Check if junction length is 0
        if ig.junction_length == 0:
            return None

        junc = re.sub('[\.-]','N', str(ig.junction))
        if model == 'aa':  junc = translate(junc)

        junc_map.setdefault(junc, []).append(ig)

    # Process records
    if len(junc_map) == 1:
        return {1:records}

    # Define junction sequences
    junctions = list(junc_map.keys())

    # Calculate pairwise distance matrix
    dists = calcDistances(junctions, nmer_len, dist_mat, norm)

    # Perform hierarchical clustering
    clusters = formClusters(dists, linkage, distance)

    # Turn clusters into clone dictionary
    clone_dict = {}
    for i, c in enumerate(clusters):
        clone_dict.setdefault(c, []).extend(junc_map[junctions[i]])

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
    score_dict = getScoreDict()
    
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
    score_dict = getScoreDict()
    
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
        records = OrderedDict(sorted(records.items(), key=lambda i: i[1].junction_length))
        dist_dict = {}
        for __ in range(len(records)):
            k,v = records.popitem(last=False)
            dist_dict[k] = [v].append(records.values())
    except:
        #sys.stderr.write('Exception in feeder grouping step\n')
        alive.value = False
        raise
    
    # Add groups to data queue
    try:
        # print 'START FEED', alive.value
        # Iterate over groups and feed data queue
        dist_iter = dist_dict.iteritems()
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
                print 'PROGRESS> Assigning clones'
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
        records = OrderedDict(sorted(records.items(), key=lambda i: i[1].junction_length))
                
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
                print 'PROGRESS> Assigning clones'
            printProgress(row_count, result_count, 0.05, start_time)
            
            # Update counts for iteration
            row_count += 1
            rec_count += len(result)
            
            # Add result row to distance matrix
            if result:
                dist_mat[range(result_count-len(result),result_count),result_count-len(result)] = result.results
                
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
            clones.setdefault(c, []).append(records[records.keys()[i]])
        
        # Write passed and failed records
        clone_count = pass_count = fail_count = 0
        if clones:
            for clone in clones.itervalues():
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
    fields = textwrap.dedent(
             '''
             output files:
               clone-pass     database with assigned clonal group numbers.
               clone-fail     database with records failing clonal grouping.

             required fields:
               SEQUENCE_ID
               V_CALL or V_CALL_GENOTYPED
               D_CALL
               J_CALL
               JUNCTION_LENGTH
               JUNCTION

             output fields:
               CLONE
              ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,  
                            version='%(prog)s:' + ' v%s-%s' %(__version__, __date__),
                            formatter_class=CommonHelpFormatter)
    subparsers = parser.add_subparsers(title='subcommands', dest='command', metavar='',
                                       help='Cloning method')
    
    # Parent parser    
    parser_parent = getCommonArgParser(seq_in=False, seq_out=False, db_in=True, 
                                       multiproc=True)
    
    # Distance cloning method
    parser_bygroup = subparsers.add_parser('bygroup', parents=[parser_parent],
                                        formatter_class=CommonHelpFormatter,
                                        help='''Defines clones as having same V assignment,
                                              J assignment, and junction length with
                                              specified substitution distance model.''')
    parser_bygroup.add_argument('-f', nargs='+', action='store', dest='fields', default=None,
                             help='Additional fields to use for grouping clones (non VDJ)')
    parser_bygroup.add_argument('--mode', action='store', dest='mode', 
                             choices=('allele', 'gene'), default='gene', 
                             help='''Specifies whether to use the V(D)J allele or gene for
                                  initial grouping.''')
    parser_bygroup.add_argument('--act', action='store', dest='action', default='set',
                             choices=('first', 'set'),
                             help='''Specifies how to handle multiple V(D)J assignments
                                  for initial grouping.''')
    parser_bygroup.add_argument('--model', action='store', dest='model', 
                             choices=('aa', 'ham', 'm1n', 'hs1f', 'hs5f'),
                             default=default_bygroup_model,
                             help='''Specifies which substitution model to use for
                                  calculating distance between junctions. Where m1n is the
                                  mouse single nucleotide transition/trasversion model
                                  of Smith et al, 1996; hs1f is the human single
                                  nucleotide model derived from Yaari et al, 2013; hs5f
                                  is the human S5F model of Yaari et al, 2013; ham is
                                  nucleotide Hamming distance; and aa is amino acid
                                  Hamming distance. The hs5f data should be
                                  considered experimental.''')
    parser_bygroup.add_argument('--dist', action='store', dest='distance', type=float, 
                             default=default_distance,
                             help='The junction distance threshold for clonal grouping')
    parser_bygroup.add_argument('--norm', action='store', dest='norm',
                             choices=('len', 'mut', 'none'), default=default_norm,
                             help='''Specifies how to normalize distances. One of none
                                  (do not normalize), len (normalize by junction length,
                                  or mut (normalize by number of mutations between sequences.''')
    parser_bygroup.add_argument('--link', action='store', dest='linkage',
                             choices=('single', 'average', 'complete'), default=default_linkage,
                             help='''Type of linkage to use for hierarchical clustering.''')
    parser_bygroup.set_defaults(feed_func=feedQueue)
    parser_bygroup.set_defaults(work_func=processQueue)
    parser_bygroup.set_defaults(collect_func=collectQueue)  
    parser_bygroup.set_defaults(group_func=indexJunctions)  
    parser_bygroup.set_defaults(clone_func=distanceClones)
    
    
    # Hierarchical clustering cloning method
    parser_hclust = subparsers.add_parser('hclust', parents=[parser_parent],
                                        formatter_class=CommonHelpFormatter,
                                        help='Defines clones by specified distance metric on CDR3s and \
                                              cutting of hierarchical clustering tree')
#     parser_hclust.add_argument('-f', nargs='+', action='store', dest='fields', default=None,
#                              help='Fields to use for grouping clones (non VDJ)')
    parser_hclust.add_argument('--method', action='store', dest='method', 
                             choices=('chen2010', 'ademokun2011'), default=default_hclust_model, 
                             help='Specifies which cloning method to use for calculating distance \
                                   between CDR3s, computing linkage, and cutting clusters')
    parser_hclust.set_defaults(feed_func=feedQueueClust)
    parser_hclust.set_defaults(work_func=processQueueClust)
    parser_hclust.set_defaults(collect_func=collectQueueClust)
    parser_hclust.set_defaults(cluster_func=hierClust)
        
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
    if args.command == 'bygroup':
        args_dict['group_args'] = {'fields': args_dict['fields'],
                                   'action': args_dict['action'], 
                                   'mode':args_dict['mode']}
        args_dict['clone_args'] = {'model':  args_dict['model'],
                                   'distance':  args_dict['distance'],
                                   'norm': args_dict['norm'],
                                   'linkage': args_dict['linkage']}

        # TODO:  can be cleaned up with abstract model class
        args_dict['clone_args']['dist_mat'] = getModelMatrix(args_dict['model'])

        del args_dict['fields']
        del args_dict['action']
        del args_dict['mode']
        del args_dict['model']
        del args_dict['distance']
        del args_dict['norm']
        del args_dict['linkage']

    # Define clone_args
    if args.command == 'hclust':
        dist_funcs = {'chen2010':distChen2010, 'ademokun2011':distAdemokun2011}
        args_dict['clone_func'] = dist_funcs[args_dict['method']]
        args_dict['cluster_args'] = {'method':  args_dict['method']}
        #del args_dict['fields']
        del args_dict['method']
    
    # Call defineClones
    del args_dict['command']
    del args_dict['db_files']
    for f in args.__dict__['db_files']:
        args_dict['db_file'] = f
        defineClones(**args_dict)