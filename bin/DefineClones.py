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
from argparse import ArgumentParser
from collections import OrderedDict
from itertools import chain
from textwrap import dedent
from time import time
from Bio.Seq import translate

# Presto and changeo imports
from presto.Defaults import default_out_args
from presto.IO import printLog, printProgress, printCount, printWarning, printError
from presto.Multiprocessing import manageProcesses
from changeo.Defaults import default_format, default_v_field, default_j_field, default_junction_field
from changeo.Commandline import CommonHelpFormatter, checkArgs, getCommonArgParser, parseCommonArgs
from changeo.Distance import distance_models, calcDistances, formClusters
from changeo.IO import countDbFile, getDbFields, getFormatOperators, getOutputHandle, \
                       AIRRWriter, ChangeoWriter
from changeo.Multiprocessing import DbResult, feedDbQueue, processDbQueue

# Defaults
default_translate = False
default_distance = 0.0
default_index_mode = 'gene'
default_index_action = 'set'
default_distance_model = 'ham'
default_norm = 'len'
default_sym = 'avg'
default_linkage = 'single'
default_max_missing=0
choices_distance_model = ('ham', 'aa', 'hh_s1f', 'hh_s5f',
                          'mk_rs1nf', 'mk_rs5nf',
                          'hs1f_compat', 'm1n_compat')


def filterMissing(data, seq_field=default_junction_field, v_field=default_v_field,
                  j_field=default_j_field, max_missing=default_max_missing):
    """
    Splits a set of sequence into passed and failed groups based on the number
    of missing characters in the sequence

    Arguments:
        data : changeo.Multiprocessing.DbData object.
        seq_field : sequence field to filter on.
        v_field : field containing the V call.
        j_field : field containing the J call.
        max_missing : maximum number of missing characters (non-ACGT) to permit before failing the record.

    Returns:
        changeo.Multiprocessing.DbResult : objected containing filtered records.
    """
    # Function to validate the sequence string
    def _pass(seq):
        if len(seq) > 0 and len(re.findall(r'[^ACGT]', seq)) <= max_missing:
            return True
        else:
            return False

    # Define result object for iteration and get data records
    result = DbResult(data.id, data.data)

    if not data:
        result.data_pass = []
        result.data_fail = data.data
        return result

    result.data_pass = []
    result.data_fail = []
    for rec in data.data:
        seq = rec.getField(seq_field)
        if _pass(seq):  result.data_pass.append(rec)
        else:  result.data_fail.append(rec)

    # Add V(D)J to log
    result.log['ID'] = ','.join([str(x) for x in data.id])
    result.log['VCALL'] = ','.join(set([(r.getVAllele(field=v_field) or '') for r in data.data]))
    result.log['JCALL'] = ','.join(set([(r.getJAllele(field=j_field) or '') for r in data.data]))
    result.log['JUNCLEN'] = ','.join(set([(str(len(r.junction)) or '0') for r in data.data]))
    result.log['CLONED'] = len(result.data_pass)
    result.log['FILTERED'] = len(result.data_fail)

    return result


def indexByIdentity(index, key, rec, group_fields=None):
    """
    Updates a preclone index with a simple key

    Arguments:
      index : preclone index from groupByGene
      key : index key
      rec : Receptor to add to the index
      group_fields : additional annotation fields to use to group preclones;
                     if None use only V, J and junction length

    Returns:
      None : Updates index with new key and records.
    """
    index.setdefault(tuple(key), []).append(rec)


def indexByUnion(index, key, rec, group_fields=None):
    """
    Updates a preclone index with the union of nested keys

    Arguments:
      index : preclone index from groupByGene
      key : index key
      rec : Receptor to add to the index
      group_fields : additional annotation fields to use to group preclones;
                     if None use only V, J and junction length

    Returns:
      None : Updates index with new key and records.
    """
    # List of values for this/new key
    val = [rec]
    f_range = list(range(2, 3 + (len(group_fields) if group_fields else 0)))

    # See if field/junction length combination exists in index
    outer_dict = index
    for field in f_range:
        try:
            outer_dict = outer_dict[key[field]]
        except KeyError:
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


def groupByGene(db_iter, group_fields=None, v_field=default_v_field, j_field=default_j_field,
                mode=default_index_mode, action=default_index_action):
    """
    Identifies preclonal groups by V, J and junction length

    Arguments: 
      db_iter : an iterator of Receptor objects defined by ChangeoReader
      group_fields : additional annotation fields to use to group preclones;
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
    if mode == 'allele' and group_fields is None:
        def _get_key(rec, act):
            return [rec.getVAllele(act, field=v_field), rec.getJAllele(act, field=j_field),
                    None if rec.junction is None else len(rec.junction)]
    elif mode == 'gene' and group_fields is None:
        def _get_key(rec, act):  
            return [rec.getVGene(act, field=v_field), rec.getJGene(act, field=j_field),
                    None if rec.junction is None else len(rec.junction)]
    elif mode == 'allele' and group_fields is not None:
        def _get_key(rec, act):
            vdj = [rec.getVAllele(act, field=v_field), rec.getJAllele(act, field=j_field),
                    None if rec.junction is None else len(rec.junction)]
            ann = [rec.getField(k) for k in group_fields]
            return list(chain(vdj, ann))
    elif mode == 'gene' and group_fields is not None:
        def _get_key(rec, act):
            vdj = [rec.getVGene(act, field=v_field), rec.getJGene(act, field=j_field),
                    None if rec.junction is None else len(rec.junction)]
            ann = [rec.getField(k) for k in group_fields]
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
        printCount(rec_count, step=1000, start_time=start_time, task='Grouping sequences')
        rec_count += 1

        # Assigned passed preclone records to key and failed to index None
        if all([k is not None and k != '' for k in key]):
            # Update index dictionary
            index_func(clone_index, key, rec, group_fields)
        else:
            clone_index.setdefault(None, []).append(rec)

    printCount(rec_count, step=1000, start_time=start_time, task='Grouping sequences', end=True)

    if action == 'set':
        clone_index = _flatten_dict(clone_index)

    return clone_index


def distanceClones(result, seq_field=default_junction_field, model=default_distance_model,
                   distance=default_distance, dist_mat=None, norm=default_norm, sym=default_sym,
                   linkage=default_linkage):
    """
    Separates a set of Receptor objects into clones

    Arguments: 
      result : a changeo.Multiprocessing.DbResult object with filtered records to clone
      seq_field : sequence field used to calculate distance between records
      model : substitution model used to calculate distance
      distance : the distance threshold to assign clonal groups
      dist_mat : pandas DataFrame of pairwise nucleotide or amino acid distances
      norm : normalization method
      sym : symmetry method
      linkage : type of linkage

    Returns: 
      changeo.Multiprocessing.DbResult : an updated DbResult object
    """
    # Get distance matrix if not provided
    if dist_mat is None:
        try:
            dist_mat = distance_models[model]
        except KeyError:
            printError('Unrecognized distance model: %s' % args_dict['model'])

    # TODO:  can be cleaned up with abstract model class
    # Determine length of n-mers
    if model in ['hs1f_compat', 'm1n_compat', 'aa', 'ham', 'hh_s1f', 'mk_rs1nf']:
        nmer_len = 1
    elif model in ['hh_s5f', 'mk_rs5nf']:
        nmer_len = 5
    else:
        printError('Unrecognized distance model: %s.\n' % model)

    # Define unique junction mapping
    seq_map = {}
    for rec in result.data_pass:
        seq = rec.getField(seq_field)
        seq = re.sub('[\.-]', 'N', seq)
        if model == 'aa':  seq = translate(seq)
        seq_map.setdefault(seq, []).append(rec)

    # Define sequences
    sequences = list(seq_map.keys())

    # Zero record case
    if not sequences:
        result.valid = False
        result.log['CLONES'] = 0
        return result

    # Single record case
    if len(sequences) == 1:
        result.results = {1: result.data_pass}
        result.valid = True
        result.log['CLONES'] = 1
        return result

    # Calculate pairwise distance matrix
    dists = calcDistances(sequences, nmer_len, dist_mat, sym=sym, norm=norm)

    # Perform hierarchical clustering
    clusters = formClusters(dists, linkage, distance)

    # Turn clusters into clone dictionary
    clone_dict = {}
    for i, c in enumerate(clusters):
        clone_dict.setdefault(c, []).extend(seq_map[sequences[i]])

    if clone_dict:
        result.results = clone_dict
        result.valid = True
        result.log['CLONES'] = len(clone_dict)
    else:
        result.log['CLONES'] = 0

    return result


def collectQueue(alive, result_queue, collect_queue, db_file, fields,
                 writer=ChangeoWriter, out_file=None, out_args=default_out_args):
    """
    Assembles results from a queue of individual sequence results and manages log/file I/O

    Arguments: 
      alive = a multiprocessing.Value boolean controlling whether processing continues
              if False exit process
      result_queue : a multiprocessing.Queue holding processQueue results
      collect_queue : a multiprocessing.Queue to store collector return values
      db_file : the input database file name
      fields : list of output field names
      writer : writer class.
      out_file : output file name. Automatically generated from the input file if None.
      out_args : common output argument dictionary from parseCommonArgs
    
    Returns:
       None : Adds a dictionary with key value pairs to collect_queue containing
            'log' defining a log object along with the 'pass' and 'fail' output file names.
    """
    # Wrapper for opening handles and writers
    def _open(x, f, writer=writer, out_file=out_file):
        if out_file is not None and x == 'pass':
            handle = open(out_file, 'w')
        else:
            handle = getOutputHandle(db_file,
                                     out_label='clone-%s' % x,
                                     out_dir=out_args['out_dir'],
                                     out_name=out_args['out_name'],
                                     out_type=out_args['out_type'])
        return handle, writer(handle, fields=f)

    # Open log file
    try:
        # Count input records
        result_count = countDbFile(db_file)

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
        # Initialize handles, writers and counters
        pass_handle, pass_writer = None, None
        fail_handle, fail_writer = None, None
        rec_count, clone_count, pass_count, fail_count = 0, 0, 0, 0
        start_time = time()

        # Iterator over results queue until sentinel object reached
        while alive.value:
            # Get result from queue
            if result_queue.empty():  continue
            else:  result = result_queue.get()
            # Exit upon reaching sentinel
            if result is None:  break

            # Print progress for previous iteration and update record count
            printProgress(rec_count, result_count, 0.05, start_time=start_time, task='Assigning clones')
            rec_count += len(result.data)
            
            # Write passed and failed records
            if result:
                # Writing passing sequences
                for clone in result.results.values():
                    clone_count += 1
                    for i, rec in enumerate(clone, start=1):
                        pass_count += 1
                        rec.setField('clone', str(clone_count))
                        result.log['CLONE%i-%i' % (clone_count, i)] = rec.junction
                        try:
                            pass_writer.writeReceptor(rec)
                        except AttributeError:
                            # Open pass file and define writer object
                            pass_handle, pass_writer = _open('pass', fields)
                            pass_writer.writeReceptor(rec)

                # Write failed sequences from passing sets
                if result.data_fail:
                    # Write failed sequences
                    for i, rec in enumerate(result.data_fail, start=1):
                        fail_count += 1
                        result.log['FAIL%i-%i' % (clone_count, i)] = rec.junction
                        if out_args['failed']:
                            try:
                                fail_writer.writeReceptor(rec)
                            except AttributeError:
                                # Open fail file and define writer object
                                fail_handle, fail_writer = _open('fail', fields)
                                fail_writer.writeReceptor(rec)
            else:
                # Write failing records
                for i, rec in enumerate(result.data, start=1):
                    fail_count += 1
                    result.log['CLONE0-%i' % (i)] = rec.junction
                    if out_args['failed']:
                        try:
                            fail_writer.writeReceptor(rec)
                        except AttributeError:
                            # Open fail file and define writer object
                            fail_handle, fail_writer = _open('fail', fields)
                            fail_writer.writeReceptor(rec)
                    
            # Write log
            printLog(result.log, handle=log_handle)
        else:
            sys.stderr.write('PID %s>  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
        
        # Print total counts
        printProgress(rec_count, result_count, 0.05, start_time=start_time, task='Assigning clones')

        # Update return list
        log = OrderedDict()
        log['OUTPUT'] = os.path.basename(pass_handle.name) if pass_handle is not None else None
        log['CLONES'] = clone_count
        log['RECORDS'] = rec_count
        log['PASS'] = pass_count
        log['FAIL'] = fail_count

        # Close file handles and generate return data
        collect_dict = {'log': log, 'pass': None, 'fail': None}
        if pass_handle is not None:
            collect_dict['pass'] = pass_handle.name
            pass_handle.close()
        if fail_handle is not None:
            collect_dict['fail'] = fail_handle.name
            fail_handle.close()
        if log_handle is not None:
            log_handle.close()
        collect_queue.put(collect_dict)
    except:
        alive.value = False
        raise

    return None


def defineClones(db_file, seq_field=default_junction_field, v_field=default_v_field,
                 j_field=default_j_field, max_missing=default_max_missing,
                 group_fields=None, group_func=groupByGene, group_args={},
                 clone_func=distanceClones, clone_args={},
                 format=default_format, out_file=None, out_args=default_out_args,
                 nproc=None, queue_size=None):
    """
    Define clonally related sequences
    
    Arguments:
      db_file : filename of input database.
      seq_field : sequence field used to determine clones.
      v_field : field containing the V call.
      j_field : field containing the J call.
      max_missing : maximum number of non-ACGT characters to allow in the junction sequence.
      group_fields : additional annotation fields to use to group preclones;
                     if None use only V and J.
      group_func : the function to use for assigning preclones.
      group_args : a dictionary of arguments to pass to group_func.
      clone_func : the function to use for determining clones within preclonal groups.
      clone_args : a dictionary of arguments to pass to clone_func.
      format : input and output format.
      out_file : output file name. Automatically generated from the input file if None.
      out_args : common output argument dictionary from parseCommonArgs.
      nproc : the number of processQueue processes;
              if None defaults to the number of CPUs.
      queue_size : maximum size of the argument queue;
                   if None defaults to 2*nproc.
    
    Returns:
      dict: dictionary of output pass and fail files.
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'DefineClones'
    log['FILE'] = os.path.basename(db_file)
    log['SEQ_FIELD'] = seq_field
    log['V_FIELD'] = v_field
    log['J_FIELD'] = j_field
    log['MAX_MISSING'] = max_missing
    log['GROUP_FIELDS'] = ','.join(group_fields) if group_fields is not None else None
    for k in sorted(group_args):
        log[k.upper()] = group_args[k]
    for k in sorted(clone_args):
        if k != 'dist_mat':  log[k.upper()] = clone_args[k]
    log['NPROC'] = nproc
    printLog(log)

    # Define format operators
    try:
        reader, writer, schema = getFormatOperators(format)
    except ValueError:
        printError('Invalid format %s.' % format)

    # Translate to Receptor attribute names
    seq_field = schema.toReceptor(seq_field)
    v_field = schema.toReceptor(v_field)
    j_field = schema.toReceptor(j_field)
    if group_fields is not None:
        group_fields = [schema.toReceptor(f) for f in group_fields]

    # Define feeder function and arguments
    group_args['group_fields'] = group_fields
    group_args['v_field'] = v_field
    group_args['j_field'] = j_field
    feed_args = {'db_file': db_file,
                 'reader': reader,
                 'group_func': group_func, 
                 'group_args': group_args}

    # Define worker function and arguments
    filter_args = {'seq_field': seq_field,
                   'v_field': v_field,
                   'j_field': j_field,
                   'max_missing': max_missing}
    clone_args['seq_field'] = seq_field
    work_args = {'process_func': clone_func,
                 'process_args': clone_args,
                 'filter_func': filterMissing,
                 'filter_args': filter_args}

    # Define collector function and arguments
    out_fields = getDbFields(db_file, add=schema.fromReceptor('clone'), reader=reader)
    out_args['out_type'] = schema.out_type
    collect_args = {'db_file': db_file,
                    'fields': out_fields,
                    'writer': writer,
                    'out_file': out_file,
                    'out_args': out_args}

    # Call process manager
    result = manageProcesses(feed_func=feedDbQueue, work_func=processDbQueue, collect_func=collectQueue,
                             feed_args=feed_args, work_args=work_args, collect_args=collect_args,
                             nproc=nproc, queue_size=queue_size)
        
    # Print log
    result['log']['END'] = 'DefineClones'
    printLog(result['log'])
    output = {k: v for k, v in result.items() if k in ('pass', 'fail')}

    return output


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
                 SEQUENCE_ID, V_CALL, J_CALL, JUNCTION
                
             output fields:
                 CLONE
              ''')
    # Define argument parser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            parents=[getCommonArgParser(format=False, multiproc=True)],
                            formatter_class=CommonHelpFormatter, add_help=False)

    # Distance cloning method
    group = parser.add_argument_group('cloning arguments')
    group.add_argument('--sf', action='store', dest='seq_field', default=default_junction_field,
                        help='Field to be used to calculate distance between records.')
    group.add_argument('--vf', action='store', dest='v_field', default=default_v_field,
                        help='Field containing the germline V segment call.')
    group.add_argument('--jf', action='store', dest='j_field', default=default_j_field,
                        help='Field containing the germline J segment call.')
    group.add_argument('--gf', nargs='+', action='store', dest='group_fields', default=None,
                        help='Additional fields to use for grouping clones aside from V, J and junction length.')
    group.add_argument('--mode', action='store', dest='mode',
                        choices=('allele', 'gene'), default=default_index_mode,
                        help='''Specifies whether to use the V(D)J allele or gene for
                             initial grouping.''')
    group.add_argument('--act', action='store', dest='action',
                        choices=('first', 'set'), default=default_index_action,
                        help='''Specifies how to handle multiple V(D)J assignments for initial grouping. 
                             The "first" action will use only the first gene listed.
                             The "set" action will use all gene assignments and construct a larger gene
                             grouping composed of any sequences sharing an assignment or linked to another
                             sequence by a common assignment (similar to single-linkage).''')
    group.add_argument('--model', action='store', dest='model',
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
    group.add_argument('--dist', action='store', dest='distance', type=float,
                        default=default_distance,
                        help='The distance threshold for clonal grouping')
    group.add_argument('--norm', action='store', dest='norm',
                        choices=('len', 'mut', 'none'), default=default_norm,
                        help='''Specifies how to normalize distances. One of none
                             (do not normalize), len (normalize by length),
                             or mut (normalize by number of mutations between sequences).''')
    group.add_argument('--sym', action='store', dest='sym',
                        choices=('avg', 'min'), default=default_sym,
                        help='''Specifies how to combine asymmetric distances. One of avg
                             (average of A->B and B->A) or min (minimum of A->B and B->A).''')
    group.add_argument('--link', action='store', dest='linkage',
                        choices=('single', 'average', 'complete'), default=default_linkage,
                        help='''Type of linkage to use for hierarchical clustering.''')
    group.add_argument('--maxmiss', action='store', dest='max_missing', type=int,
                        default=default_max_missing,
                        help='''The maximum number of non-ACGT characters (gaps or Ns) to 
                             permit in the junction sequence before excluding the record 
                             from clonal assignment. Note, under single linkage 
                             non-informative positions can create artifactual links 
                             between unrelated sequences. Use with caution.''')
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

    # # Set default fields if not specified.
    # default_fields = {'seq_field': default_junction_field,
    #                   'v_field': default_v_field,
    #                   'j_field': default_j_field}
    #
    # # Default Change-O fields
    # if args_dict['format'] == 'changeo':
    #     for f in default_fields:
    #         if args_dict[f] is None:  args_dict[f] = default_fields[f]
    #         else: args_dict[f] = args_dict[f].upper()
    #
    # # Default AIRR fields
    # if args_dict['format'] == 'airr':
    #     for f in default_fields:
    #         if args_dict[f] is None:  args_dict[f] = ChangeoSchema.toAIRR(default_fields[f])
    #         else: args_dict[f] = args_dict[f].lower()

    # Define grouping and cloning function arguments
    args_dict['group_args'] = {'action': args_dict['action'],
                               'mode':args_dict['mode']}
    args_dict['clone_args'] = {'model':  args_dict['model'],
                               'distance':  args_dict['distance'],
                               'norm': args_dict['norm'],
                               'sym': args_dict['sym'],
                               'linkage': args_dict['linkage']}

    # Get distance matrix
    try:
        args_dict['clone_args']['dist_mat'] = distance_models[args_dict['model']]
    except KeyError:
        printError('Unrecognized distance model: %s' % args_dict['model'])

    # Clean argument dictionary
    del args_dict['action']
    del args_dict['mode']
    del args_dict['model']
    del args_dict['distance']
    del args_dict['norm']
    del args_dict['sym']
    del args_dict['linkage']

    # Clean arguments dictionary
    del args_dict['db_files']
    if 'out_files' in args_dict: del args_dict['out_files']

    # Call main function for each input file
    for i, f in enumerate(args.__dict__['db_files']):
        args_dict['db_file'] = f
        args_dict['out_file'] = args.__dict__['out_files'][i] \
            if args.__dict__['out_files'] else None
        defineClones(**args_dict)
