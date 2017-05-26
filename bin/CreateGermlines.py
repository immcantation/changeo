#!/usr/bin/env python3
"""
Reconstructs germline sequences from alignment data
"""
# Info
__author__ = 'Namita Gupta, Jason Anthony Vander Heiden'
from changeo import __version__, __date__

# Imports
import os
import sys
from argparse import ArgumentParser
from collections import OrderedDict
from textwrap import dedent
from time import time

# Presto and change imports
from presto.Defaults import default_out_args
from presto.IO import getOutputHandle, printLog, printProgress
from changeo.Commandline import CommonHelpFormatter, checkArgs, getCommonArgParser, parseCommonArgs
from changeo.IO import getDbWriter, readDbFile, countDbFile, readRepo
from changeo.Receptor import allele_regex, parseAllele

# Defaults
default_germ_types = 'dmask'
default_v_field = 'V_CALL'
default_seq_field = 'SEQUENCE_IMGT'

    
def joinGermline(align, repo_dict, germ_types, v_field, seq_field):
    """
    Join gapped germline sequences aligned with sample sequences
    
    Arguments:
    align = iterable yielding dictionaries of sample sequence data
    repo_dict = dictionary of IMGT gapped germline sequences
    germ_types = types of germline sequences to be output
                     (full germline, D-region masked, only V-region germline)
    v_field = field in which to look for V call
    seq_field = field in which to look for sequence
    
    Returns:
    dictionary of germline_type: germline_sequence
    """
    j_field = 'J_CALL'
    germlines = {'full': '', 'dmask': '', 'vonly': '', 'regions': ''}
    result_log = OrderedDict()
    result_log['ID'] = align['SEQUENCE_ID']

    # Find germline V-region gene
    if v_field == 'V_CALL_GENOTYPED':
        vgene = parseAllele(align[v_field], allele_regex, 'list')
        vkey = vgene
    else:
        vgene = parseAllele(align[v_field], allele_regex, 'first')
        vkey = (vgene, )

    try:
        int(align['P3V_LENGTH'])
        int(align['N1_LENGTH'])
        int(align['P5D_LENGTH'])
        int(align['P3D_LENGTH'])
        int(align['N2_LENGTH'])
        int(align['P5J_LENGTH'])
    except:
        regions_style = 'IgBLAST'
    else:
        regions_style = 'IMGT'

    # Build V-region germline
    if vgene is not None:
        result_log['V_CALL'] = ','.join(vkey)
        if vkey in repo_dict:
            vseq = repo_dict[vkey]
            # Germline start
            try: vstart = int(align['V_GERM_START_IMGT']) - 1
            except (TypeError, ValueError): vstart = 0
            # Germline length
            try: vlen = int(align['V_GERM_LENGTH_IMGT'])
            except (TypeError, ValueError): vlen = 0
            # TODO:  not sure what this line is doing here. it no make no sense.
            vpad = vlen - len(vseq[vstart:])
            if vpad < 0: vpad = 0
            germ_vseq = vseq[vstart:(vstart + vlen)] + ('N' * vpad)
        else:
            result_log['ERROR'] = 'Germline %s not in repertoire' % ','.join(vkey)
            return result_log, germlines
    else:
        result_log['V_CALL'] = None
        try: vlen = int(align['V_GERM_LENGTH_IMGT'])
        except (TypeError, ValueError): vlen = 0
        germ_vseq = 'N' * vlen

    # Find germline D-region gene
    dgene = parseAllele(align['D_CALL'], allele_regex, 'first')

    # Build D-region germline
    if dgene is not None:
        result_log['D_CALL'] = dgene
        dkey = (dgene, )
        if dkey in repo_dict:
            dseq = repo_dict[dkey]
            # Germline start
            try: dstart = int(align['D_GERM_START']) - 1
            except (TypeError, ValueError): dstart = 0
            # Germline length
            try: dlen = int(align['D_GERM_LENGTH'])
            except (TypeError, ValueError): dlen = 0
            germ_dseq = repo_dict[dkey][dstart:(dstart + dlen)]
        else:
            result_log['ERROR'] = 'Germline %s not in repertoire' % dgene
            return result_log, germlines
    else:
        result_log['D_CALL'] = None
        germ_dseq = ''

    # Find germline J-region gene
    jgene = parseAllele(align[j_field], allele_regex, 'first')

    # Build D-region germline
    if jgene is not None:
        result_log['J_CALL'] = jgene
        jkey = (jgene, )
        if jkey in repo_dict:
            jseq = repo_dict[jkey]
            # Germline start
            try: jstart = int(align['J_GERM_START']) - 1
            except (TypeError, ValueError): jstart = 0
            # Germline length
            try: jlen = int(align['J_GERM_LENGTH'])
            except (TypeError, ValueError): jlen = 0
            # TODO:  not sure what this line is doing either
            jpad = jlen - len(jseq[jstart:])
            if jpad < 0: jpad = 0
            germ_jseq = jseq[jstart:(jstart + jlen)] + ('N' * jpad)
        else:
            result_log['ERROR'] = 'Germline %s not in repertoire' % jgene
            return result_log, germlines
    else:
        result_log['J_CALL'] = None
        try: jlen = int(align['J_GERM_LENGTH'])
        except (TypeError, ValueError): jlen = 0
        germ_jseq = 'N' * jlen

    # Assemble pieces starting with V-region
    germ_seq = germ_vseq
    regions = 'V' * len(germ_vseq)

    try:
        np1_len = int(align['NP1_LENGTH'])
    except (TypeError, ValueError):
        np1_len = 0

    # NP nucleotide additions after V
    if regions_style == 'IMGT':
        # P nucleotide additions
        try:
            p3v_len = int(align['P3V_LENGTH'])
        except (TypeError, ValueError):
            p3v_len = 0
        if p3v_len < 0:
            result_log['ERROR'] = 'P3V_LENGTH is negative'
            return result_log, germlines

        regions += 'P' * p3v_len

        # N1 nucleotide additions
        try:
            n1_len = int(align['N1_LENGTH'])
        except (TypeError, ValueError):
            n1_len = 0
        if n1_len < 0:
            result_log['ERROR'] = 'N1_LENGTH is negative'
            return result_log, germlines

        regions += 'N' * n1_len

        # P nucleotide additions before D
        try: p5d_len = int(align['P5D_LENGTH'])
        except (TypeError, ValueError): p5d_len = 0
        if p5d_len < 0:
            result_log['ERROR'] = 'P5D_LENGTH is negative'
            return result_log, germlines

        regions += 'P' * p5d_len
    else:
        # IgBLAST style
        # PNP nucleotide additions after V
        if np1_len < 0:
            result_log['ERROR'] = 'NP1_LENGTH is negative'
            return result_log, germlines

        regions += 'N' * np1_len

    germ_seq += 'N' * np1_len

    # Add D-region
    germ_seq += germ_dseq
    regions += 'D' * len(germ_dseq)

    #print 'VD>', germ_seq, '\nVD>', regions

    try:
        np2_len = int(align['NP2_LENGTH'])
    except (TypeError, ValueError):
        np2_len = 0

    # NP nucleotide additions before J
    if regions_style == 'IMGT':
        # P nucleotide additions
        try:
            p3d_len = int(align['P3D_LENGTH'])
        except (TypeError, ValueError):
            p3d_len = 0
        if p3d_len < 0:
            result_log['ERROR'] = 'P3D_LENGTH is negative'
            return result_log, germlines

        regions += 'P' * p3d_len

        # N2 nucleotide additions
        try:
            n2_len = int(align['N2_LENGTH'])
        except (TypeError, ValueError):
            n2_len = 0
        if n2_len < 0:
            result_log['ERROR'] = 'N2_LENGTH is negative'
            return result_log, germlines

        regions += 'N' * n2_len

        # P nucleotide additions
        try:
            p5j_len = int(align['P5J_LENGTH'])
        except (TypeError, ValueError):
            p5j_len = 0
        if p5j_len < 0:
            result_log['ERROR'] = 'P5J_LENGTH is negative'
            return result_log, germlines

        regions += 'P' * p5j_len
    else:
        # IgBLAST style
        # NP nucleotide additions
        if np2_len < 0:
            result_log['ERROR'] = 'NP2_LENGTH is negative'
            return result_log, germlines

        regions += 'N' * np2_len

    germ_seq += 'N' * np2_len

    # Add J-region
    germ_seq += germ_jseq
    regions += 'J' * len(germ_jseq)

    #print('\nREGIONS>',regions,'\n')

    # Define return germlines
    germlines['full'] = germ_seq
    germlines['regions'] = regions

    if 'dmask' in germ_types:
        germlines['dmask'] = germ_seq[:len(germ_vseq)] + \
                             'N' * (len(germ_seq) - len(germ_vseq) - len(germ_jseq)) + \
                             germ_seq[-len(germ_jseq):]
    if 'vonly' in germ_types:
        germlines['vonly'] = germ_vseq

    # Check that input and germline sequence match
    if len(align[seq_field]) == 0:
        result_log['ERROR'] = 'Sequence is missing from %s column' % seq_field
    elif len(germlines['full']) != len(align[seq_field]):
        result_log['ERROR'] = 'Germline sequence is %d nucleotides longer than input sequence' % \
                              (len(germlines['full']) - len(align[seq_field]))

    # Convert to uppercase
    for k, v in germlines.items():  germlines[k] = v.upper()

    return result_log, germlines


def assembleEachGermline(db_file, repo, germ_types, v_field, seq_field, out_args=default_out_args):
    """
    Write germline sequences to tab-delimited database file

    Arguments:
    db_file = input tab-delimited database file
    repo = folder with germline repertoire files
    germ_types = types of germline sequences to be output
                     (full germline, D-region masked, only V-region germline)
    v_field = field in which to look for V call
    seq_field = field in which to look for sequence
    out_args = arguments for output preferences

    Returns:
    None
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'CreateGermlines'
    log['DB_FILE'] = os.path.basename(db_file)
    log['GERM_TYPES'] = germ_types if isinstance(germ_types, str) else ','.join(germ_types)
    log['CLONED'] = 'False'
    log['V_FIELD'] = v_field
    log['SEQ_FIELD'] = seq_field
    printLog(log)

    # Get repertoire and open Db reader
    repo_dict = readRepo(repo)
    reader = readDbFile(db_file, ig=False)

    # Exit if V call field does not exist in reader
    if v_field not in reader.fieldnames:
        sys.exit('Error: V field does not exist in input database file.')

    # Define log handle
    if out_args['log_file'] is None:
        log_handle = None
    else:
        log_handle = open(out_args['log_file'], 'w')

    add_fields = []
    seq_type = seq_field.split('_')[-1]
    if 'full' in germ_types: add_fields +=  ['GERMLINE_' + seq_type]
    if 'dmask' in germ_types: add_fields += ['GERMLINE_' + seq_type + '_D_MASK']
    if 'vonly' in germ_types: add_fields += ['GERMLINE_' + seq_type + '_V_REGION']
    if 'regions' in germ_types: add_fields += ['GERMLINE_REGIONS']

    # Create output file handle and Db writer
    pass_handle = getOutputHandle(db_file, 'germ-pass',
                                  out_dir=out_args['out_dir'],
                                  out_name=out_args['out_name'],
                                  out_type=out_args['out_type'])
    pass_writer = getDbWriter(pass_handle, db_file, add_fields=add_fields)

    if out_args['failed']:
        fail_handle = getOutputHandle(db_file, 'germ-fail',
                                      out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'],
                                      out_type=out_args['out_type'])
        fail_writer = getDbWriter(fail_handle, db_file, add_fields=add_fields)
    else:
        fail_handle = None
        fail_writer = None

    # Initialize time and total count for progress bar
    start_time = time()
    rec_count = countDbFile(db_file)
    pass_count = fail_count = 0
    # Iterate over rows
    for i,row in enumerate(reader):
        # Print progress
        printProgress(i, rec_count, 0.05, start_time)

        result_log, germlines = joinGermline(row, repo_dict, germ_types, v_field, seq_field)

        # Add germline field(s) to dictionary
        if 'full' in germ_types: row['GERMLINE_' + seq_type] = germlines['full']
        if 'dmask' in germ_types: row['GERMLINE_' + seq_type + '_D_MASK'] = germlines['dmask']
        if 'vonly' in germ_types: row['GERMLINE_' + seq_type + '_V_REGION'] = germlines['vonly']
        if 'regions' in germ_types: row['GERMLINE_REGIONS'] = germlines['regions']

        # Write row to pass or fail file
        if 'ERROR' in result_log:
            fail_count += 1
            if fail_writer is not None: fail_writer.writerow(row)
        else:
            result_log['SEQUENCE'] = row[seq_field]
            result_log['GERMLINE'] = germlines['full']
            result_log['REGIONS'] = germlines['regions']

            pass_count += 1
            pass_writer.writerow(row)
        printLog(result_log, handle=log_handle)

    # Print log
    printProgress(i+1, rec_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    log['END'] = 'CreateGermlines'
    printLog(log)

    # Close file handles
    pass_handle.close()
    if fail_handle is not None: fail_handle.close()
    if log_handle is not None:  log_handle.close()


def makeCloneGermline(clone, clone_dict, repo_dict, germ_types, v_field, seq_field, counts, writers, out_args):
    """
    Determine consensus clone sequence and create germline for clone

    Arguments:
    clone = clone ID
    clone_dict = iterable yielding dictionaries of sequence data from clone
    repo_dict = dictionary of IMGT gapped germline sequences
    germ_types = types of germline sequences to be output
                     (full germline, D-region masked, only V-region germline)
    v_field = field in which to look for V call
    seq_field = field in which to look for sequence
    counts = dictionary of pass counter and fail counter
    writers = dictionary with pass and fail DB writers
    out_args = arguments for output preferences

    Returns:
    None
    """
    seq_type = seq_field.split('_')[-1]
    j_field = 'J_CALL'

    # Create dictionaries to count observed V/J calls
    v_dict = OrderedDict()
    j_dict = OrderedDict()

    # Find longest sequence in clone
    max_length = 0
    for val in clone_dict.values():
        v = val[v_field]
        v_dict[v] = v_dict.get(v,0) + 1
        j = val[j_field]
        j_dict[j] = j_dict.get(j,0) + 1
        if len(val[seq_field]) > max_length: max_length = len(val[seq_field])

    # Consensus V and J having most observations
    v_cons = [k for k in list(v_dict.keys()) if v_dict[k] == max(v_dict.values())]
    j_cons = [k for k in list(j_dict.keys()) if j_dict[k] == max(j_dict.values())]
    # Consensus sequence(s) with consensus V/J calls and longest sequence
    cons = [val for val in list(clone_dict.values()) if val.get(v_field,'') in v_cons and \
                                                  val.get(j_field,'') in j_cons and \
                                                  len(val[seq_field])==max_length]
    # Sequence(s) with consensus V/J are not longest
    if not cons:
        # Sequence(s) with consensus V/J (not longest)
        cons = [val for val in list(clone_dict.values()) if val.get(v_field,'') in v_cons and val.get(j_field,'') in j_cons]

        # No sequence has both consensus V and J call
        if not cons:
            result_log = OrderedDict()
            result_log['ID'] = clone
            result_log['V_CALL'] = ','.join(v_cons)
            result_log['J_CALL'] = ','.join(j_cons)
            result_log['ERROR'] = 'No consensus sequence for clone found'
        else:
            # Pad end of consensus sequence with gaps to make it the max length
            cons = cons[0]
            cons['J_GERM_LENGTH'] = str(int(cons['J_GERM_LENGTH'] or 0) + max_length - len(cons[seq_field]))
            cons[seq_field] += '.'*(max_length - len(cons[seq_field]))
            result_log, germlines = joinGermline(cons, repo_dict, germ_types, v_field, seq_field)
            result_log['ID'] = clone
            result_log['CONSENSUS'] = cons['SEQUENCE_ID']
    else:
        cons = cons[0]
        result_log, germlines = joinGermline(cons, repo_dict, germ_types, v_field, seq_field)
        result_log['ID'] = clone
        result_log['CONSENSUS'] = cons['SEQUENCE_ID']

    # Write sequences of clone
    for val in clone_dict.values():
        if 'ERROR' not in result_log:
            # Update lengths padded to longest sequence in clone
            val['J_GERM_LENGTH'] = str(int(val['J_GERM_LENGTH'] or 0) + max_length - len(val[seq_field]))
            val[seq_field] += '.'*(max_length - len(val[seq_field]))

            # Add column(s) to tab-delimited database file
            if 'full' in germ_types: val['GERMLINE_' + seq_type] = germlines['full']
            if 'dmask' in germ_types: val['GERMLINE_' + seq_type + '_D_MASK'] = germlines['dmask']
            if 'vonly' in germ_types: val['GERMLINE_' + seq_type + '_V_REGION'] = germlines['vonly']
            if 'regions' in germ_types: val['GERMLINE_REGIONS'] = germlines['regions']

            # Add field
            val['GERMLINE_V_CALL'] = result_log['V_CALL']
            val['GERMLINE_D_CALL'] = result_log['D_CALL']
            val['GERMLINE_J_CALL'] = result_log['J_CALL']
            
            result_log['SEQUENCE'] = cons[seq_field]
            result_log['GERMLINE'] = germlines['full']
            result_log['REGIONS'] = germlines['regions']

            # Write to pass file
            counts['pass'] += 1
            writers['pass'].writerow(val)
        else:
            # Write to fail file
            counts['fail'] += 1
            if writers['fail'] is not None: writers['fail'].writerow(val)
    # Return log
    return result_log


def assembleCloneGermline(db_file, repo, germ_types, v_field, seq_field, out_args=default_out_args):
    """
    Assemble one germline sequence for each clone in a tab-delimited database file

    Arguments:
    db_file = input tab-delimited database file
    repo = folder with germline repertoire files
    germ_types = types of germline sequences to be output
                     (full germline, D-region masked, only V-region germline)
    v_field = field in which to look for V call
    seq_field = field in which to look for sequence
    out_args = arguments for output preferences

    Returns:
    None
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'CreateGermlines'
    log['DB_FILE'] = os.path.basename(db_file)
    log['GERM_TYPES'] = germ_types if isinstance(germ_types, str) else ','.join(germ_types)
    log['CLONED'] = 'True'
    log['V_FIELD'] = v_field
    log['SEQ_FIELD'] = seq_field
    printLog(log)

    # Get repertoire and open Db reader
    repo_dict = readRepo(repo)
    reader = readDbFile(db_file, ig=False)

    # Exit if V call field does not exist in reader
    if v_field not in reader.fieldnames:
        sys.exit('Error: V field does not exist in input database file.')

    # Define log handle
    if out_args['log_file'] is None:
        log_handle = None
    else:
        log_handle = open(out_args['log_file'], 'w')

    add_fields = []
    seq_type = seq_field.split('_')[-1]
    if 'full' in germ_types: add_fields +=  ['GERMLINE_' + seq_type]
    if 'dmask' in germ_types: add_fields += ['GERMLINE_' + seq_type + '_D_MASK']
    if 'vonly' in germ_types: add_fields += ['GERMLINE_' + seq_type + '_V_REGION']
    if 'regions' in germ_types: add_fields += ['GERMLINE_REGIONS']

    add_fields += ['GERMLINE_V_CALL']
    add_fields += ['GERMLINE_D_CALL']
    add_fields += ['GERMLINE_J_CALL']
    
    # Create output file handle and Db writer
    writers = {}
    pass_handle = getOutputHandle(db_file, 'germ-pass', out_dir=out_args['out_dir'],
                                 out_name=out_args['out_name'], out_type=out_args['out_type'])
    writers['pass'] = getDbWriter(pass_handle, db_file, add_fields=add_fields)

    if out_args['failed']:
        fail_handle = getOutputHandle(db_file, 'germ-fail', out_dir=out_args['out_dir'],
                                     out_name=out_args['out_name'], out_type=out_args['out_type'])
        writers['fail'] = getDbWriter(fail_handle, db_file, add_fields=add_fields)
    else:
        fail_handle = None
        writers['fail'] = None

    # Initialize time and total count for progress bar
    start_time = time()
    rec_count = countDbFile(db_file)
    counts = {}
    clone_count = counts['pass'] = counts['fail'] = 0
    # Iterate over rows
    clone = 'initial'
    clone_dict = OrderedDict()
    for i,row in enumerate(reader):
        # Print progress
        printProgress(i, rec_count, 0.05, start_time)

        # Clone isn't over yet
        if row.get('CLONE','') == clone:
            clone_dict[row["SEQUENCE_ID"]] = row
        # Clone just finished
        elif clone_dict:
            clone_count += 1
            result_log = makeCloneGermline(clone, clone_dict, repo_dict, germ_types,
                                           v_field, seq_field, counts, writers, out_args)
            printLog(result_log, handle=log_handle)
            # Now deal with current row (first of next clone)
            clone = row['CLONE']
            clone_dict = OrderedDict([(row['SEQUENCE_ID'],row)])
        # Last case is only for first row of file
        else:
            clone = row['CLONE']
            clone_dict = OrderedDict([(row['SEQUENCE_ID'],row)])
    clone_count += 1
    result_log = makeCloneGermline(clone, clone_dict, repo_dict, germ_types, v_field,
                                   seq_field, counts, writers, out_args)
    printLog(result_log, handle=log_handle)

    # Print log
    printProgress(i+1, rec_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['CLONES'] = clone_count
    log['RECORDS'] = rec_count
    log['PASS'] = counts['pass']
    log['FAIL'] = counts['fail']
    log['END'] = 'CreateGermlines'
    printLog(log)

    # Close file handles
    pass_handle.close()
    if fail_handle is not None: fail_handle.close()
    if log_handle is not None:  log_handle.close()


def getArgParser():
    """
    Defines the ArgumentParser

    Arguments:
    None

    Returns:
    an ArgumentParser object
    """
    # Define input and output field help message
    fields = dedent(
             '''
             output files:
                 germ-pass
                    database with assigned germline sequences.
                 germ-fail
                    database with records failing germline assignment.

             required fields:
                 SEQUENCE_ID, SEQUENCE_VDJ or SEQUENCE_IMGT,
                 V_CALL or V_CALL_GENOTYPED, D_CALL, J_CALL,
                 V_SEQ_START, V_SEQ_LENGTH, V_GERM_START_IMGT, V_GERM_LENGTH_IMGT,
                 D_SEQ_START, D_SEQ_LENGTH, D_GERM_START, D_GERM_LENGTH,
                 J_SEQ_START, J_SEQ_LENGTH, J_GERM_START, J_GERM_LENGTH,
                 NP1_LENGTH, NP2_LENGTH

             optional fields:
                 N1_LENGTH, N2_LENGTH, P3V_LENGTH, P5D_LENGTH, P3D_LENGTH, P5J_LENGTH,
                 CLONE


             output fields:
                 GERMLINE_VDJ, GERMLINE_VDJ_D_MASK, GERMLINE_VDJ_V_REGION,
                 GERMLINE_IMGT, GERMLINE_IMGT_D_MASK, GERMLINE_IMGT_V_REGION,
                 GERMLINE_V_CALL, GERMLINE_D_CALL, GERMLINE_J_CALL,
                 GERMLINE_REGIONS
              ''')

    # Parent parser
    parser_parent = getCommonArgParser(seq_in=False, seq_out=False, db_in=True,
                                       annotation=False)
    # Define argument parser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            parents=[parser_parent],
                            formatter_class=CommonHelpFormatter)
    parser.add_argument('--version', action='version',
                        version='%(prog)s:' + ' %s-%s' %(__version__, __date__))

    parser.add_argument('-r', nargs='+', action='store', dest='repo', required=True,
                        help='''List of folders and/or fasta files (with .fasta, .fna or .fa
                         extension) with germline sequences.''')
    parser.add_argument('-g', action='store', dest='germ_types', default=default_germ_types,
                        nargs='+', choices=('full', 'dmask', 'vonly', 'regions'),
                        help='''Specify type(s) of germlines to include full germline,
                             germline with D-region masked, or germline for V region only.''')
    parser.add_argument('--cloned', action='store_true', dest='cloned',
                        help='''Specify to create only one germline per clone. Assumes input file is
                             sorted by clone column, and will not yield correct results if the data
                             is unsorted. Note, if allele calls are ambiguous within a clonal group,
                             this will place the germline call used for the entire clone within the
                             GERMLINE_V_CALL, GERMLINE_D_CALL and GERMLINE_J_CALL fields.''')
    parser.add_argument('--vf', action='store', dest='v_field', default=default_v_field,
                        help='Specify field to use for germline V call')
    parser.add_argument('--sf', action='store', dest='seq_field', default=default_seq_field,
                        help='Specify field to use for sequence')

    return parser


if __name__ == "__main__":
    """
    Parses command line arguments and calls main
    """

    # Parse command line arguments
    parser = getArgParser()
    checkArgs(parser)
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)
    del args_dict['db_files']
    del args_dict['cloned']
    args_dict['v_field'] = args_dict['v_field'].upper()
    args_dict['seq_field'] = args_dict['seq_field'].upper()
    
    for f in args.__dict__['db_files']:
        args_dict['db_file'] = f
        if args.__dict__['cloned']:
            assembleCloneGermline(**args_dict)
        else:
            assembleEachGermline(**args_dict)
