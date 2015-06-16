#!/usr/bin/env python
"""
Reconstructs germline sequences from alignment data
"""
__author__    = 'Namita Gupta, Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.2.0'
__date__      = '2015.05.30'

# Imports
import os, sys, textwrap
from Bio import SeqIO
from argparse import ArgumentParser
from collections import OrderedDict
from time import time

# IgCore and DbCore imports
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from IgCore import default_out_args 
from IgCore import getOutputHandle, printLog, printProgress
from IgCore import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from DbCore import readDbFile, getDbWriter, countDbFile
from DbCore import allele_regex, IgRecord, parseAllele

# Defaults
default_repo = 'germlines'
default_germ_types = 'dmask'
default_v_field = 'V_CALL'
default_seq_field = 'SEQUENCE_IMGT'

def getRepo(repo):
    """
    Parses germline repositories

    Arguments: 
    repo = a string specifying either a directory or file 
           from which to read repository files

    Returns: 
    a dictionary of {allele: sequence} germlines
    """
    if os.path.isdir(repo):
        repo_files = [os.path.join(repo, f) for f in os.listdir(repo)]
    if os.path.isfile(repo):
        with open(repo, 'rU') as file_handle:
            repo_files = [f.strip() for f in file_handle]
    
    repo_dict = {}
    for file_name in repo_files:
        with open(file_name, "rU") as file_handle:
            germlines = SeqIO.parse(file_handle, "fasta")
            for g in germlines:
                germ_key = parseAllele(g.description, allele_regex, 'list')
                repo_dict[germ_key] = str(g.seq).upper() # @UndefinedVariable
    return repo_dict

    
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
    germs = {'full': '', 'dmask':'', 'vonly':''}
    result_log = OrderedDict()
    result_log['ID'] = align['SEQUENCE_ID']

    # Find germline V-Region
    if v_field == 'V_CALL_GENOTYPED':
        vgene = parseAllele(align[v_field], allele_regex, 'list')
        vkey = vgene
    else:
        vgene = parseAllele(align[v_field], allele_regex, 'first')
        vkey = (vgene, )
    if vgene is not None:
        result_log['V_CALL'] = ','.join(vkey)
        if vkey in repo_dict:
            x = int(align['V_GERM_START'] or 1) - 1
            y = x + int(align['V_GERM_LENGTH'] or 0)
            # TODO:  not sure what this line is doing
            z = int(align['V_GERM_LENGTH'] or 0) - len(repo_dict[vkey][int(align['V_GERM_START'] or 1) - 1:])
            germ_vseq = repo_dict[vkey][x:y] + ('N' * z)
        else:
            result_log['ERROR'] = 'Germline %s not in repertoire' % result_log['V_CALL']
            return result_log, germs
    else:
        result_log['V_CALL'] = None
        germ_vseq = 'N' * int(align['V_GERM_LENGTH'] or 0)

    # Find germline D-Region
    dgene = parseAllele(align['D_CALL'], allele_regex, 'first')
    result_log['D_CALL'] = dgene
    if dgene is not None:
        dkey = (dgene, )
        if dkey in repo_dict:
            x = int(align['D_GERM_START'] or 1) - 1
            y = x + int(align['D_GERM_LENGTH'] or 0)
            germ_dseq = repo_dict[dkey][x:y]
        else:
            result_log['ERROR'] = 'Germline %s not in repertoire' % dgene
            return result_log, germs
    else:
        germ_dseq = ''

    # Find germline J-Region
    jgene = parseAllele(align[j_field], allele_regex,'first')
    result_log['J_CALL'] = jgene
    if jgene is not None:
        jkey = (jgene, )
        if jkey in repo_dict:
            x = int(align['J_GERM_START'] or 1) - 1
            y = x + int(align['J_GERM_LENGTH'] or 0)
            # TODO:  not sure what this line is doing either
            z = int(align['V_GERM_LENGTH'] or 0) - len(germ_vseq[int(align['V_GERM_START'] or 1) - 1:])
            germ_jseq = repo_dict[jkey][x:y] + ('N' * z)
        else:
            result_log['ERROR'] = 'Germline %s not in repertoire' % jgene
            return result_log, germs
    else: 
        germ_jseq = 'N' * int(align['J_GERM_LENGTH'] or 0)
    
    germ_seq = germ_vseq
    regions = 'V' * len(germ_vseq)
    # Nucleotide additions before D (before J for light chains)
    # TODO: HACK, the 1 is suppposed to be 'V_SEQ_START' but that isn't working!
    germ_seq += 'N' * (int(align['D_SEQ_START'] or 0) - int(align['V_SEQ_LENGTH'] or 0) - 1)
    regions += 'N' * (int(align['D_SEQ_START'] or 0) - int(align['V_SEQ_LENGTH'] or 0) - 1)
    germ_seq += germ_dseq
    regions += 'D' * len(germ_dseq)
    # Nucleotide additions after D (heavy chains only)
    germ_seq += 'N' * (int(align['J_SEQ_START'] or 0) - int(align['D_SEQ_LENGTH'] or 0) - int(align['D_SEQ_START'] or 0))
    regions += 'N' * (int(align['J_SEQ_START'] or 0) - int(align['D_SEQ_LENGTH'] or 0) - int(align['D_SEQ_START'] or 0))
    germ_seq += germ_jseq
    regions += 'J' * len(germ_jseq)
    germs['full'] = germ_seq.upper()
    germs['regions'] = regions
    if 'dmask' in germ_types: germs['dmask'] = germ_seq[:len(germ_vseq)] + \
                                  "N" * (len(germ_seq) - len(germ_vseq) - len(germ_jseq)) + \
                                  germ_seq[-len(germ_jseq):]
    if 'vonly' in germ_types: germs['vonly'] = germ_vseq

    if len(align[seq_field]) == 0:
        result_log['ERROR'] = 'Gapped sequence is missing from SEQUENCE_GAP column'
    elif len(germs['full']) != len(align[seq_field]):
        result_log['ERROR'] = 'Germline sequence is %d nucleotides longer than input sequence' % (len(germs['full'])-len(align[seq_field]))
        
    for v in germs.itervalues(): v = v.upper()
    
    return result_log, germs


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
    log['GERM_TYPES'] = germ_types if isinstance(germ_types, basestring) else ','.join(germ_types)
    log['CLONED'] = 'False'
    log['V_FIELD'] = v_field
    log['SEQ_FIELD'] = seq_field
    printLog(log)
    
    # Get repertoire and open Db reader
    repo_dict = getRepo(repo)
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
        
        result_log, germs = joinGermline(row, repo_dict, germ_types, v_field, seq_field)
        
        # Add germline field(s) to dictionary
        if 'full' in germ_types: row['GERMLINE_' + seq_type] = germs['full']
        if 'dmask' in germ_types: row['GERMLINE_' + seq_type + '_D_MASK'] = germs['dmask']
        if 'vonly' in germ_types: row['GERMLINE_' + seq_type + '_V_REGION'] = germs['vonly']

        # Write row to pass or fail file
        if 'ERROR' in result_log:
            fail_count += 1
            if fail_writer is not None: fail_writer.writerow(row)
        else:
            result_log['SEQUENCE'] = row[seq_field]
            result_log['GERMLINE'] = germs['full']
            result_log['REGIONS'] = germs['regions']
            
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
    for val in clone_dict.itervalues():
        v = val[v_field]
        v_dict[v] = v_dict.get(v,0) + 1
        j = val[j_field]
        j_dict[j] = j_dict.get(j,0) + 1
        if len(val[seq_field]) > max_length: max_length = len(val[seq_field])
    
    # Consensus V and J having most observations
    v_cons = [k for k in v_dict.keys() if v_dict[k] == max(v_dict.values())]
    j_cons = [k for k in j_dict.keys() if j_dict[k] == max(j_dict.values())]
    # Consensus sequence(s) with consensus V/J calls and longest sequence
    cons = [val for val in clone_dict.values() if val.get(v_field,'') in v_cons and \
                                                  val.get(j_field,'') in j_cons and \
                                                  len(val[seq_field])==max_length]
    # Sequence(s) with consensus V/J are not longest
    if not cons:
        # Sequence(s) with consensus V/J (not longest)
        cons = [val for val in clone_dict.values() if val.get(v_field,'') in v_cons and val.get(j_field,'') in j_cons]
        
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
            result_log, germs = joinGermline(cons, repo_dict, germ_types, v_field, seq_field)
            result_log['ID'] = clone
            result_log['CONSENSUS'] = cons['SEQUENCE_ID']
    else:
        cons = cons[0]
        result_log, germs = joinGermline(cons, repo_dict, germ_types, v_field, seq_field)
        result_log['ID'] = clone
        result_log['CONSENSUS'] = cons['SEQUENCE_ID']

    # Write sequences of clone
    for val in clone_dict.itervalues():
        if 'ERROR' not in result_log:
            # Update lengths padded to longest sequence in clone
            val['J_GERM_LENGTH'] = str(int(val['J_GERM_LENGTH'] or 0) + max_length - len(val[seq_field]))
            val[seq_field] += '.'*(max_length - len(val[seq_field]))
            
            # Add column(s) to tab-delimited database file
            if 'full' in germ_types: val['GERMLINE_' + seq_type] = germs['full']
            if 'dmask' in germ_types: val['GERMLINE_' + seq_type + '_D_MASK'] = germs['dmask']
            if 'vonly' in germ_types: val['GERMLINE_' + seq_type + '_V_REGION'] = germs['vonly']
            
            result_log['SEQUENCE'] = cons
            result_log['GERMLINE'] = germs['full']
            result_log['REGIONS'] = germs['regions']
            
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
    log['GERM_TYPES'] = germ_types if isinstance(germ_types, basestring) else ','.join(germ_types)
    log['CLONED'] = 'True'
    log['V_FIELD'] = v_field
    log['SEQ_FIELD'] = seq_field
    printLog(log)
    
    # Get repertoire and open Db reader
    repo_dict = getRepo(repo)
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
    fields = textwrap.dedent(
             '''
             output files:
                 germ-pass           database with assigned germline sequences.
                 germ-fail           database with records failing germline assignment.

             required fields:
                 SEQUENCE_ID 
                 SEQUENCE_INPUT
                 SEQUENCE_VDJ or SEQUENCE_IMGT
                 V_CALL or V_CALL_GENOTYPED 
                 D_CALL
                 J_CALL
                 V_SEQ_START
                 V_SEQ_LENGTH
                 V_GERM_START
                 V_GERM_LENGTH
                 D_SEQ_START
                 D_SEQ_LENGTH
                 D_GERM_START
                 D_GERM_LENGTH
                 J_SEQ_START
                 J_SEQ_LENGTH
                 J_GERM_START
                 J_GERM_LENGTH
              
              optional fields:
                 CLONE
                
              output fields:
                 GERMLINE_VDJ
                 GERMLINE_VDJ_D_MASK
                 GERMLINE_VDJ_V_REGION
                 GERMLINE_IMGT
                 GERMLINE_IMGT_D_MASK
                 GERMLINE_IMGT_V_REGION
              ''')

    # Parent parser
    parser_parent = getCommonArgParser(seq_in=False, seq_out=False, db_in=True,
                                       annotation=False)
    # Define argument parser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            version='%(prog)s:' + ' v%s-%s' %(__version__, __date__), 
                            parents=[parser_parent],
                            formatter_class=CommonHelpFormatter)
                                     
    parser.add_argument('-r', action='store', dest='repo', default=default_repo,
                        help='Folder where repertoire fasta files are located')
    parser.add_argument('-g', action='store', dest='germ_types', default=default_germ_types,
                        nargs='+', choices=('full','dmask','vonly'),
                        help='Specify type(s) of germlines to include full germline, \
                              germline with D-region masked, or germline for V region only.')
    parser.add_argument('--cloned', action='store_true', dest='cloned',
                        help='Specify to create only one germline per clone \
                             (assumes input file is sorted by clone column)')
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
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)
    del args_dict['db_files']
    del args_dict['cloned']
    
    for f in args.__dict__['db_files']:
        args_dict['db_file'] = f
        if args.__dict__['cloned']:
            assembleCloneGermline(**args_dict)
        else:
            assembleEachGermline(**args_dict)
