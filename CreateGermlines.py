#!/usr/bin/env python
"""
Reconstructs germline sequences from alignment data
"""
__author__    = 'Namita Gupta, Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.0'
__date__      = '2014.10.2'

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
from DbCore import default_repo, IgRecord

# Defaults
default_germ_types = 'dmask'
default_v_field = 'V_CALL'

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
                repo_dict[tuple(IgRecord.allele_regex.findall(g.description))] = str(g.seq).upper() # @UndefinedVariable
    return repo_dict

    
def joinGermline(align, repo_dict, germ_types, v_field):
    """
    Join gapped germline sequences aligned with sample sequences
    
    Arguments:
    align = iterable yielding dictionaries of sample sequence data
    repo_dict = dictionary of IMGT gapped germline sequences
    germ_types = types of germline sequences to be output
                     (full germline, D-region masked, only V-region germline)
    v_field = field in which to look for V call
    
    Returns:
    dictionary of germline_type: germline_sequence
    """
    j_field = 'J_CALL'
    germs = {'full': '', 'dmask':'', 'vonly':''}
    result_log = OrderedDict()
    result_log['ID'] = align['SEQUENCE_ID']
    
    # Find germline V-Region
    v = align[v_field]
    v_match = IgRecord.allele_regex.findall(v) # @UndefinedVariable
    if v_match:
        vgene = tuple(IgRecord.allele_regex.findall(v)) if v_field == 'V_CALL_GENOTYPED' else (IgRecord.allele_regex.findall(v)[0],) # @UndefinedVariable
        if vgene in repo_dict:
            result_log['V_CALL'] = ','.join(vgene)
            germ_vseq = repo_dict[vgene]
            germ_vseq = germ_vseq[int(align['V_GERM_START'] or 1)-1:int(align['V_GERM_START'] or 1)-1+int(align['V_GERM_LENGTH'] or 0)] + \
                        'N' * (int(align['V_GERM_LENGTH'] or 0) - len(germ_vseq[int(align['V_GERM_START'] or 1)-1:]))
        else:
            result_log['ERROR'] = 'Germline %s not in repertoire' % ','.join(vgene)
            return result_log, germs
    else:
        result_log['V_CALL'] = 'None'
        germ_vseq = 'N' * int(align['V_GERM_LENGTH'] or 0)

    # Find germline D-Region
    if align['D_CALL']:
        dgene = (IgRecord.allele_regex.search(align['D_CALL']).group(0), ) # @UndefinedVariable
        if dgene in repo_dict:
            result_log['D_CALL'] = ','.join(dgene)
            germ_dseq = repo_dict[dgene]
            germ_dseq = germ_dseq[int(align['D_GERM_START'] or 1)-1:int(align['D_GERM_START'] or 1)-1+int(align['D_GERM_LENGTH'] or 0)]
        else:
            result_log['ERROR'] = 'Germline %s not in repertoire' % dgene
            return result_log, germs
    else:
        result_log['D_CALL'] = 'None'
        germ_dseq = ''
    
    # Find germline J-Region
    j = align[j_field]
    j_match = IgRecord.allele_regex.search(j) # @UndefinedVariable
    if j_match is not None:
        jgene = (j_match.group(0), )
        if jgene in repo_dict:
            result_log['J_CALL'] = ','.join(jgene)
            germ_jseq = repo_dict[jgene]
            germ_jseq = germ_jseq[int(align['J_GERM_START'] or 1)-1:int(align['J_GERM_START'] or 1)-1+int(align['J_GERM_LENGTH'] or 0)] + \
                        'N' * (int(align['V_GERM_LENGTH'] or 0) - len(germ_vseq[int(align['V_GERM_START'] or 1)-1:]))
        else:
            result_log['ERROR'] = 'Germline %s not in repertoire' % jgene
            return result_log, germs
    else: 
        result_log['J_CALL'] = 'None'
        germ_jseq = 'N' * int(align['J_GERM_LENGTH'] or 0)
    
    germ_seq = germ_vseq
    regions = 'V' * len(germ_vseq)
    # Nucleotide additions before D (before J for light chains)
    germ_seq += 'N' * (int(align['D_SEQ_START'] or 0) - int(align['V_SEQ_LENGTH'] or 0) - int(align['V_SEQ_START'] or 0))
    regions += 'N' * (int(align['D_SEQ_START'] or 0) - int(align['V_SEQ_LENGTH'] or 0) - int(align['V_SEQ_START'] or 0))
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

    if len(align['SEQUENCE_GAP']) == 0:
        result_log['ERROR'] = 'Gapped sequence is missing from SEQUENCE_GAP column'
    elif len(germs['full']) != len(align['SEQUENCE_GAP']):
        result_log['ERROR'] = 'Germline sequence is %d nucleotides longer than input sequence' % (len(germs['full'])-len(align['SEQUENCE_GAP']))
        
    for v in germs.itervalues(): v = v.upper()
    
    return result_log, germs


def assembleEachGermline(db_file, repo, germ_types, v_field, out_args=default_out_args):
    """
    Write germline sequences to tab-delimited database file
    
    Arguments:
    db_file = input tab-delimited database file
    repo = folder with germline repertoire files
    germ_types = types of germline sequences to be output
                     (full germline, D-region masked, only V-region germline)
    v_field = field in which to look for V call
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
    printLog(log)
    
    # Get repertoire and open Db reader
    repo_dict = getRepo(repo)
    reader = readDbFile(db_file, ig=False)

    # Exit if V call field does not exist in reader
    if v_field not in reader.fieldnames:
        sys.stderr.write('Error: V field does not exist in input database file.')
        exit
    
    # Define log handle
    if out_args['log_file'] is None:  
        log_handle = None
    else:  
        log_handle = open(out_args['log_file'], 'w')
    
    # Create output file handle and Db writer
    pass_handle = getOutputHandle(db_file, 'germ-pass', out_dir=out_args['out_dir'],
                                 out_name=out_args['out_name'], out_type=out_args['out_type'])
    fail_handle = getOutputHandle(db_file, 'germ-fail', out_dir=out_args['out_dir'],
                                 out_name=out_args['out_name'], out_type=out_args['out_type']) if not out_args['clean'] else None
    add_fields = []
    if 'full' in germ_types: add_fields +=  ['GERMLINE_GAP']
    if 'dmask' in germ_types: add_fields += ['GERMLINE_GAP_D_MASK']
    if 'vonly' in germ_types: add_fields += ['GERMLINE_GAP_V_REGION']
    pass_writer = getDbWriter(pass_handle, db_file, add_fields=add_fields)
    fail_writer = getDbWriter(fail_handle, db_file, add_fields=add_fields) if not out_args['clean'] else None
    
    # Initialize time and total count for progress bar
    start_time = time()
    rec_count = countDbFile(db_file)
    pass_count = fail_count = 0
    # Iterate over rows
    for i,row in enumerate(reader):
        # Print progress
        printProgress(i, rec_count, 0.05, start_time)
        
        result_log, germs = joinGermline(row, repo_dict, germ_types, v_field)
        
        # Add germline field(s) to dictionary
        if 'full' in germ_types: row['GERMLINE_GAP'] = germs['full']
        if 'dmask' in germ_types: row['GERMLINE_GAP_D_MASK'] = germs['dmask']
        if 'vonly' in germ_types: row['GERMLINE_GAP_V_REGION'] = germs['vonly']

        # Write row to pass or fail file
        if 'ERROR' in result_log:
            fail_count += 1
            if fail_writer is not None: fail_writer.writerow(row)
        else:
            result_log['SEQUENCE'] = row['SEQUENCE_GAP']
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


def makeCloneGermline(clone, clone_dict, repo_dict, germ_types, v_field, counts, writers, out_args):
    """
    Determine consensus clone sequence and create germline for clone

    Arguments:
    clone = clone ID
    clone_dict = iterable yielding dictionaries of sequence data from clone
    repo_dict = dictionary of IMGT gapped germline sequences
    germ_types = types of germline sequences to be output
                     (full germline, D-region masked, only V-region germline)
    v_field = field in which to look for V call
    counts = dictionary of pass counter and fail counter
    writers = dictionary with pass and fail DB writers
    out_args = arguments for output preferences

    Returns:
    None
    """
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
        if len(val['SEQUENCE_GAP']) > max_length: max_length = len(val['SEQUENCE_GAP'])
    
    # Consensus V and J having most observations
    v_cons = [k for k in v_dict.keys() if v_dict[k] == max(v_dict.values())]
    j_cons = [k for k in j_dict.keys() if j_dict[k] == max(j_dict.values())]
    # Consensus sequence(s) with consensus V/J calls and longest sequence
    cons = [val for val in clone_dict.values() if val.get(v_field,'') in v_cons and \
                                                  val.get(j_field,'') in j_cons and \
                                                  len(val['SEQUENCE_GAP'])==max_length]
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
            cons['J_GERM_LENGTH'] = str(int(cons['J_GERM_LENGTH'] or 0) + max_length - len(cons['SEQUENCE_GAP']))
            cons['SEQUENCE_GAP'] += '.'*(max_length - len(cons['SEQUENCE_GAP']))
            result_log, germs = joinGermline(cons, repo_dict, germ_types, v_field)
            result_log['ID'] = clone
            result_log['CONSENSUS'] = cons['SEQUENCE_ID']
    else:
        cons = cons[0]
        result_log, germs = joinGermline(cons, repo_dict, germ_types, v_field)
        result_log['ID'] = clone
        result_log['CONSENSUS'] = cons['SEQUENCE_ID']

    # Write sequences of clone
    for val in clone_dict.itervalues():
        if 'ERROR' not in result_log:
            # Update lengths padded to longest sequence in clone
            val['J_GERM_LENGTH'] = str(int(val['J_GERM_LENGTH'] or 0) + max_length - len(val['SEQUENCE_GAP']))
            val['SEQUENCE_GAP'] += '.'*(max_length - len(val['SEQUENCE_GAP']))
            
            # Add column(s) to tab-delimited database file
            if 'full' in germ_types: val['GERMLINE_GAP'] = germs['full']
            if 'dmask' in germ_types: val['GERMLINE_GAP_D_MASK'] = germs['dmask']
            if 'vonly' in germ_types: val['GERMLINE_GAP_V_REGION'] = germs['vonly']
            
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
        
        
def assembleCloneGermline(db_file, repo, germ_types, v_field, out_args=default_out_args):
    """
    Assemble one germline sequence for each clone in a tab-delimited database file
    
    Arguments:
    db_file = input tab-delimited database file
    repo = folder with germline repertoire files
    germ_types = types of germline sequences to be output
                     (full germline, D-region masked, only V-region germline)
    v_field = field in which to look for V call
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
    printLog(log)
    
    # Get repertoire and open Db reader
    repo_dict = getRepo(repo)
    reader = readDbFile(db_file, ig=False)

    # Exit if V call field does not exist in reader
    if v_field not in reader.fieldnames:
        sys.stderr.write('Error: V field does not exist in input database file.')
        exit
    
    # Define log handle
    if out_args['log_file'] is None:  
        log_handle = None
    else:  
        log_handle = open(out_args['log_file'], 'w')
    
    # Create output file handle and Db writer
    pass_handle = getOutputHandle(db_file, 'germ-pass', out_dir=out_args['out_dir'],
                                 out_name=out_args['out_name'], out_type=out_args['out_type'])
    fail_handle = getOutputHandle(db_file, 'germ-fail', out_dir=out_args['out_dir'],
                                 out_name=out_args['out_name'], out_type=out_args['out_type']) if not out_args['clean'] else None
    add_fields = []
    if 'full' in germ_types: add_fields +=  ['GERMLINE_GAP']
    if 'dmask' in germ_types: add_fields += ['GERMLINE_GAP_D_MASK']
    if 'vonly' in germ_types: add_fields += ['GERMLINE_GAP_V_REGION']
    writers = {}
    writers['pass'] = getDbWriter(pass_handle, db_file, add_fields=add_fields)
    writers['fail'] = getDbWriter(fail_handle, db_file, add_fields=add_fields) if not out_args['clean'] else None
    
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
            result_log = makeCloneGermline(clone, clone_dict, repo_dict, germ_types, v_field, counts, writers, out_args)
            printLog(result_log, handle=log_handle)
            # Now deal with current row (first of next clone)
            clone = row['CLONE']
            clone_dict = OrderedDict([(row['SEQUENCE_ID'],row)])
        # Last case is only for first row of file
        else:
            clone = row['CLONE']
            clone_dict = OrderedDict([(row['SEQUENCE_ID'],row)])
    clone_count += 1
    result_log = makeCloneGermline(clone, clone_dict, repo_dict, germ_types, v_field, counts, writers, out_args)
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
             required fields:
                 SEQUENCE_ID 
                 SEQUENCE
                 SEQUENCE_GAP
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
                 GERMLINE_GAP
                 GERMLINE_GAP_D_MASK
                 GERMLINE_GAP_V_REGION
              ''')

    # Parent parser
    parser_parent = getCommonArgParser(seq_in=False, seq_out=False, db_in=True, db_out=True, 
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
    parser.add_argument('--cloned', action='store_true', dest='byClone', 
                        help='Specify to create only one germline per clone \
                             (assumes input file is sorted by clone column)')
    parser.add_argument('--vfield', action='store', dest='v_field', default=default_v_field,
                        help='Specify field to use for germline V call')

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
    del args_dict['byClone']
    
    for f in args.__dict__['db_files']:
        args_dict['db_file'] = f
        if args.__dict__['byClone']: assembleCloneGermline(**args_dict)
        else: assembleEachGermline(**args_dict)
