#!/usr/bin/env python
"""
Create tab-delimited database file to store sequence alignment information
"""

__author__    = 'Namita Gupta, Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.2.1'
__date__      = '2015.06.18'

# Imports
import csv, os, re, sys, textwrap
import pandas as pd
from zipfile import ZipFile, is_zipfile
from tempfile import mkdtemp
from shutil import rmtree
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from argparse import ArgumentParser
from itertools import izip, groupby
from collections import OrderedDict
from time import time


# IgCore and DbCore imports 
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from IgCore import default_out_args, countSeqFile, readSeqFile
from IgCore import parseAnnotation, printLog, printProgress
from IgCore import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from DbCore import getDbWriter, countDbFile
from DbCore import IgRecord, parseAllele
from DbCore import v_allele_regex, d_allele_regex, j_allele_regex

# Default parameters
default_delimiter = ('\t', ',', '-')


def getSeqforIgBlast(seq_file):
    """
    Fetch input sequences for IgBlast queries

    Arguments:
    seq_file = a fasta file of sequences input to IgBlast

    Returns:
    a dictionary of {ID:Seq}
    """

    seq_dict = SeqIO.index(seq_file, "fasta", IUPAC.ambiguous_dna)

    # Create a seq_dict ID translation using IDs truncate up to space or 50 chars
    seqs = {}
    for seq in seq_dict.itervalues():
        seqs.update({seq.description:str(seq.seq)})

    return seqs


def findLine(handle, query):
    """
    Finds line with query string in file

    Arguments:
    handle = file handle in which to search for line
    query = query string for which to search in file

    Returns:
    line from handle in which query string was found
    """
    for line in handle:
        if(re.match(query, line)):
            return line


def extractIMGT(imgt_output):
    """
    Extract necessary files from IMGT results, zipped or unzipped
    
    Arguments:
    imgt_output = zipped file or unzipped folder output by IMGT
    
    Returns:
    sorted list of filenames from which information will be read
    """
    #file_ext = os.path.splitext(imgt_output)[1].lower()
    imgt_flags = ('1_Summary', '2_IMGT-gapped', '3_Nt-sequences', '6_Junction')
    temp_dir = mkdtemp()
    if is_zipfile(imgt_output):
        # Open zip file
        # TODO: upon switch to Python 3 we should be able to use compression=zipfile.ZIP_LZMA for .txz files
        imgt_zip = ZipFile(imgt_output, 'r')
        # Extract required files
        imgt_files = sorted([n for n in imgt_zip.namelist() \
                             if os.path.basename(n).startswith(imgt_flags)])
        imgt_zip.extractall(temp_dir, imgt_files)
        # Define file list
        imgt_files = [os.path.join(temp_dir, f) for f in imgt_files]
    elif os.path.isdir(imgt_output):
        # Find required files in folder
        folder_files = []
        for root, dirs, files in os.walk(imgt_output):
            folder_files.extend([os.path.join(os.path.abspath(root), f) for f in files])
        # Define file list
        imgt_files = sorted([n for n in folder_files \
                             if os.path.basename(n).startswith(imgt_flags)])
    else:
        sys.exit('ERROR: Unsupported IGMT output file. Must be either a zipped file (.zip) or a folder.')
    
    if len(imgt_files) > len(imgt_flags): # e.g. multiple 1_Summary files
        sys.exit('ERROR: Wrong files in IMGT output %s.' % imgt_output)
    elif len(imgt_files) < len(imgt_flags):
        sys.exit('ERROR: Missing necessary file IMGT output %s.' % imgt_output)
        
    return temp_dir, imgt_files


# TODO: return a dictionary with keys determined by the comment strings in the blocks, thus avoiding problems with missing blocks
def readOneIgBlastResult(block):
    """
    Parse a single IgBLAST query result

    Arguments:
    block =  itertools groupby object of single result

    Returns:
    None if no results, otherwise list of DataFrames for each result block
    """
    results = list()
    i = 0
    for match, subblock in groupby(block, lambda l: l=='\n'):
        if not match:
            # Strip whitespace and comments
            sub = [s.strip() for s in subblock if not s.startswith('#')]

            # Continue on empty block
            if not sub:  continue
            else:  i += 1

            # Split by tabs
            sub = [s.split('\t') for s in sub]

            # Append list for "V-(D)-J rearrangement summary" (i == 1)
            # And "V-(D)-J junction details" (i == 2)
            # Otherwise append DataFrame of subblock
            if i == 1 or i == 2:
                results.append(sub[0])
            else:
                df = pd.DataFrame(sub)
                if not df.empty: results.append(df)

    return results if results else None


# TODO:  needs more speeds. pandas is probably to blame.
def readIgBlast(igblast_output, seq_dict):
    """
    Reads IgBlast output

    Arguments:
    igblast_output = IgBlast output file (format 7)
    seq_dict = a dictionary of {ID:Seq} from input fasta file

    Returns:
    a generator of dictionaries containing alignment data
    """

    # Open IgBlast output file
    with open(igblast_output) as f:
        # Iterate over individual results (separated by # IGBLASTN)
        for k1, block in groupby(f, lambda x: re.match('# IGBLASTN', x)):
            if not k1:
                # Extract sequence ID
                query_name = ' '.join(block.next().strip().split(' ')[2:])
                # Initialize db_gen to have ID and input sequence
                db_gen = {'SEQUENCE_ID':     query_name,
                          'SEQUENCE_INPUT':  seq_dict[query_name]}

                # Parse further sub-blocks
                block_list = readOneIgBlastResult(block)

                # If results exist, parse further to obtain full db_gen
                if block_list is not None:
                    # Parse quality information
                    db_gen['STOP'] = 'T' if block_list[0][-4] == 'Yes' else 'F'
                    db_gen['IN_FRAME'] = 'T' if block_list[0][-3] == 'In-frame' else 'F'
                    db_gen['FUNCTIONAL'] = 'T' if block_list[0][-2] == 'Yes' else 'F'

                    # Parse V, D, and J calls
                    call_str = ' '.join(block_list[0])
                    v_call = parseAllele(call_str, v_allele_regex, action='list')
                    d_call = parseAllele(call_str, d_allele_regex, action='list')
                    j_call = parseAllele(call_str, j_allele_regex, action='list')
                    db_gen['V_CALL'] = ','.join(v_call) if v_call is not None else 'None'
                    db_gen['D_CALL'] = ','.join(d_call) if d_call is not None else 'None'
                    db_gen['J_CALL'] = ','.join(j_call) if j_call is not None else 'None'

                    # Parse junction sequence
                    db_gen['JUNCTION'] = re.sub('(N/A)|\[|\(|\)|\]', '', ''.join(block_list[1]))
                    db_gen['JUNCTION_LENGTH'] = len(db_gen['JUNCTION'])

                    # TODO:  IgBLAST does a stupid and doesn't output block #3 sometimes. why?
                    # TODO:  maybe we should fail these. they look craptastic.
                    #pd.set_option('display.width', 500)
                    #print query_name, len(block_list), hit_idx
                    #for i, x in enumerate(block_list):
                    #    print '[%i]' % i
                    #    print x

                    # Parse segment start and stop positions
                    hit_df = block_list[-1]

                    # If V call exists, parse V alignment information
                    vdj_start, vdj_end = None, None
                    if v_call is not None:
                        v_align = hit_df[hit_df[0] == 'V'].iloc[0]
                        db_gen['V_SEQ_START'] = v_align[8]
                        db_gen['V_SEQ_LENGTH'] = int(v_align[9]) - int(v_align[8]) + 1
                        db_gen['V_GERM_START'] = v_align[10]
                        db_gen['V_GERM_LENGTH'] = int(v_align[11]) - int(v_align[10]) + 1
                        db_gen['INDELS'] = 'F' if int(v_align[6]) == 0 else 'T'

                        # Update input sequence positions
                        vdj_start = int(v_align[8]) - 1
                        vdj_end = int(v_align[9])

                    # TODO:  needs to check that the V results are present before trying to determine N1_LENGTH from them.
                    # If D call exists, parse D alignment information
                    if d_call is not None:
                        d_align = hit_df[hit_df[0] == 'D'].iloc[0]
                        db_gen['D_SEQ_START'] = d_align[8]
                        db_gen['N1_LENGTH'] = int(d_align[8]) - int(db_gen['V_SEQ_LENGTH']) - int(db_gen['V_SEQ_START'])
                        db_gen['D_SEQ_LENGTH'] = int(d_align[9]) - int(d_align[8]) + 1
                        db_gen['D_GERM_START'] = d_align[10]
                        db_gen['D_GERM_LENGTH'] = int(d_align[11]) - int(d_align[10]) + 1

                        # Update input sequence positions
                        if vdj_start is None:  vdj_start = int(d_align[8]) - 1
                        vdj_end = int(d_align[9])

                    # TODO:  needs to check that the V results are present before trying to determine N1_LENGTH from them.
                    # If J call exists, parse J alignment information
                    if j_call is not None:
                        j_align = hit_df[hit_df[0] == 'J'].iloc[0]
                        db_gen['J_SEQ_START'] = j_align[8]
                        if d_call is not None:
                            db_gen['N2_LENGTH'] = int(j_align[8]) - \
                                                  int(db_gen['D_SEQ_LENGTH']) - \
                                                  int(db_gen['D_SEQ_START'])
                        else:
                            db_gen['N1_LENGTH'] = int(j_align[8]) - \
                                                  int(db_gen['V_SEQ_LENGTH']) - \
                                                  int(db_gen['V_SEQ_START'])
                        db_gen['J_SEQ_LENGTH'] = int(j_align[9]) - int(j_align[8]) + 1
                        db_gen['J_GERM_START'] = j_align[10]
                        db_gen['J_GERM_LENGTH'] = int(j_align[11]) - int(j_align[10]) + 1

                        # Update input sequence positions
                        if vdj_start is None:  vdj_start = int(j_align[8]) - 1
                        vdj_end = int(j_align[9])

                    # Set VDJ sequence
                    if vdj_start is not None and vdj_end is not None:
                        db_gen['SEQUENCE_VDJ'] = db_gen['SEQUENCE_INPUT'][vdj_start:vdj_end]
                    else:
                        db_gen['SEQUENCE_VDJ'] = 'None'

                yield IgRecord(db_gen)


# TODO:  should be more readable
def readIMGT(imgt_files):
    """
    Reads IMGT/HighV-Quest output

    Arguments: 
    imgt_files = IMGT/HighV-Quest output files 1, 2, 3, and 6
        
    Returns: 
    a generator of dictionaries containing alignment data
    """
    imgt_iters = [csv.DictReader(open(f, 'rU'), delimiter='\t') for f in imgt_files]
    # Create a dictionary for each sequence alignment and yield its generator
    for sm, gp, nt, jn in izip(*imgt_iters):
        if "No results" not in sm['Functionality']:
            db_gen = {'SEQUENCE_ID':       sm['Sequence ID'],
                      'SEQUENCE_INPUT':    sm['Sequence'],
                      'FUNCTIONAL':        ['?','T','F'][('productive' in sm['Functionality']) +
                                                         ('unprod' in sm['Functionality'])],
                      'IN_FRAME':          ['?','T','F'][('in-frame' in sm['JUNCTION frame']) +
                                                         ('out-of-frame' in sm['JUNCTION frame'])],
                      'STOP':              ['F','?','T'][('stop codon' in sm['Functionality comment']) +
                                                         ('unprod' in sm['Functionality'])],
                      'MUTATED_INVARIANT': ['F','?','T'][(any(('missing' in sm['Functionality comment'],
                                                               'missing' in sm['V-REGION potential ins/del']))) +
                                                         ('unprod' in sm['Functionality'])],
                      'INDELS':            ['F','T'][any((sm['V-REGION potential ins/del'],
                                                          sm['V-REGION insertions'],
                                                          sm['V-REGION deletions']))],
                      'V_CALL':            re.sub( '\sor\s', ',', re.sub(',','',gp['V-GENE and allele']) ),
                      'D_CALL':            re.sub( '\sor\s', ',', re.sub(',','',gp['D-GENE and allele']) ),
                      'J_CALL':            re.sub( '\sor\s', ',', re.sub(',','',gp['J-GENE and allele']) ),
                      'SEQUENCE_VDJ':      nt['V-D-J-REGION'] if nt['V-D-J-REGION'] else nt['V-J-REGION'],
                      'SEQUENCE_IMGT':     gp['V-D-J-REGION'] if gp['V-D-J-REGION'] else gp['V-J-REGION'],
                      'V_SEQ_START':       nt['V-REGION start'],
                      'V_SEQ_LENGTH':      len(nt['V-REGION']) if nt['V-REGION'] else 0,
                      'V_GERM_START':      1,
                      'V_GERM_LENGTH':     len(gp['V-REGION']) if gp['V-REGION'] else 0,
                      'N1_LENGTH':         sum(int(i) for i in [jn["P3'V-nt nb"],
                                                                jn['N-REGION-nt nb'],
                                                                jn['N1-REGION-nt nb'],
                                                                jn["P5'D-nt nb"]] if i),
                      'D_SEQ_START':       sum(int(i) for i in [1, len(nt['V-REGION']),
                                                                jn["P3'V-nt nb"],
                                                                jn['N-REGION-nt nb'],
                                                                jn['N1-REGION-nt nb'],
                                                                jn["P5'D-nt nb"]] if i),
                      'D_SEQ_LENGTH':      int(jn["D-REGION-nt nb"] or 0),
                      'D_GERM_START':      int(jn["5'D-REGION trimmed-nt nb"] or 0) + 1,
                      'D_GERM_LENGTH':     int(jn["D-REGION-nt nb"] or 0),
                      'N2_LENGTH':         sum(int(i) for i in [jn["P3'D-nt nb"],
                                                                jn['N2-REGION-nt nb'],
                                                                jn["P5'J-nt nb"]] if i),
                      'J_SEQ_START':       sum(int(i) for i in [1, len(nt['V-REGION']),
                                                                jn["P3'V-nt nb"],
                                                                jn['N-REGION-nt nb'],
                                                                jn['N1-REGION-nt nb'],
                                                                jn["P5'D-nt nb"],
                                                                jn["D-REGION-nt nb"],
                                                                jn["P3'D-nt nb"],
                                                                jn['N2-REGION-nt nb'],
                                                                jn["P5'J-nt nb"]] if i),
                      'J_SEQ_LENGTH':      len(nt['J-REGION']) if nt['J-REGION'] else 0,
                      'J_GERM_START':      int(jn["5'J-REGION trimmed-nt nb"] or 0) + 1,
                      'J_GERM_LENGTH':     len(gp['J-REGION']) if gp['J-REGION'] else 0,
                      'JUNCTION_LENGTH':   len(jn['JUNCTION']) if jn['JUNCTION'] else 0,
                      'JUNCTION':          jn['JUNCTION']}
        else:
            db_gen = {'SEQUENCE_ID':sm['Sequence ID'],
                      'SEQUENCE_INPUT':sm['Sequence'],
                      'V_CALL':'None', 'D_CALL':'None', 'J_CALL':'None'}
        yield IgRecord(db_gen)

    
def getIDforIMGT(seq_file):
    """
    Create a sequence ID translation using IMGT truncation
    
    Arguments: 
    seq_file = a fasta file of sequences input to IMGT
                    
    Returns: 
    a dictionary of {truncated ID: full seq description} 
    """
    
    # Create a seq_dict ID translation using IDs truncate up to space or 50 chars
    ids = {}
    for i, rec in enumerate(SeqIO.parse(seq_file, 'fasta', IUPAC.ambiguous_dna)):
        if len(rec.description) <= 50:
            id_key = rec.description
        else:
            id_key = re.sub('\||\s|!|&|\*|<|>|\?','_',rec.description[:50])
        ids.update({id_key:rec.description})

    return ids


def writeDb(db_gen, file_prefix, total_count, id_dict={}, no_parse=True,
            out_args=default_out_args):
    """
    Writes tab-delimited database file in output directory
    
    Arguments:
    db_gen = a generator of IgRecord objects containing alignment data
    file_prefix = directory and prefix for CLIP tab-delim file
    total_count = number of records (for progress bar)
    id_dict = a dictionary of {IMGT ID: full seq description}
    no_parse = if ID is to be parsed for pRESTO output with default delimiters
    out_args = common output argument dictionary from parseCommonArgs

    Returns:
    None
    """
    pass_file = "%s_db-pass.tab" % file_prefix
    fail_file = "%s_db-fail.tab" % file_prefix
    ordered_fields = ['SEQUENCE_ID',
                      'SEQUENCE_INPUT',
                      'FUNCTIONAL',
                      'IN_FRAME',
                      'STOP',
                      'MUTATED_INVARIANT',
                      'INDELS',
                      'V_CALL',
                      'D_CALL',
                      'J_CALL',
                      'SEQUENCE_VDJ',
                      'SEQUENCE_IMGT',
                      'V_SEQ_START',
                      'V_SEQ_LENGTH',
                      'V_GERM_START',
                      'V_GERM_LENGTH',
                      'N1_LENGTH',
                      'D_SEQ_START',
                      'D_SEQ_LENGTH',
                      'D_GERM_START',
                      'D_GERM_LENGTH',
                      'N2_LENGTH',
                      'J_SEQ_START',
                      'J_SEQ_LENGTH',
                      'J_GERM_START',
                      'J_GERM_LENGTH',
                      'JUNCTION_LENGTH',
                      'JUNCTION']

    # TODO:  This is not the best approach. should pass in output fields.
    # Open passed file
    #pass_handle = open(pass_file, 'wb')
    #pass_writer = getDbWriter(pass_handle, add_fields=ordered_fields)

    # Open failed file
    if out_args['failed']:
        fail_handle = open(fail_file, 'wb')
        fail_writer = getDbWriter(fail_handle, add_fields=['SEQUENCE_ID','SEQUENCE_INPUT'])
    else:
        fail_handle = None
        fail_writer = None

    # Initialize counters and file
    pass_writer = None
    start_time = time()
    rec_count = pass_count = fail_count = 0
    for record in db_gen:
        #printProgress(i + (total_count/2 if id_dict else 0), total_count, 0.05, start_time)
        printProgress(rec_count, total_count, 0.05, start_time)
        rec_count += 1

        # Count pass or fail
        if (record.v_call == 'None' and record.j_call == 'None') or \
                        record.functional is None or not record.seq_vdj or not record.junction:
            fail_count += 1
            if fail_writer is not None: fail_writer.writerow(record.toDict())
            continue
        else: 
            pass_count += 1
            record.seq_input = (record.seq_input.upper() if record.seq_input else '')
            record.seq_imgt = (record.seq_imgt.upper() if record.seq_imgt else '')
            record.seq_vdj = (record.seq_vdj.upper() if record.seq_vdj else '')
            record.junction = (record.junction.upper() if record.junction else '')
            
        # Build sample sequence description
        if record.id in id_dict:
            record.id = id_dict[record.id]

        # Parse sequence description into new columns
        if not no_parse:
            record.annotations = parseAnnotation(record.id, delimiter=out_args['delimiter'])
            record.id = record.annotations['ID']
            del record.annotations['ID']

        # TODO:  This is not the best approach. should pass in output fields.
        # If first sequence, use parsed description to create new columns and initialize writer
        if pass_writer is None:
            if not no_parse:  ordered_fields.extend(record.annotations.keys())
            pass_handle = open(pass_file, 'wb')
            pass_writer = getDbWriter(pass_handle, add_fields=ordered_fields)

        # Write row to tab-delim CLIP file
        pass_writer.writerow(record.toDict())
    
    # Print log
    #printProgress(i+1 + (total_count/2 if id_dict else 0), total_count, 0.05, start_time)
    printProgress(rec_count, total_count, 0.05, start_time)

    log = OrderedDict()
    log['OUTPUT'] = pass_file
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    log['END'] = 'MakeDb'
    printLog(log)
    
    pass_handle.close()
    if fail_handle is not None: fail_handle.close()


# TODO:  may be able to merge with parseIMGT
def parseIgBlast(igblast_output, seq_file, no_parse=True, out_args=default_out_args):
    """
    Main for IgBlast aligned sample sequences

    Arguments:
    igblast_output = IgBlast output file to process
    seq_file = fasta file input to IgBlast (from which to get sequence)
    no_parse = if ID is to be parsed for pRESTO output with default delimiters
    out_args = common output argument dictionary from parseCommonArgs

    Returns:
    None
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'MakeDB'
    log['ALIGNER'] = 'IgBlast'
    log['ALIGN_RESULTS'] = os.path.basename(igblast_output)
    log['SEQ_FILE'] = os.path.basename(seq_file)
    log['NO_PARSE'] = no_parse
    printLog(log)

    # Get input sequence dictionary
    seq_dict = getSeqforIgBlast(seq_file)

    # Formalize out_dir and file-prefix
    if not out_args['out_dir']:
        out_dir = os.path.split(igblast_output)[0]
    else:
        out_dir = os.path.abspath(out_args['out_dir'])
        if not os.path.exists(out_dir):  os.mkdir(out_dir)
    if out_args['out_name']:
        file_prefix = out_args['out_name']
    else:
        file_prefix = os.path.basename(os.path.splitext(igblast_output)[0])
    file_prefix = os.path.join(out_dir, file_prefix)

    total_count = countSeqFile(seq_file)

    # Create
    igblast_dict = readIgBlast(igblast_output, seq_dict)
    writeDb(igblast_dict, file_prefix, total_count,
            no_parse=no_parse, out_args=out_args)


# TODO:  may be able to merge with parseIgBlast
def parseIMGT(imgt_output, seq_file=None, no_parse=True, out_args=default_out_args):
    """
    Main for IMGT aligned sample sequences

    Arguments:
    imgt_output = zipped file or unzipped folder output by IMGT
    seq_file = FASTA file input to IMGT (from which to get seqID)
    no_parse = if ID is to be parsed for pRESTO output with default delimiters
    out_args = common output argument dictionary from parseCommonArgs
        
    Returns: 
    None
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'MakeDb'
    log['ALIGNER'] = 'IMGT'
    log['ALIGN_RESULTS'] = imgt_output
    log['SEQ_FILE'] = os.path.basename(seq_file) if seq_file else ''
    log['NO_PARSE'] = no_parse
    printLog(log)
    
    # Get individual IMGT result files
    temp_dir, imgt_files = extractIMGT(imgt_output)
        
    # Formalize out_dir and file-prefix
    if not out_args['out_dir']:
        out_dir = os.path.dirname(os.path.abspath(imgt_output))
    else:
        out_dir = os.path.abspath(out_args['out_dir'])
        if not os.path.exists(out_dir):  os.mkdir(out_dir)
    if out_args['out_name']:
        file_prefix = out_args['out_name']
    else:
        file_prefix = os.path.splitext(os.path.split(os.path.abspath(imgt_output))[1])[0]
    file_prefix = os.path.join(out_dir, file_prefix)

    total_count = countDbFile(imgt_files[0])
    
    # Get (parsed) IDs from fasta file submitted to IMGT
    id_dict = getIDforIMGT(seq_file) if seq_file else {}
    
    # Create
    imgt_dict = readIMGT(imgt_files)
    writeDb(imgt_dict, file_prefix, total_count, id_dict=id_dict,
            no_parse=no_parse, out_args=out_args)

    # Delete temp directory
    rmtree(temp_dir)


def getArgParser():
    """
    Defines the ArgumentParser

    Arguments: 
    None
                      
    Returns: 
    an ArgumentParser object
    """
    fields = textwrap.dedent(
             '''
              output files:
                db-pass        database of parsed alignment records.
                db-fail        database with records failing alignment.

              output fields:
                SEQUENCE_ID
                SEQUENCE_INPUT
                FUNCTIONAL
                IN_FRAME
                STOP
                MUTATED_INVARIANT
                INDELS
                V_CALL
                D_CALL
                J_CALL
                SEQUENCE_VDJ and/or SEQUENCE_IMGT
                V_SEQ_START
                V_SEQ_LENGTH
                V_GERM_START
                V_GERM_LENGTH
                N1_LENGTH
                D_SEQ_START
                D_SEQ_LENGTH
                D_GERM_START
                D_GERM_LENGTH
                N2_LENGTH
                J_SEQ_START
                J_SEQ_LENGTH
                J_GERM_START
                J_GERM_LENGTH
                JUNCTION_LENGTH
                JUNCTION
              ''')
                
    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields, 
                            version='%(prog)s:' + ' v%s-%s' %(__version__, __date__),  
                            formatter_class=CommonHelpFormatter)
    subparsers = parser.add_subparsers(title='subcommands', dest='command', 
                                       help='Aligner used', metavar='')
    
    # Parent parser    
    parser_parent = getCommonArgParser(seq_in=False, seq_out=False, log=False)

    # IgBlast Aligner
    parser_igblast = subparsers.add_parser('igblast', help='Process IgBlast output',
                                           parents=[parser_parent],
                                           formatter_class=CommonHelpFormatter)
    parser_igblast.set_defaults(func=parseIgBlast)
    parser_igblast.add_argument('-i', nargs='+', action='store', dest='aligner_files',
                                required=True,
                                help='''IgBLAST output files in format 7 (IgBLAST
                                     argument \'-outfmt 7\').''')
    parser_igblast.add_argument('-s', action='store', nargs='+', dest='seq_files',
                                required=True,
                                help='List of input FASTA files containing sequences')
    parser_igblast.add_argument('--noparse', action='store_true', dest='no_parse',
                                help='''Specify if input IDs should not be parsed to add
                                     new columns to database.''')
    
    # IMGT aligner
    parser_imgt = subparsers.add_parser('imgt', help='Process IMGT/HighV-Quest output', 
                                        parents=[parser_parent], 
                                        formatter_class=CommonHelpFormatter)
    imgt_arg_group =  parser_imgt.add_mutually_exclusive_group(required=True)
    imgt_arg_group.add_argument('-i', nargs='+', action='store', dest='aligner_files',
                                help='''Either zipped IMGT output files (.zip) or a folder
                                     containing unzipped IMGT output files (which must
                                     include 1_Summary, 2_IMGT-gapped, 3_Nt-sequences,
                                     and 6_Junction).''')
    parser_imgt.add_argument('-s', nargs='*', action='store', dest='seq_files',
                             required=False,
                             help='List of input FASTA files containing sequences')
    parser_imgt.add_argument('--noparse', action='store_true', dest='no_parse', 
                             help='''Specify if input IDs should not be parsed to add new
                                  columns to database.''')
    parser_imgt.set_defaults(func=parseIMGT)

    return parser
    
    
if __name__ == "__main__":
    """
    Parses command line arguments and calls main
    """
    parser = getArgParser()    
    args = parser.parse_args()
    args_dict = parseCommonArgs(args, in_arg='aligner_files')

    # Set no ID parsing if sequence files are not provided
    if 'seq_files' in args_dict and not args_dict['seq_files']:
        args_dict['no_parse'] = True

    # Delete
    if 'seq_files' in args_dict: del args_dict['seq_files']
    if 'aligner_files' in args_dict: del args_dict['aligner_files']
    if 'command' in args_dict: del args_dict['command']
    if 'func' in args_dict: del args_dict['func']           
    
    # IMGT parser
    if args.command == 'imgt':
        for i in range(len(args.__dict__['aligner_files'])):
            args_dict['imgt_output'] = args.__dict__['aligner_files'][i]
            args_dict['seq_file'] = args.__dict__['seq_files'][i] \
                                    if args.__dict__['seq_files'] else None
            args.func(**args_dict)
    elif args.command == 'igblast':
        for i in range(len(args.__dict__['aligner_files'])):
            args_dict['igblast_output'] =  args.__dict__['aligner_files'][i]
            args_dict['seq_file'] = args.__dict__['seq_files'][i]
            args.func(**args_dict)
