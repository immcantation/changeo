#!/usr/bin/env python
"""
Create tab-delimited database file to store sequence alignment information
"""

__author__    = 'Namita Gupta, Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.0'
__date__      = '2015.04.06'

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
from DbCore import getDbWriter, IgRecord, countDbFile

# Default parameters
default_V_regex = re.compile(r'((IG[HLK]|TR[ABGD])V\d+[-/\w]*[-\*][\.\w]+)')
default_D_regex = re.compile(r'((IG[HLK]|TR[ABGD])D\d+[-/\w]*[-\*][\.\w]+)')
default_J_regex = re.compile(r'((IG[HLK]|TR[ABGD])J\d+[-/\w]*[-\*][\.\w]+)')
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
    imgt_flags = ["1_Summary", "2_IMGT-gapped", "3_Nt-sequences", "6_Junction"]
    temp_dir = mkdtemp()
    if is_zipfile(imgt_output):
        # Extract selected files from the IMGT zip file to temp directory
        imgt_zip = ZipFile(imgt_output, 'r')
        imgt_files = sorted([n for n in imgt_zip.namelist() for i in imgt_flags if i in n])
        imgt_zip.extractall(temp_dir, imgt_files)
        imgt_files = [os.path.join(temp_dir, f) for f in imgt_files]
    elif os.path.isdir(imgt_output):
        imgt_files = sorted([ (imgt_output + ('/' if imgt_output[-1]!='/' else '') + n) \
                             for n in os.listdir(imgt_output) for j in imgt_flags if j in n ])
    else:
        sys.exit('ERROR: Unsupported IGMT output file-type - must be zipped file (-z) or folder (-f)')
    
    if len(imgt_files) > len(imgt_flags): # e.g. multiple 1_Summary files
        sys.exit('ERROR: Wrong files in folder %s' % imgt_output)
    elif len(imgt_files) < len(imgt_flags):
        sys.exit('ERROR: Missing necessary file in folder %s' % imgt_output)
        
    return temp_dir, imgt_files


def readOneIgBlastResult(block):
    """
    Parse a single IgBLAST query result

    :param block: itertools groupby object of single result
    :return: None if no results, otherwise list of DataFrames for each result block
    """
    results = list()
    for i, (match, subblock) in enumerate(groupby(block, lambda l: l=='\n')):
        if not match:
            sub = [s.strip() for s in subblock if not s.startswith('#')]
            if i==2 or i==4:
                results.append(sub[0])
            else:
                sub = [s.split('\t') for s in sub]
                df = pd.DataFrame(sub)
                if not df.empty: results.append(df)
    if not results: return None
    return results


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
                    # Parse V, D, and J calls
                    v_call = IgRecord._parseAllele(block_list[0], default_V_regex, action='list')
                    d_call = IgRecord._parseAllele(block_list[0], default_D_regex, action='list')
                    j_call = IgRecord._parseAllele(block_list[0], default_J_regex, action='list')
                    db_gen['V_CALL'] = ','.join(v_call) if v_call is not None else 'None'
                    db_gen['D_CALL'] = ','.join(d_call) if d_call is not None else 'None'
                    db_gen['J_CALL'] = ','.join(j_call) if j_call is not None else 'None'

                    # Parse quality information
                    quals = block_list[0].split()
                    db_gen['STOP'] = 'T' if quals[-4] == 'Yes' else 'F'
                    db_gen['IN_FRAME'] = 'T' if quals[-3] == 'In-frame' else 'F'
                    db_gen['FUNCTIONAL'] = 'T' if quals[-2] == 'Yes' else 'F'

                    # Parse junction sequence
                    db_gen['JUNCTION'] = re.sub('(N/A)|\[|\(|\)|\]', '',
                                                ''.join(block_list[1].split()))
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


def writeDb(db_gen, no_parse, file_prefix, total_count, out_args,
            id_dict={}):
    """
    Writes tab-delimited database file in output directory
    
    Arguments:
    db_gen = a generator of IgRecord objects containing alignment data
    no_parse = if ID is to be parsed for pRESTO output with default delimiters
    file_prefix = directory and prefix for CLIP tab-delim file
    total_count = number of records (for progress bar)
    out_args = common output argument dictionary from parseCommonArgs
    id_dict = a dictionary of {IMGT ID: full seq description}
    
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

    # TODO:  pass_writer is overwritten later if no_parse=True, which is not the best approach. should pass in output fields.
    # Open passed file
    pass_handle = open(pass_file, 'wb')
    pass_writer = getDbWriter(pass_handle, add_fields=ordered_fields)

    # Open failed file
    if out_args['failed']:
        fail_handle = open(fail_file, 'wb')
        fail_writer = getDbWriter(fail_handle, add_fields=['SEQUENCE_ID','SEQUENCE_INPUT'])
    else:
        fail_handle = None
        fail_writer = None

    # Initialize counters and file
    start_time = time()
    rec_count = pass_count = fail_count = 0
    for i, record in enumerate(db_gen):
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
            record.seq_in = (record.seq_in.upper() if record.seq_in else '')
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

            # If first sequence, use parsed description to create new columns and re-initialize writer
            if i == 0:
                ordered_fields.extend(record.annotations.keys())
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


def parseIgBlast(seq_file, igblast_output, no_parse, out_args=default_out_args):
    """
    Main for IgBlast aligned sample sequences

    Arguments:
    seq_file = fasta file input to IgBlast (from which to get sequence)
    igblast_output = IgBlast output file to process
    no_parse = if ID is to be parsed for pRESTO output with default delimiters
    out_args = common output argument dictionary from parseCommonArgs

    Returns:
    None
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'MakeDB'
    log['ALIGNER'] = 'IgBlast'
    log['SEQ_FILE'] = os.path.basename(seq_file)
    log['ALIGN_RESULTS'] = os.path.basename(igblast_output)
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
    file_prefix = os.path.join( out_dir, file_prefix)

    total_count = countSeqFile(seq_file)
    start_time = time()

    # Create
    igblast_dict = readIgBlast(igblast_output, seq_dict)
    writeDb(igblast_dict, no_parse, file_prefix, start_time, total_count, out_args)
    

def parseIMGT(seq_file, imgt_output, no_parse, out_args=default_out_args):
    """
    Main for IMGT aligned sample sequences

    Arguments: 
    seq_file = FASTA file input to IMGT (from which to get seqID)
    imgt_output = zipped file or unzipped folder output by IMGT
    no_parse = if ID is to be parsed for pRESTO output with default delimiters
    out_args = common output argument dictionary from parseCommonArgs
        
    Returns: 
    None
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'MakeDb'
    log['ALIGNER'] = 'IMGT'
    log['SEQ_FILE'] = os.path.basename(seq_file) if seq_file else ''
    log['ALIGN_RESULTS'] = imgt_output
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
    writeDb(imgt_dict, no_parse, file_prefix, total_count, out_args, id_dict=id_dict)

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
    parser_igblast.add_argument('-i', nargs='+', action='store', dest='aligner_output', required=True,
                                help='IgBLAST output files in format 7 (IgBLAST argument "-outfmt 7").')
    parser_igblast.add_argument('-s', action='store', nargs='+', dest='seq_files', required=True,
                                help='List of input FASTA files containing sequences')
    parser_igblast.add_argument('--noparse', action='store_true', dest='no_parse',
                                help='Specify if input IDs should not be parsed to add new columns to database')
    
    # IMGT aligner
    parser_imgt = subparsers.add_parser('imgt', help='Process IMGT/HighV-Quest output', 
                                        parents=[parser_parent], 
                                        formatter_class=CommonHelpFormatter)
    parser_imgt.set_defaults(func=parseIMGT)
    imgt_arg_group =  parser_imgt.add_mutually_exclusive_group(required=True)
    imgt_arg_group.add_argument('-z', nargs='+', action='store', dest='aligner_output',
                                help='Zipped IMGT output files')
    imgt_arg_group.add_argument('-f', nargs='+', action='store', dest='aligner_output', 
                                help='Folder with unzipped IMGT output files \
                                     (must have 1_Summary, 2_IMGT-gapped, 3_Nt-sequences, and 6_Junction)')
    parser_imgt.add_argument('-s', action='store', nargs='+', dest='seq_files',
                             help='List of input FASTA files containing sequences')
    parser_imgt.add_argument('--noparse', action='store_true', dest='no_parse', 
                             help='Specify if input IDs should not be parsed to add new columns to database')
    
    return parser
    
    
if __name__ == "__main__":
    """
    Parses command line arguments and calls main
    """
    parser = getArgParser()    
    args = parser.parse_args()
    args_dict = parseCommonArgs(args, in_arg='aligner_output')
    
    if 'seq_files' in args_dict: del args_dict['seq_files']
    else: args_dict['no_parse'] = True
    if 'aligner_output' in args_dict: del args_dict['aligner_output']
    if 'command' in args_dict: del args_dict['command']
    if 'func' in args_dict: del args_dict['func']           
    
    # IMGT parser
    if args.command == 'imgt':
        if args.__dict__['aligner_output']:
            for i in range(len(args.__dict__['aligner_output'])):
                args_dict['seq_file'] = args.__dict__['seq_files'][i] if args.__dict__['seq_files'] else None
                args_dict['imgt_output'] = args.__dict__['aligner_output'][i]
                args.func(**args_dict)
                # TODO: figure out how to delete extracted zip files safely
                # if is_zipfile(args_dict['imgt_output']):
                #     rmtree(os.path.splitext(args_dict['imgt_output'])[0])
        else:
            parser.error('Must include either (-z) zipped IMGT files or \
                         (-f) folder with individual files 1_, 2_, 3_, and 6_')
    elif args.command == 'igblast':
        if args.__dict__['aligner_output']:
            for i in range(len(args.__dict__['aligner_output'])):
                args_dict['seq_file'] = args.__dict__['seq_files'][i] if args.__dict__['seq_files'] else \
                                        parser.error('Must include fasta file input to IgBlast')
                args_dict['igblast_output'] =  args.__dict__['aligner_output'][i]
                args.func(**args_dict)
        else:
            parser.error('Must include IgBlast output file (-o)')