#!/usr/bin/env python
"""
Create tab-delimited database file to store sequence alignment information
Output Columns: SEQUENCE_ID, SEQUENCE, FUNCTIONAL, IN_FRAME, STOP, MUTATED_INVARIANT, INDELS,
                V_MATCH, V_LENGTH, J_MATCH, J_LENGTH, V_CALL, D_CALL, J_CALL, SEQUENCE_GAP,
                V_SEQ_START, V_SEQ_LENGTH, V_GERM_START, V_GERM_LENGTH, N1_LENGTH, D_SEQ_START,
                D_SEQ_LENGTH, D_GERM_START, D_GERM_LENGTH, N2_LENGTH, J_SEQ_START, J_SEQ_LENGTH,
                J_GERM_START, J_GERM_LENGTH, JUNCTION_GAP_LENGTH, JUNCTION 
"""

__author__    = 'Namita Gupta, Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.0'
__date__      = '2014.9.4'

# Imports
import re
import csv
import os, sys
from zipfile import ZipFile, is_zipfile
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from itertools import izip
from collections import OrderedDict
from time import time

# IgCore and DbCore imports 
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from IgCore import default_out_args, parseAnnotation, printLog, printProgress
from IgCore import getCommonArgParser, parseCommonArgs, countSeqFile
from DbCore import getDbWriter, IgRecord, countDbFile

# Default parameters
default_V_regex = re.compile(r"(IG[HLK][V]\d+[-/\w]*[-\*][\.\w]+)")
default_D_regex = re.compile(r"(IG[HLK][D]\d+[-/\w]*[-\*][\.\w]+)")
default_J_regex = re.compile(r"(IG[HLK][J]\d+[-/\w]*[-\*][\.\w]+)")
default_delimiter = ('\t', ',', '-')


def extractIMGT(imgt_output):
    """
    Extract necessary files from IMGT results, zipped or unzipped
    
    Arguments:
    imgt_output = zipped file or unzipped folder output by IMGT
    
    Returns:
    sorted list of filenames from which information will be read
    """
    imgt_flags = ["1_Summary", "2_IMGT-gapped", "3_Nt-sequences", "6_Junction"]
    if is_zipfile(imgt_output):
        # Extract selected files from the IMGT zip file
        # Try to retain compressed format
        imgt_zip = ZipFile(imgt_output, 'r')
        imgt_files = sorted([n for n in imgt_zip.namelist() for i in imgt_flags if i in n])
        for f in imgt_files:  imgt_zip.extract(f)
    elif os.path.isdir(imgt_output):
        imgt_files = sorted([ (imgt_output + ('/' if imgt_output[-1]!='/' else '') + n) \
                             for n in os.listdir(imgt_output) for j in imgt_flags if j in n ])
    else:
        sys.exit('ERROR: Unsupported IGMT output file-type - must be zipped file (-z) or folder (-f)')
    
    if len(imgt_files) > len(imgt_flags): # e.g. multiple 1_Summary files
        sys.exit('ERROR: Wrong files in folder %s' % imgt_output)
    elif len(imgt_files) < len(imgt_flags):
        sys.exit('ERROR: Missing necessary file in folder %s' % imgt_output)
        
    return imgt_files


def readIMGT(imgt_files):
    """
    Reads IMGT/HighV-Quest output

    Arguments: 
    imgt_files = IMGT/HighV-Quest output files 1, 2, 3, and 6
        
    Returns: 
    a generator of dictionaries containing alignment data
    """
    imgt_iters = [csv.DictReader(open(f, 'rU'), delimiter='\t') for f in imgt_files]
    # Create a generator of dictionaries for each sequence alignment
    db_gen = (IgRecord(
              {'SEQUENCE_ID':        sm['Sequence ID'],
               'SEQUENCE':           sm['Sequence'],
               'FUNCTIONAL':         ['?','T','F'][('productive' in sm['Functionality']) + ('unprod' in sm['Functionality'])],
               'IN_FRAME':           ['?','T','F'][('in-frame' in sm['JUNCTION frame']) + ('out-of-frame' in sm['JUNCTION frame'])],
               'STOP':               ['F','?','T'][('stop codon' in sm['Functionality comment']) + ('unprod' in sm['Functionality'])],
               'MUTATED_INVARIANT':  ['F','?','T'][(any(('missing' in sm['Functionality comment'],
                                                       'missing' in sm['V-REGION potential ins/del']))) + ('unprod' in sm['Functionality'])],
               'INDELS':             ['F','T'][any((sm['V-REGION potential ins/del'], 
                                                    sm['V-REGION insertions'], 
                                                    sm['V-REGION deletions']))],
               'V_MATCH':            0 if sm['V-REGION identity nt'] == 'null' or not sm['V-REGION identity nt'] \
                                        else int(sm['V-REGION identity nt'].split('/')[0] or 0) + nt['V-REGION'].count('n'),
               'V_LENGTH':           0 if sm['V-REGION identity nt'] == 'null' or not sm['V-REGION identity nt'] \
                                        else int(sm['V-REGION identity nt'].split('/')[1].split()[0]),
               'J_MATCH':            0 if sm['J-REGION identity nt'] == 'null' or not sm['J-REGION identity nt'] \
                                        else int(sm['J-REGION identity nt'].split('/')[0] or 0) + nt['J-REGION'].count('n'),
               'J_LENGTH':           0 if sm['J-REGION identity nt'] == 'null' or not sm['J-REGION identity nt'] \
                                        else int(sm['J-REGION identity nt'].split('/')[1].split()[0]),
               'V_CALL':             re.sub( '\sor\s', ',', re.sub(',','',gp['V-GENE and allele']) ), # replace or with comma
               'D_CALL':             re.sub( '\sor\s', ',', re.sub(',','',gp['D-GENE and allele']) ),
               'J_CALL':             re.sub( '\sor\s', ',', re.sub(',','',gp['J-GENE and allele']) ),
               'SEQUENCE_GAP':       gp['V-D-J-REGION'] if gp['V-D-J-REGION'] else gp['V-J-REGION'],
               'V_SEQ_START':        1,
               'V_SEQ_LENGTH':       len(nt['V-REGION']) if nt['V-REGION'] else 0,
               'V_GERM_START':       1,
               'V_GERM_LENGTH':      len(gp['V-REGION']) if gp['V-REGION'] else 0,
               'N1_LENGTH':          sum(int(i) for i in [jn["P3'V-nt nb"], 
                                                          jn['N-REGION-nt nb'], 
                                                          jn['N1-REGION-nt nb'], 
                                                          jn["P5'D-nt nb"]] if i),
               #'D_5_TRIM':           int(jn["5'D-REGION trimmed-nt nb"] or 0),
               'D_SEQ_START':        sum(int(i) for i in [1, len(nt['V-REGION']),
                                                          jn["P3'V-nt nb"], 
                                                          jn['N-REGION-nt nb'], 
                                                          jn['N1-REGION-nt nb'], 
                                                          jn["P5'D-nt nb"]] if i),
               'D_SEQ_LENGTH':       int(jn["D-REGION-nt nb"] or 0),
               'D_GERM_START':       int(jn["5'D-REGION trimmed-nt nb"] or 0) + 1,
               'D_GERM_LENGTH':      int(jn["D-REGION-nt nb"] or 0),
               'N2_LENGTH':          sum(int(i) for i in [jn["P3'D-nt nb"],
                                                          jn['N2-REGION-nt nb'],
                                                          jn["P5'J-nt nb"]] if i),   
               #'J_5_TRIM':           int(jn["5'J-REGION trimmed-nt nb"] or 0),
               'J_SEQ_START':        sum(int(i) for i in [1, len(nt['V-REGION']),
                                                          jn["P3'V-nt nb"], 
                                                          jn['N-REGION-nt nb'], 
                                                          jn['N1-REGION-nt nb'], 
                                                          jn["P5'D-nt nb"],
                                                          jn["D-REGION-nt nb"],
                                                          jn["P3'D-nt nb"],
                                                          jn['N2-REGION-nt nb'],
                                                          jn["P5'J-nt nb"]] if i),
               'J_SEQ_LENGTH':       len(nt['J-REGION']) if nt['J-REGION'] else 0,
               'J_GERM_START':       int(jn["5'J-REGION trimmed-nt nb"] or 0) + 1,
               'J_GERM_LENGTH':      len(gp['J-REGION']) if gp['J-REGION'] else 0,
               'JUNCTION_GAP_LENGTH': len(jn['JUNCTION']) if jn['JUNCTION'] else 0,
               'JUNCTION':           jn['JUNCTION']}) if "No results" not in sm['Functionality'] else \
              IgRecord({'SEQUENCE_ID':sm['Sequence ID'], 'SEQUENCE':sm['Sequence'],
                        'V_CALL':'None', 'D_CALL':'None', 'J_CALL':'None'})
              for sm, gp, nt, jn in izip(*imgt_iters) )
    
    return db_gen

    
def getIDforIMGT(seq_file, start_time, total_count, out_args, id_only=False):
    """
    Create a sequence ID translation using IMGT truncation
    
    Arguments: 
    seq_file = a fasta file of sequences input to IMGT
    start_time = time from which to count elapsed time
    total_count = number of records (for progress bar)
    out_args = common output argument dictionary from parseCommonArgs
    id_only = flag whether only sequence ID 
              (not full description) was used for IMGT input
                    
    Returns: 
    a dictionary of {truncated ID: full seq description} 
    """
    
    seq_dict = SeqIO.index(seq_file, "fasta", IUPAC.ambiguous_dna)
    
    # Create a seq_dict ID translation using IDs truncate up to space or 50 chars
    ids = {}
    for i,seq in enumerate(seq_dict.itervalues()):
        if id_only:
            id_key = parseAnnotation(seq.description, fields=['ID'], delimiter=out_args['delimiter'])['ID']
        else:
            id_key = re.sub('\||\s','_',seq.description[:50])
        ids.update({id_key:seq.description})
        printProgress(i, total_count, 0.05, start_time)
    return ids


def getGapDict(vgap_file, delimiter=default_delimiter):
    """
    Gets dictionary of gap locations and length by V allele from file
    
    Arguments:
    vgap_file = file with gapping information
    delimiter = a tuple of delimiters for between (ID and gaps, gaps, and gap start and gap length)
    
    Returns:
    Dictionary of {V_GL: [list of (gap_start, gap_length),...] } 
    """
    gap_dict = {}
    gap_handle = open(vgap_file, 'rU')
    line = gap_handle.readline()
    while line:
        w = line.split(delimiter[0])
        gaps = w[1].split(delimiter[1])
        
        # Make each gap a tuple of (start, length)
        gaps = [tuple(gap.split(delimiter[2])) for gap in gaps]
        gap_dict[tuple(IgRecord.allele_regex.findall(w[0]))] = gaps
        
        line = gap_handle.readline()
        
    return gap_dict


def getJaaDict(jaa_file, delimiter=default_delimiter):
    """
    Gets dictionary of location of conserved aa by J allele from file
    
    Arguments:
    jaa_file = file with conserved aa information for each J-allele
    delimiter = a tuple of delimiters for between (ID and gaps, gaps, and gap start and gap length)
    
    Returns:
    Dictionary of {J_GL: conserved aa position}
    """
    j_dict = {}
    j_handle = open(jaa_file, 'rU')
    line = j_handle.readline()
    while line:
        w = line.split(delimiter[0])
        j_dict[tuple(IgRecord.allele_regex.findall(w[0]))] = w[1]
        line = j_handle.readline()
        
    return j_dict


def gapIgBlastQuery(record, gap_dict, j_dict):
    """
    Adds gapped sequences and IMGT junctions to IgBlast aligned database files

    Arguments:
    record = dictionary with one IgBlast result
    gap_dict = dictionary of {IgBlast_GL: [list of (gap_start, gap_length),...] }
    j_dict = dictionary of {J_GL: conserved aa start}

    Returns:
    dictionary with gapped IgBlast result and IMGT junction
    """
    v_call_col = 'V_CALL'
    if (record.v_call == 'None' and record.j_call == 'None') or record.functional is None:
        return record
    else:
        record = record.toDict()
 
    seq = (int(record['V_GERM_START'] or 0)-1)*'.' + \
           record['SEQUENCE'][(int(record['V_SEQ_START'] or 0)-1):(int(record['V_SEQ_LENGTH'] or 0)-int(record['V_SEQ_START'] or 0))]
    v_call = record[v_call_col]
    if v_call in gap_dict:
        gap_seq = seq
        for gap_event in gap_dict[v_call]:
            gap_seq = gap_seq[:gap_event[0]] + gap_event[1]*'.' + gap_seq[gap_event[0]:]
        record['V_GERM_LENGTH'] = len(gap_seq)
        gap_seq = gap_seq + record['JUNCTION']
        pre_j_len = len(gap_seq)
        gap_seq = gap_seq + seq[int(record['J_START_SEQ'] or 0):(int(record['J_SEQ_LENGTH'] or 0)-int(record['J_SEQ_START'] or 0))]
    else: 
        gap_seq = ''
        record['V_GERM_LENGTH'] = int(record['V_SEQ_GERM'] or 0)

    record['N1_LENGTH'] = int(record['D_SEQ_START'] or 0) - \
                        int(record['V_SEQ_END'] or 0) - 1
    record['N2_LENGTH'] = int(record['J_SEQ_START'] or 0) - \
                        int(record['D_SEQ_END'] or 0) - 1
                        
    if(record['J_CALL'] in j_dict):
        record['JUNCTION'] = gap_seq[312:(pre_j_len+j_dict[record['J_CALL']]-record['J_GERM_START']+3)]
        record['JUNCTION_GAP_LENGTH'] = len(record['JUNCTION'])
                            
    return IgRecord(record)


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


def readIgBlast(igblast_output, seq_dict):
    """
    Reads IgBlast output

    Arguments: ma
    igblast_output = IgBlast output file (format 7)
    seq_dict = a dictionary of {ID:Seq} from input fasta file
        
    Returns: 
    a generator of dictionaries containing alignment data
    """
    db_gen = {}
    igblast_handle = open(igblast_output, 'rU')
    
    for line in igblast_handle:
        if(re.match("# Query:", line)):
            if 'SEQUENCE_ID' in db_gen and db_gen['FUNCTIONAL'] != 'No Results': yield db_gen
            words = line.split()
            db_gen = {}
            db_gen['SEQUENCE_ID'] = ' '.join(words[2:])
            db_gen['SEQUENCE'] = seq_dict[db_gen['SEQUENCE_ID']]
        if(re.match("# 0 hits found", line)):
            db_gen = IgRecord({'SEQUENCE_ID':db_gen['SEQUENCE_ID'], 'SEQUENCE':db_gen['SEQUENCE'],
                               'V_CALL':'None', 'D_CALL':'None', 'J_CALL':'None'})
            yield db_gen
        if(re.match(r"# V-\(D\)-J rearrangement summary", line)):
            line = next(igblast_handle)
            words = line.split()
            db_gen['V_CALL'] = ','.join(re.findall(default_V_regex, line))
            db_gen['D_CALL'] = ','.join(re.findall(default_D_regex, line))
            db_gen['J_CALL'] = ','.join(re.findall(default_J_regex, line))
            cnt = 4 if db_gen['D_CALL'] else 3
            if(words[cnt]=='No'): 
                db_gen['STOP'] = 'F'
            elif(words[cnt]=='Yes'): 
                db_gen['STOP'] = 'T'
            else: '?'
            cnt += 1
            if(words[cnt]=='In-frame'): 
                db_gen['IN_FRAME'] = 'T'
            elif(words[cnt]=='Out-of-frame'): 
                db_gen['IN_FRAME'] = 'F'
            else: db_gen['IN_FRAME'] = '?'
            cnt += 1
            if(words[cnt]=='No'): 
                db_gen['FUNCTIONAL'] = 'F'
            elif(words[cnt]=='Yes'): 
                db_gen['FUNCTIONAL'] = 'T'
            else: db_gen['FUNCTIONAL'] = '?'
        if(re.match(r"# V-\(D\)-J junction", line)):
            line = next(igblast_handle)
            db_gen['JUNCTION'] = re.sub("(N/A)|\[|\(|\)|\]",'',''.join(line.split()))
            db_gen['JUNCTION_GAP_LENGTH'] = len(db_gen['JUNCTION'])
        if(re.match(r"# Hit table", line)):
            if('V_CALL' in db_gen and db_gen['V_CALL']):
                    line = findLine(igblast_handle,"V")
                    words = line.split()
                    db_gen['V_LENGTH'] = words[4]
                    db_gen['V_MATCH'] = str(int(words[4]) - int(words[5]))
                    db_gen['V_SEQ_START'] =  words[8]
                    db_gen['V_SEQ_LENGTH'] =  words[9] - words[8] + 1
                    db_gen['V_GERM_START'] =  words[10]
                    db_gen['V_GERM_LENGTH'] =  words[11] - words[10] + 1
                    if(words[6] == '0'):
                        db_gen['INDELS'] = 'F'
                    else: db_gen['INDELS'] = 'T'
            if('D_CALL' in db_gen and db_gen['D_CALL']):
                    line = findLine(igblast_handle,"D")
                    words = line.split()
                    db_gen['D_LENGTH'] = words[4]
                    db_gen['D_MATCH'] = str(int(words[4]) - int(words[5]))
                    db_gen['D_SEQ_START'] = words[8]
                    db_gen['D_SEQ_LENGTH'] = words[9] - words[8] + 1
                    db_gen['D_GERM_START'] = words[10]
                    db_gen['D_GERM_LENGTH'] = words[11] - words[10] + 1
            if('J_CALL' in db_gen and db_gen['J_CALL']):
                    line = findLine(igblast_handle,"J")
                    words = line.split()
                    db_gen['J_LENGTH'] = words[4]
                    db_gen['J_MATCH'] = str(int(words[4]) - int(words[5]))
                    db_gen['J_SEQ_START'] = words[8] - words[8] + 1
                    db_gen['J_SEQ_LENGTH'] = words[9]
                    db_gen['J_GERM_START'] = words[10]
                    db_gen['J_GERM_LENGTH'] = words[11] - words[10] + 1
    igblast_handle.close()
    yield IgRecord(db_gen)


def writeDb(db_gen, no_parse, file_prefix, aligner, start_time, total_count, out_args, id_dict={}, seq_dict={}, gap_dict={}, j_dict={}):
    """
    Writes tab-delimited database file in output directory
    
    Arguments:
    db_gen = a generator of IgRecord objects containing alignment data
    no_parse = if ID is to be parsed for pRESTO output with default delimiters
    file_prefix = directory and prefix for CLIP tab-delim file
    aligner = which aligner was used
    start_time = time from which to count elapsed time
    total_count = number of records (for progress bar)
    out_args = common output argument dictionary from parseCommonArgs
    id_dict = a dictionary of {truncated ID: full seq description}
    seq_dict = a dictionary of {ID:Seq} from input fasta file
    gap_dict = dictionary of {IgBlast_GL: [list of (gap_start, gap_length),...] }
    j_dict = dictionary of {J_GL: conserved aa start}
    
    Returns:
    None
    """
    pass_file = "%s_db-pass.tab" % file_prefix
    fail_file = "%s_db-fail.tab" % file_prefix
    ordered_fields = ['SEQUENCE_ID','SEQUENCE','FUNCTIONAL','IN_FRAME','STOP','MUTATED_INVARIANT','INDELS',
                      'V_MATCH','V_LENGTH','J_MATCH','J_LENGTH','V_CALL','D_CALL','J_CALL','SEQUENCE_GAP',
                      'V_SEQ_START','V_SEQ_LENGTH','V_GERM_START','V_GERM_LENGTH','N1_LENGTH','D_SEQ_START',
                      'D_SEQ_LENGTH','D_GERM_START','D_GERM_LENGTH','N2_LENGTH','J_SEQ_START','J_SEQ_LENGTH',
                      'J_GERM_START','J_GERM_LENGTH','JUNCTION_GAP_LENGTH','JUNCTION']
   
    pass_handle = open(pass_file, 'wb')
    fail_handle = open(fail_file, 'wb') if not out_args['clean'] else None
    # Create DbWriter
    if not no_parse:
        for v in id_dict.itervalues():
            tmp = parseAnnotation(v, delimiter=out_args['delimiter'])
            del tmp['ID']
            ordered_fields.extend(tmp.keys())
            break
    pass_writer = getDbWriter(pass_handle, add_fields=ordered_fields)
    fail_writer = getDbWriter(fail_handle, add_fields=['SEQUENCE_ID','SEQUENCE']) if not out_args['clean'] else None
    # Initialize counters
    pass_count = fail_count = 0
    
    for i,record in enumerate(db_gen):
        printProgress(i + (total_count/2 if id_dict else 0), total_count, 0.05, start_time)
        
        # Gap sequence and form IMGT-defined junction region
        if aligner=='igblast':
            db_gen =  gapIgBlastQuery(db_gen, gap_dict, j_dict)
        
        # Count pass or fail
        if (record.v_call == 'None' and record.j_call == 'None') or record.functional is None or not record.seq_gap or not record.junction: 
            fail_count += 1
            if fail_writer is not None: fail_writer.writerow(record.toDict())
            continue
        else: 
            pass_count += 1
            record.seq = record.seq.upper()
            record.seq_gap = record.seq_gap.upper()
            record.junction = record.junction.upper() 
            
            
        # Build sample sequence description
        if record.id.split(' ')[0] in id_dict:
            record.id = id_dict[record.id]
        # Parse sequence description into new columns
        if not no_parse:
            record.annotations = parseAnnotation(record.id, delimiter=out_args['delimiter'])
            record.id = record.annotations['ID']
            del record.annotations['ID']
            
        # Write row to tab-delim CLIP file
        pass_writer.writerow(record.toDict())
    
    # Print log
    printProgress(i+1 + (total_count/2 if id_dict else 0), total_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT'] = pass_file
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    log['END'] = 'MakeDb'
    printLog(log)
    
    pass_handle.close()
    if fail_handle is not None: fail_handle.close()
    

def parseIMGT(seq_file, imgt_output, id_only, no_parse, out_args=default_out_args):
    """
    Main for IMGT aligned sample sequences

    Arguments: 
    seq_file = FASTA file input to IMGT (from which to get seqID)
    imgt_output = zipped file or unzipped folder output by IMGT
    id_only = whether only the sequence ID (with no pRESTO information) was passed to IMGT
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
    log['ID_ONLY'] = id_only 
    log['NO_PARSE'] = no_parse
    printLog(log)
    
    # Get individual IMGT result files
    imgt_files = extractIMGT(imgt_output)
        
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
    
    total_count = countDbFile(imgt_files[0]) * (2 if not no_parse else 1)
    start_time = time()
    
    # Get (parsed) IDs from fasta file submitted to IMGT
    id_dict = getIDforIMGT(seq_file, start_time, total_count, out_args, id_only) if seq_file else {}
    
    # Create
    imgt_dict = readIMGT(imgt_files)
    writeDb(imgt_dict, no_parse, file_prefix, 'imgt', start_time, total_count, out_args, id_dict=id_dict)
    

def parseIgBlast(seq_file, vgap_file, jaa_file, igblast_output, no_parse, out_args=default_out_args):
    """
    Main for IgBlast aligned sample sequences

    Arguments: 
    seq_file = fasta file input to IgBlast (from which to get sequence)
    vgap_file = file with gapping information for V germlines
    jaa_file = file with conserved amino acid position for J germlines
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
    log['V_GAP_FILE'] = os.path.basename(vgap_file)
    log['J_AA_FILE'] = os.path.basename(jaa_file)
    log['ALIGN_RESULTS'] = os.path.basename(igblast_output)
    log['NO_PARSE'] = no_parse
    printLog(log)
    
    # Get various dictionaries for sequence, gapping, and junction determination
    seq_dict = getSeqforIgBlast(seq_file)
    gap_dict = getGapDict(vgap_file)
    j_dict = getJaaDict(jaa_file)
    
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
    writeDb(igblast_dict, no_parse, file_prefix, 'igblast', start_time, total_count, out_args, seq_dict=seq_dict, gap_dict=gap_dict, j_dict=j_dict)
    

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
                                       help='Aligner used', metavar='')
    
    # Parent parser    
    parser_parent = getCommonArgParser(seq_in=False, seq_out=False, db_out=True, log=False)
    
    # IMGT aligner
    parser_imgt = subparsers.add_parser('imgt', help='Process IMGT/HighV-Quest output', 
                                        parents=[parser_parent], 
                                        formatter_class=ArgumentDefaultsHelpFormatter)
    parser_imgt.set_defaults(func=parseIMGT)
    imgt_arg_group =  parser_imgt.add_mutually_exclusive_group(required=True)
    imgt_arg_group.add_argument('-z', nargs='+', action='store', dest='aligner_output',
                                help='Zipped IMGT output files')
    imgt_arg_group.add_argument('-f', nargs='+', action='store', dest='aligner_output', 
                                help='Folder with unzipped IMGT output files \
                                     (must have 1_Summary, 2_IMGT-gapped, 3_Nt-sequences, and 6_Junction)')
    parser_imgt.add_argument('-s', action='store', nargs='+', dest='in_files',
                             help='List of input FASTA files containing sequences')
    parser_imgt.add_argument('--id', action='store_true', dest='id_only', 
                             help='Specify if only sequence ID passed to IMGT')
    parser_imgt.add_argument('--noparse', action='store_true', dest='no_parse', 
                             help='Specify if input IDs should not be parsed to add new columns to database')
    
    # IgBlast Aligner
    parser_igblast = subparsers.add_parser('igblast', help='Process IgBlast output',
                                           parents=[parser_parent],
                                           formatter_class=ArgumentDefaultsHelpFormatter)
    parser_igblast.set_defaults(func=parseIgBlast)
    parser_igblast.add_argument('-o', nargs='+', action='store', dest='aligner_output', required=True,
                                help='IgBlast output files')
    parser_igblast.add_argument('-s', action='store', nargs='+', dest='in_files', required=True,
                                help='List of input FASTA files containing sequences')
    parser_igblast.add_argument('-j', action='store', nargs='+', dest='jaa_files', required=True,
                                help='List of files with conserved amino acid position for J germlines')
    parser_igblast.add_argument('--vgap', action='store', nargs='+', dest='vgap_files', required=True,
                                help='List of files with gapping information for V germlines')
    parser_igblast.add_argument('--noparse', action='store_true', dest='no_parse', 
                                help='Specify if input IDs should not be parsed to add new columns to database')
    
    return parser


if __name__ == "__main__":
    """
    Parses command line arguments and calls main
    """
    parser = getArgParser()    
    args = parser.parse_args()
    args_dict = parseCommonArgs(args, in_arg='aligner_output')
    
    if 'in_files' in args_dict: del args_dict['in_files']
    else: args_dict['no_parse'] = True
    if 'aligner_output' in args_dict: del args_dict['aligner_output']
    if 'jaa_files' in args_dict: del args_dict['jaa_files']
    if 'vgap_files' in args_dict: del args_dict['vgap_files']
    if 'command' in args_dict: del args_dict['command']
    if 'func' in args_dict: del args_dict['func']           
    
    # IMGT parser
    if args.command == 'imgt':
        if args.__dict__['aligner_output']:
            for i in range(len(args.__dict__['aligner_output'])):
                args_dict['seq_file'] = args.__dict__['in_files'][i] if args.__dict__['in_files'] else None
                args_dict['imgt_output'] = args.__dict__['aligner_output'][i]
                args.func(**args_dict)
        else:
            parser.error('Must include either (-z) zipped IMGT files or \
                         (-f) folder with individual files 1_, 2_, 3_, and 6_')
    elif args.command == 'igblast':
        if args.__dict__['aligner_output']:
            for i in range(len(args.__dict__['aligner_output'])):
                args_dict['seq_file'] = args.__dict__['in_files'][i] if args.__dict__['in_files'] else \
                                        parser.error('Must include fasta file input to IgBlast')
                args_dict['vgap_file'] = args.__dict__['vgap_files'][i] if args.__dict__['vgap_files'] else \
                                         parser.error('Must include file with gapping information for V germlines')
                args_dict['jaa_file'] = args.__dict__['jaa_files'][i] if args.__dict__['jaa_files'] else \
                                        parser.error('Must include fasta file input to IgBlast')
                args_dict['igblast_output'] =  args.__dict__['aligner_output'][i]
        else:
            parser.error('Must include IgBlast output file (-o)')