#!/usr/bin/env python
"""
Create tab-delimited database file to store sequence alignment information
"""

__author__    = 'Namita Gupta'
__copyright__ = 'Copyright 2014 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.0'
__date__      = '2014.4.10'

# Imports
import re
import csv
import os, sys
from zipfile import ZipFile
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from itertools import izip
from collections import OrderedDict
from time import time

# IgCore and DbCore imports 
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from IgCore import default_out_args, parseAnnotation, printLog, printProgress
from IgCore import getCommonArgParser, parseCommonArgs
from DbCore import getDbWriter, IgRecord, countDbFile

# Default parameters
default_V_regex = re.compile(r"(IG[HLK][V]\d+[-/\w]*[-\*][\.\w]+)")
default_D_regex = re.compile(r"(IG[HLK][D]\d+[-/\w]*[-\*][\.\w]+)")
default_J_regex = re.compile(r"(IG[HLK][J]\d+[-/\w]*[-\*][\.\w]+)")


def extractIMGT(imgt_zipfile):
    """
    Extract necessary files from IMGT zipped results
    
    Arguments:
    imgt_zipfile = zipped file output by IMGT
    
    Returns:
    sorted list of filenames from which information will be read
    """
    # Extract selected files from the IMGT zip file
    # Try to retain compressed format
    imgt_zip = ZipFile(imgt_zipfile, 'r')
    imgt_flags = ["1_Summary", "2_IMGT-gapped", "3_Nt-sequences", "6_Junction"]
    imgt_files = sorted([n for n in imgt_zip.namelist() for i in imgt_flags if i in n])
    for f in imgt_files:  imgt_zip.extract(f)
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
               'D_SEQ_START':        sum(int(i) for i in [len(nt['V-REGION']),
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
               'J_SEQ_START':        sum(int(i) for i in [len(nt['V-REGION']),
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


def writeCLIP(db_gen, parse_id, file_prefix, aligner, start_time, total_count, out_args, id_dict={}, seq_dict={}):
    """
    Writes CLIP intermediate tab-delim file in current directory
    
    Arguments:
    db_gen = a generator of IgRecord objects containing alignment data
    id_dict = a dictionary of {truncated ID: full seq description}
    file_prefix = directory and prefix for CLIP tab-delim file
    
    Returns:
    None
    """
    pass_file = "%s_CLIP-pass.tab" % file_prefix
    fail_file = "%s_CLIP-fail.tab" % file_prefix
    if aligner=='imgt':
        ordered_fields = ['SEQUENCE_ID','SEQUENCE','FUNCTIONAL','IN_FRAME','STOP','MUTATED_INVARIANT','INDELS',
                          'V_MATCH','V_LENGTH','J_MATCH','J_LENGTH','V_CALL','D_CALL','J_CALL','SEQUENCE_GAP',
                          'V_SEQ_START','V_SEQ_LENGTH','V_GERM_START','V_GERM_LENGTH','N1_LENGTH','D_SEQ_START',
                          'D_SEQ_LENGTH','D_GERM_START','D_GERM_LENGTH','N2_LENGTH','J_SEQ_START','J_SEQ_LENGTH',
                          'J_GERM_START','J_GERM_LENGTH','JUNCTION_GAP_LENGTH','JUNCTION']
   
    pass_handle = open(pass_file, 'wb')
    fail_handle = open(fail_file, 'wb')
    # Create DbWriter
    if parse_id:
        for v in id_dict.itervalues():
            tmp = parseAnnotation(v, delimiter=out_args['delimiter'])
            del tmp['ID']
            ordered_fields.extend(tmp.keys())
            break
    pass_writer = getDbWriter(pass_handle, add_fields=ordered_fields)
    fail_writer = getDbWriter(fail_handle, add_fields=['SEQUENCE_ID','SEQUENCE'])
    # Initialize counters
    pass_count = fail_count = 0
    
    for i,record in enumerate(db_gen):
        printProgress(i+1 + (total_count/2 if parse_id else 0), total_count, 0.05, start_time)
        # Count pass or fail
        if record.v_call == 'None' and record.j_call == 'None': 
            fail_count += 1
            fail_writer.writerow(record.toDict())
            continue
        else: 
            pass_count += 1
        # Build sample sequence description
        if record.id.split(' ')[0] in id_dict:
            record.id = id_dict[record.id]
        # Parse sequence description into new columns
        if parse_id:
            record.annotations = parseAnnotation(record.id, delimiter=out_args['delimiter'])
            record.id = record.annotations['ID']
            del record.annotations['ID']
        # Write row to tab-delim CLIP file
        pass_writer.writerow(record.toDict())
    
    # Print log
    printProgress(i+1 + (total_count/2 if parse_id else 0), total_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT'] = pass_file
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    log['END'] = 'MakeClip'
    printLog(log)
    

def parseIMGT(seq_file, zip_file, imgt_files, id_only, parse_id, out_args=default_out_args):
    """
    Main for IMGT aligned sample sequences

    Arguments: 
    seq_file = FASTA file input to IMGT (from which to get seqID)
    zip_file = zipped IMGT output file to process
    imgt_files = list of 1_Summary, 2_IMGT-gapped, 3_Nt-sequences, 6_Junction filenames
        
    Returns: 
    None
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'MakeClip'
    log['ALIGNER'] = 'IMGT'
    log['SEQ_FILE'] = os.path.basename(seq_file)
    log['ALIGN_RESULTS'] = zip_file if zip_file is not None else \
    '_'.join( filter( None, os.path.basename(imgt_files[0]).split('_') )[2:-1] )
    log['ID_ONLY'] = id_only 
    log['PARSE_ID'] = parse_id
    printLog(log)
    
    # Unzip zipped file
    if zip_file:
        imgt_files = extractIMGT(zip_file)
    file_prefix = '_'.join( filter( None, os.path.basename(imgt_files[0]).split('_') )[2:-1] )
        
    # Formalize out_dir and file-prefix
    if not out_args['out_dir']:
        out_dir = os.path.split(seq_file)[0]
    else:
        out_dir = os.path.abspath(out_args['out_dir'])
        if not os.path.exists(out_dir):  os.mkdir(out_dir)
    file_prefix = os.path.join(out_dir, file_prefix)
    
    total_count = countDbFile(imgt_files[0]) * (2 if parse_id else 1)
    start_time = time()
    
    # Get (parsed) IDs from fasta file submitted to IMGT
    id_dict = getIDforIMGT(seq_file, start_time, total_count, out_args, id_only) if seq_file else {}
    
    # Create
    imgt_dict = readIMGT(imgt_files)
    writeCLIP(imgt_dict, parse_id, file_prefix, 'imgt', start_time, total_count, out_args, id_dict)
    

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
    parser_parent = getCommonArgParser(seq_in=False, seq_out=False, log=False)
    
    # IMGT aligner
    parser_imgt = subparsers.add_parser('imgt', help='Process IMGT/HighV-Quest output', 
                                        parents=[parser_parent], 
                                        formatter_class=ArgumentDefaultsHelpFormatter)
    parser_imgt.set_defaults(func=parseIMGT)
    imgt_arg_group =  parser_imgt.add_mutually_exclusive_group(required=True)
    imgt_arg_group.add_argument('-z', nargs='+', action='store', dest='zip_files',
                                help='Zipped IMGT output files')
    imgt_arg_group.add_argument('-f', nargs='+', action='store', dest='al_folders', 
                                help='Folder with unzipped IMGT files \
                                     (must have 1_Summary, 2_IMGT-gapped, 3_Nt-sequences, and 6_Junction)')
    parser_imgt.add_argument('-s', action='store', nargs='+', dest='seq_files',
                             help='List of input FASTA files containing sequences')
    parser_imgt.add_argument('--id', action='store_true', dest='id_only', default=False,
                             help='Specify if only sequence ID passed to IMGT')
    parser_imgt.add_argument('--noParse', action='store_false', dest='parse_id', default=True,
                             help='Specify if input IDs should not be parsed to add new columns to CLIP-tab')
    
    return parser


if __name__ == "__main__":
    """
    Parses command line arguments and calls main
    """
    parser = getArgParser()    
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)
    
    if 'seq_files' in args_dict: del args_dict['seq_files']
    if 'zip_files' in args_dict: del args_dict['zip_files']
    if 'al_folders' in args_dict: del args_dict['al_folders']
    if 'command' in args_dict: del args_dict['command']
    if 'func' in args_dict: del args_dict['func']           
    
    # IMGT parser
    if args.command == 'imgt':
        if args.__dict__['zip_files']: # input IMGT zip-files
            for i in range(len(args.__dict__['zip_files'])):
                args_dict['zip_file'] = args.__dict__['zip_files'][i]
                args_dict['seq_file'] = args.__dict__['seq_files'][i] if args.__dict__['seq_files'] else None
                args_dict['imgt_files'] = None
                args.func(**args_dict)
        elif args.__dict__['al_folders']: # input folders with IMGT summary files
            imgt_flags = ["1_Summary", "2_IMGT-gapped", "3_Nt-sequences", "6_Junction"] # necessary files
            for i in range( len(args.__dict__['al_folders']) ):
                folder = args.__dict__['al_folders'][i]
                imgt_files = sorted([ (folder + ('/' if folder[-1]!='/' else '') + n) \
                                    for n in os.listdir(folder) for j in imgt_flags if j in n ])
                if all( j in f for j,f in zip(imgt_flags, imgt_files) ):
                    args_dict['seq_file'] = args.__dict__['seq_files'][i] if args.__dict__['seq_files'] else None
                    args_dict['zip_file'] = None
                    args_dict['imgt_files'] = imgt_files
                    args.func(**args_dict)
                elif len(imgt_files) >= len(imgt_flags): # e.g. multiple 1_Summary files
                    parser.error('Wrong files in folder %s' % folder)
                else:
                    parser.error('Missing necessary file in folder %s' % folder)
        else:
            parser.error('Must include either (-z) zipped IMGT files or \
                         (-f) folder with 1_, 2_, 3_, and 6_ individual files')