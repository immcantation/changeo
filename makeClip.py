#!/usr/bin/env python
##/usr/bin/env /usr/local/cluster/software/installation/python/Python-2.7/bin/python
"""
Create CLIP tab-delimited file to store sequence alignment information
"""

__author__    = 'Namita Gupta'
__copyright__ = 'Copyright 2013, Kleinstein Lab, Yale University School of Medicine'
__license__   = 'GPLv3'
__version__   = '0.3'
__date__      = '2014.01.24'

# Imports
import re
import csv
from os import listdir, mkdir, path
from zipfile import ZipFile
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from itertools import izip

# ChAnGEo imports
from sys import path as syspath
syspath.append(path.dirname(path.realpath(__file__)))
from IgCore import default_out_args, parseAnnotation
from IgCore import getCommonArgParser, parseCommonArgs
from IgCore import getFileType, getOutputHandle, printLog, printProgress
from DbCore import countDbFile, readDbFile, getDbWriter, IgRecord

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
    db_flags = ["1_Summary", "2_IMGT-gapped", "3_Nt-sequences", "6_Junction"]
    db_files = sorted([n for n in imgt_zip.namelist() for d in db_flags if d in n])
    for f in db_files:  imgt_zip.extract(f)
    return db_files


def readIMGT(db_files):
    """
    Reads IMGT/HighV-Quest output

    Arguments: 
    db_files = IMGT/HighV-Quest output files 1, 2, 3, and 6
        
    Returns: 
    a generator of dictionaries containing alignment data
    """
    db_iters = [csv.DictReader(open(d, 'rU'), delimiter='\t') for d in db_files]
    # Create a generator of dictionaries for each sequence alignment
    db_gen = ({'SEQUENCE_ID':              sm['Sequence ID'],
               'SEQUENCE':                sm['Sequence'],
               'FUNCTIONAL':      ['?','T','F'][('productive' in sm['Functionality']) + ('unprod' in sm['Functionality'])],
               'IN_FRAME':           ['?','T','F'][('in-frame' in sm['JUNCTION frame']) + ('out-of-frame' in sm['JUNCTION frame'])],
               'STOP':         ['F','?','T'][('stop codon' in sm['Functionality comment']) + ('unprod' in sm['Functionality'])],
               'MUTATED_INVARIANT':      ['F','?','T'][(any(('missing' in sm['Functionality comment'],
                                                       'missing' in sm['V-REGION potential ins/del']))) + ('unprod' in sm['Functionality'])],
               'INDELS':             ['F','T'][any((sm['V-REGION potential ins/del'], 
                                                    sm['V-REGION insertions'], 
                                                    sm['V-REGION deletions']))],
               'V_MATCH':            0 if sm['V-REGION identity nt'] == 'null' or not sm['V-REGION identity nt'] \
                                     else int(sm['V-REGION identity nt'].split('/')[0] or 0) + nt['V-REGION'].count('n'),
               'V_LENGTH':               0 if sm['V-REGION identity nt'] == 'null' or not sm['V-REGION identity nt'] \
                                     else int(sm['V-REGION identity nt'].split('/')[1].split()[0]),
               'J_MATCH':            0 if sm['J-REGION identity nt'] == 'null' or not sm['J-REGION identity nt'] \
                                     else int(sm['J-REGION identity nt'].split('/')[0] or 0) + nt['J-REGION'].count('n'),
               'J_LENGTH':               0 if sm['J-REGION identity nt'] == 'null' or not sm['J-REGION identity nt'] \
                                     else int(sm['J-REGION identity nt'].split('/')[1].split()[0]),
               'V_CALL':             re.sub( '\sor\s', ',', re.sub(',','',gp['V-GENE and allele']) ), # replace or with comma
               'D_CALL':             re.sub( '\sor\s', ',', re.sub(',','',gp['D-GENE and allele']) ),
               'J_CALL':             re.sub( '\sor\s', ',', re.sub(',','',gp['J-GENE and allele']) ),
               'SEQUENCE_GAP':            gp['V-D-J-REGION'] if gp['V-D-J-REGION'] else gp['V-J-REGION'],
               'V_GAP_LENGTH':               len(gp['V-REGION']) if gp['V-REGION'] else 0,
               'N1_LENGTH':              sum(int(i) for i in [jn["P3'V-nt nb"], 
                                                          jn['N-REGION-nt nb'], 
                                                          jn['N1-REGION-nt nb'], 
                                                          jn["P5'D-nt nb"]] if i),
               'D_5_TRIM':         int(jn["5'D-REGION trimmed-nt nb"] or 0),
               'D_LENGTH':        int(jn["D-REGION-nt nb"] or 0),
               'N2_LENGTH':              sum(int(i) for i in [jn["P3'D-nt nb"],
                                                          jn['N2-REGION-nt nb'],
                                                          jn["P5'J-nt nb"]] if i),   
               'J_5_TRIM':         int(jn["5'J-REGION trimmed-nt nb"] or 0),
               'J_GAP_LENGTH':               len(gp['J-REGION']) if gp['J-REGION'] else 0,
               'JUNCTION_GAP_LENGTH':        len(jn['JUNCTION']) if jn['JUNCTION'] else 0,
               'JUNCTION':           jn['JUNCTION']} if "No results" not in sm['Functionality'] else 
                                                                             { 'SEQUENCE_ID':sm['Sequence ID'],'FUNCTIONAL':'No Results' }
              for sm, gp, nt, jn in izip(*db_iters) )
    
    return db_gen

    
def getIDforIMGT(seq_file, id_only=False):
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
    for seq in seq_dict.itervalues():
        if id_only:
            id_key = parseAnnotation(seq.description, fields=['ID'])['ID']
        else:
            id_key = re.sub('\||\s','_',seq.description[:50])
        ids.update({id_key:seq.description})
    return ids


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
    for line in handle:
        if(re.match(query, line)):
            return line


def readIgBlast(igblast_file, seq_dict):
    """
    Reads IMGT/HighV-Quest output

    Arguments: 
    igblast_file = IgBlast output file (format 7)
        
    Returns: 
    a generator of dictionaries containing alignment data
    """
    db_gen = {}
    igblast_handle = open(igblast_file, 'rU')
    
    for line in igblast_handle:
        if(re.match("# Query:", line)):
            if 'SEQUENCE_ID' in db_gen and db_gen['FUNCTIONAL'] != 'No Results': yield db_gen
            words = line.split()
            db_gen = {}
            db_gen['SEQUENCE_ID'] = ' '.join(words[2:])
            db_gen['SEQUENCE'] = seq_dict[db_gen['SEQUENCE_ID']]
        if(re.match("# 0 hits found", line)):
            db_gen['FUNCTIONAL'] = 'No Results'
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
                vs = len(db_gen['V_CALL'].split(','))
                for i in range(vs):
                    line = findLine(igblast_handle,"V")
                    words = line.split()
                    db_gen['V_LENGTH'] = ','.join([db_gen.get('V_LENGTH',''),words[4]]) if('V_LENGTH' in db_gen) else words[4]
                    db_gen['V_MATCH'] = ','.join([db_gen.get('V_MATCH',''),str(int(words[4]) - int(words[5]))]) if('V_MATCH' in db_gen) else str(int(words[4]) - int(words[5]))
                    db_gen['V_START_SEQ'] = ','.join([db_gen.get('V_START_SEQ',''),words[8]]) if('V_START_SEQ' in db_gen) else words[8]
                    db_gen['V_END_SEQ'] = ','.join([db_gen.get('V_END_SEQ',''),words[9]]) if('V_END_SEQ' in db_gen) else words[9]
                    db_gen['V_START_GERM'] = ','.join([db_gen.get('V_START_GERM',''),words[10]]) if('V_START_GERM' in db_gen) else words[10]
                    db_gen['V_END_GERM'] = ','.join([db_gen.get('V_END_GERM',''),words[11]]) if('V_END_GERM' in db_gen) else words[11]
            if('D_CALL' in db_gen and db_gen['D_CALL']):
                ds = len(db_gen['D_CALL'].split(','))
                for i in range(ds):
                    line = findLine(igblast_handle,"D")
                    words = line.split()
                    db_gen['D_LENGTH'] = ','.join([db_gen.get('D_LENGTH',''),words[4]]) if 'D_LENGTH' in db_gen else words[4]
                    db_gen['D_MATCH'] = ','.join([db_gen.get('D_MATCH',''),str(int(words[4]) - int(words[5]))]) if 'D_MATCH' in db_gen else str(int(words[4]) - int(words[5]))
                    db_gen['D_START_SEQ'] = ','.join([db_gen.get('D_START_SEQ',''),words[8]]) if 'D_START_SEQ' in db_gen else words[8]
                    db_gen['D_END_SEQ'] = ','.join([db_gen.get('D_END_SEQ',''),words[9]]) if 'D_END_SEQ' in db_gen else words[9]
                    db_gen['D_START_GERM'] = ','.join([db_gen.get('D_START_GERM',''),words[10]]) if 'D_START_GERM' in db_gen else words[10]
                    db_gen['D_END_GERM'] = ','.join([db_gen.get('D_END_GERM',''),words[11]]) if 'D_END_GERM' in db_gen else words[11]
            if('J_CALL' in db_gen and db_gen['J_CALL']):
                js = len(db_gen['J_CALL'].split(','))
                for i in range(js):
                    line = findLine(igblast_handle,"J")
                    words = line.split()
                    db_gen['J_LENGTH'] = ','.join([db_gen.get('J_LENGTH',''),words[4]]) if 'J_LENGTH' in db_gen else words[4]
                    db_gen['J_MATCH'] = ','.join([db_gen.get('J_MATCH',''),str(int(words[4]) - int(words[5]))]) if 'J_MATCH' in db_gen else str(int(words[4]) - int(words[5]))
                    db_gen['J_START_SEQ'] = ','.join([db_gen.get('J_START_SEQ',''),words[8]]) if 'J_START_SEQ' in db_gen else words[8]
                    db_gen['J_END_SEQ'] = ','.join([db_gen.get('J_END_SEQ',''),words[9]]) if 'J_END_SEQ' in db_gen else words[9]
                    db_gen['J_START_GERM'] = ','.join([db_gen.get('J_START_GERM',''),words[10]]) if 'J_START_GERM' in db_gen else words[10]
                    db_gen['J_END_GERM'] = ','.join([db_gen.get('J_END_GERM',''),words[11]]) if 'J_END_GERM' in db_gen else words[11]
    igblast_handle.close()
    yield db_gen
        

def writeCLIP(align_dict, parse_id, file_prefix, aligner, id_dict={}, seq_dict={}):
    """
    Writes CLIP intermediate tab-delim file in current directory
    
    Arguments:
    align_dict = a generator of dictionaries containing alignment data
    id_dict = a dictionary of {truncated ID: full seq description}
    file_prefix = directory and prefix for CLIP tab-delim file
    
    Returns:
    None
    """
    clip_file = "%s_CLIP.tab" % file_prefix
    if aligner=='imgt':
        ordered_fields = ['SEQUENCE_ID','SEQUENCE','FUNCTIONAL','IN_FRAME','STOP','MUTATED_INVARIANT','INDELS',
                          'V_MATCH','V_LENGTH','J_MATCH','J_LENGTH','V_CALL','D_CALL','J_CALL','SEQUENCE_GAP',
                          'V_GAP_LENGTH','N1_LENGTH','D_5_TRIM','D_LENGTH','N2_LENGTH','J_5_TRIM','J_GAP_LENGTH',
                          'JUNCTION_GAP_LENGTH','JUNCTION']
    elif aligner=='igblast':
        ordered_fields = ['SEQUENCE_ID','SEQUENCE','FUNCTIONAL','IN_FRAME','STOP','V_MATCH','V_LENGTH','J_MATCH','J_LENGTH',
                          'V_START_SEQ','V_END_SEQ','V_START_GERM','V_END_GERM','D_START_SEQ','D_END_SEQ',
                          'D_START_GERM','D_END_GERM','J_START_SEQ','J_END_SEQ','J_START_GERM','J_END_GERM',
                          'V_CALL','D_CALL','J_CALL','SEQUENCE_GAP','V_GAP_LENGTH','JUNCTION_GAP_LENGTH','JUNCTION',
                          'MUTATED_INVARIANT','INDELS','N1_LENGTH','D_5_TRIM','D_LENGTH','N2_LENGTH','J_5_TRIM','J_GAP_LENGTH',
                          'SUBJECT_ID','TIME','COMPARTMENT']
    with open(clip_file, 'wb') as clip_handle:
        
        for i,align in enumerate(align_dict):
            # Build sample sequence description
            if align.get('SEQUENCE_ID','').split(' ')[0] in id_dict:
                align['SEQUENCE_ID'] = "%s" % id_dict[align['SEQUENCE_ID']]
            # print ">%s" % align['SEQUENCE_ID']
            
            if parse_id:
                # Parse sequence description into new columns
                id_field_dict = parseAnnotation(align['SEQUENCE_ID'])
                align['SEQUENCE_ID'] = id_field_dict['ID']
                del id_field_dict['ID']
                if i==0:
                    ordered_fields.extend(id_field_dict.keys())
                    dict_writer = csv.DictWriter(clip_handle, delimiter='\t', fieldnames=ordered_fields, extrasaction='ignore')
                    dict_writer.writeheader()
                # Add description field entries to align dictionary
                for k,v in id_field_dict.iteritems(): align[k] = v
            elif i==0: 
                    dict_writer = csv.DictWriter(clip_handle, delimiter='\t', fieldnames=ordered_fields, extrasaction='ignore')
                    dict_writer.writeheader()
            
            # Write row to tab-delim CLIP file
            dict_writer.writerow(align)
    print "%d records written to CLIP-tab file" % (i+1)
    

def parseIMGT(seq_file, zip_file, db_files, id_only, parse_id, out_args=default_out_args):
    """
    Main for IMGT aligned sample sequences

    Arguments: 
    seq_file = FASTA file input to IMGT (from which to get seqID)
    zip_file = zipped IMGT output file to process
    db_files = list of 1_Summary, 2_IMGT-gapped, 3_Nt-sequences, 6_Junction filenames
        
    Returns: 
    None
    """
    print "-> parseIMGT" 
    if zip_file: 
        print "Unzipping..."
        db_files = extractIMGT(zip_file)
    print "Reading IMGT files..."
    imgt_dict = readIMGT(db_files)
    # Formalize out_dir and file-prefix
    if not out_args['out_dir']:
        out_dir = path.split(db_files[0])[0]
    else:
        out_dir = path.abspath(out_args['out_dir'])
        if not path.exists(out_dir):  mkdir(out_dir)    
    file_prefix = path.join( out_dir, '_'.join( filter( None, path.basename(db_files[0]).split('_') )[2:-1] ) )
    id_dict = getIDforIMGT(seq_file, id_only) if seq_file else {}
    print "Writing CLIP-tab file..."
    writeCLIP(imgt_dict, parse_id, file_prefix, 'imgt', id_dict)
    
    
def parseIgBlast(seq_file, igblast_file, parse_id, out_args=default_out_args):
    """
    Main for IgBlast aligned sample sequences

    Arguments: 
    seq_file = FASTA file input to IgBlast (from which to get seq)
    igblast_file = IgBlast output file to process
        
    Returns: 
    None
    """
    print "-> parseIgBlast" 
    print "Getting sequences from input file"
    seq_dict = getSeqforIgBlast(seq_file)
    print "Reading IgBlast output file..."
    igblast_dict = readIgBlast(igblast_file, seq_dict)
    # Formalize out_dir and file-prefix
    if not out_args['out_dir']:
        out_dir = path.split(igblast_file)[0]
    else:
        out_dir = path.abspath(out_args['out_dir'])
        if not path.exists(out_dir):  mkdir(out_dir)
    file_prefix = path.join( out_dir, path.basename(path.splitext(igblast_file)[0]) )
    print "Writing CLIP-tab file..."
    writeCLIP(igblast_dict, parse_id, file_prefix, 'igblast')
       

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
    parser_parent = getCommonArgParser(seq_in=False, seq_out=False, db_in=True, multiproc=False)
    
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
    
    # IgBlast Aligner
    parser_igblast = subparsers.add_parser('igblast', help='Process IgBlast output',
                                           parents=[parser_parent],
                                           formatter_class=ArgumentDefaultsHelpFormatter)
    parser_igblast.set_defaults(func=parseIgBlast)
    parser_igblast.add_argument('-o', nargs='+', action='store', dest='igblast_files', required=True,
                                help='IgBlast output files')
    parser_igblast.add_argument('-s', action='store', nargs='+', dest='seq_files',
                             help='List of input FASTA files containing sequences')
    parser_igblast.add_argument('--noParse', action='store_false', dest='parse_id', default=True,
                             help='Specify if input IDs should not be parsed to add new columns to CLIP-tab')


if __name__ == "__main__":
    """
    Parses command line arguments and calls main
    """
    parser = getArgParser()    
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)
    
    if 'seq_files' in args_dict: del args_dict['seq_files']
    if 'igblast_files' in args_dict: del args_dict['igblast_files']
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
                args_dict['db_files'] = None
                args.func(**args_dict)
                print ''
        elif args.__dict__['al_folders']: # input folders with IMGT summary files
            db_flags = ["1_Summary", "2_IMGT-gapped", "3_Nt-sequences", "6_Junction"] # necessary files
            for i in range( len(args.__dict__['al_folders']) ):
                folder = args.__dict__['al_folders'][i]
                db_files = sorted([ (folder + ('/' if folder[-1]!='/' else '') + n) \
                                    for n in listdir(folder) for d in db_flags if d in n ])
                if all( d in f for d,f in zip(db_flags, db_files) ):
                    args_dict['seq_file'] = args.__dict__['seq_files'][i] if args.__dict__['seq_files'] else None
                    args_dict['zip_file'] = None
                    args_dict['db_files'] = db_files
                    args.func(**args_dict)
                elif len(db_files) >= len(db_flags): # e.g. multiple 1_Summary files
                    parser.error('Wrong files in folder %s' % folder)
                else:
                    parser.error('Missing necessary file in folder %s' % folder)
        else:
            parser.error('Must include either (-z) zipped IMGT files or \
                         (-f) folder with 1_, 2_, 3_, and 6_ individual files')
    # IgBlast parser
    elif args.command == 'igblast':
        for i in range(len(args.__dict__['igblast_files'])):
            args_dict['igblast_file'] = args.__dict__['igblast_files'][i]
            args_dict['seq_file'] = args.__dict__['seq_files'][i] if args.__dict__['seq_files'] else None
            args.func(**args_dict)
            print ''
