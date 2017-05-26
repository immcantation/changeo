#!/usr/bin/env python3
"""
Create tab-delimited database file to store sequence alignment information
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
from Bio import SeqIO

# Presto and changeo imports
from presto.Defaults import default_out_args
from presto.Annotation import parseAnnotation
from presto.IO import countSeqFile, printLog, printMessage, printProgress, readSeqFile
from changeo.Commandline import CommonHelpFormatter, checkArgs, getCommonArgParser, parseCommonArgs
from changeo.IO import countDbFile, extractIMGT, getDbWriter, readRepo
from changeo.Parsers import IgBLASTReader, IMGTReader, IHMMuneReader, getIDforIMGT


def getFilePrefix(aligner_output, out_args):
    """
    Get file name prefix and create output directory

    Arguments:
      aligner_output : aligner output file or directory.
      out_args : dictionary of output arguments.

    Returns:
        str : file name prefix.
    """
    # Determine output directory
    if not out_args['out_dir']:
        out_dir = os.path.dirname(os.path.abspath(aligner_output))
    else:
        out_dir = os.path.abspath(out_args['out_dir'])
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

    # Determine file prefix
    if out_args['out_name']:
        file_prefix = out_args['out_name']
    else:
        file_prefix = os.path.splitext(os.path.split(os.path.abspath(aligner_output))[1])[0]

    return os.path.join(out_dir, file_prefix)


def getSeqDict(seq_file):
    """
    Create a dictionary from a sequence file.

    Arguments:
      seq_file : sequence file.

    Returns:
        dict : sequence description as keys with Bio.SeqRecords as values.
    """
    seq_dict = SeqIO.to_dict(readSeqFile(seq_file),
                             key_function=lambda x: x.description)

    return seq_dict


def writeDb(db, fields, file_prefix, total_count, id_dict=None, no_parse=True, partial=False,
            out_args=default_out_args):
    """
    Writes tab-delimited database file in output directory.
    
    Arguments:
      db : a iterator of IgRecord objects containing alignment data.
      fields : a list of ordered field names to write.
      file_prefix : directory and prefix for CLIP tab-delim file.
      total_count : number of records (for progress bar).
      id_dict : a dictionary of the truncated sequence ID mapped to the full sequence ID.
      no_parse : if ID is to be parsed for pRESTO output with default delimiters.
      partial : if True put incomplete alignments in the pass file.
      out_args : common output argument dictionary from parseCommonArgs.

    Returns:
      None
    """
    # Function to check for valid records strictly
    def _pass_strict(rec):
        valid = [rec.v_call and rec.v_call != 'None',
                 rec.j_call and rec.j_call != 'None',
                 rec.functional is not None,
                 rec.seq_vdj,
                 rec.junction]
        return all(valid)

    # Function to check for valid records loosely
    def _pass_gentle(rec):
        valid = [rec.v_call and rec.v_call != 'None',
                 rec.d_call and rec.d_call != 'None',
                 rec.j_call and rec.j_call != 'None']
        return any(valid)

    # Set pass criteria
    _pass = _pass_gentle if partial else _pass_strict

    # Define output file names
    pass_file = '%s_db-pass.tab' % file_prefix
    fail_file = '%s_db-fail.tab' % file_prefix

    # Initiate handles, writers and counters
    pass_handle = None
    fail_handle = None
    pass_writer = None
    fail_writer = None
    start_time = time()
    rec_count = pass_count = fail_count = 0

    # Validate and write output
    printProgress(0, total_count, 0.05, start_time)
    for i, record in enumerate(db, start=1):

        # Replace sequence description with full string, if required
        if id_dict is not None and record.id in id_dict:
            record.id = id_dict[record.id]

        # Parse sequence description into new columns
        if not no_parse:
            try:
                record.annotations = parseAnnotation(record.id, delimiter=out_args['delimiter'])
                record.id = record.annotations['ID']
                del record.annotations['ID']

                # TODO:  This is not the best approach. should pass in output fields.
                # If first record, use parsed description to define extra columns
                if i == 1:  fields.extend(list(record.annotations.keys()))
            except IndexError:
                # Could not parse pRESTO-style annotations so fall back to no parse
                no_parse = True
                sys.stderr.write('\nWARNING: Sequence annotation format not recognized. Sequence headers will not be parsed.\n')

        # Count pass or fail and write to appropriate file
        if _pass(record):
            # Open pass file
            if pass_writer is None:
                pass_handle = open(pass_file, 'wt')
                pass_writer = getDbWriter(pass_handle, add_fields=fields)

            # Write row to pass file
            pass_count += 1
            pass_writer.writerow(record.toDict())
        else:
            # Open failed file
            if out_args['failed'] and fail_writer is None:
                fail_handle = open(fail_file, 'wt')
                fail_writer = getDbWriter(fail_handle, add_fields=fields)

            # Write row to fail file if specified
            fail_count += 1
            if fail_writer is not None:
                fail_writer.writerow(record.toDict())

        # Print progress
        printProgress(i, total_count, 0.05, start_time)

    # Print consol log
    log = OrderedDict()
    log['OUTPUT'] = pass_file
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    log['END'] = 'MakeDb'
    printLog(log)
    
    if pass_handle is not None: pass_handle.close()
    if fail_handle is not None: fail_handle.close()


# TODO:  may be able to merge with other mains
def parseIMGT(aligner_output, seq_file=None, no_parse=True, partial=False,
              parse_scores=False, parse_regions=False, parse_junction=False,
              out_args=default_out_args):
    """
    Main for IMGT aligned sample sequences.

    Arguments:
      aligner_output : zipped file or unzipped folder output by IMGT.
      seq_file : FASTA file input to IMGT (from which to get seqID).
      no_parse : if ID is to be parsed for pRESTO output with default delimiters.
      partial : If True put incomplete alignments in the pass file.
      parse_scores : if True add alignment score fields to output file.
      parse_regions : if True add FWR and CDR region fields to output file.
      out_args : common output argument dictionary from parseCommonArgs.

    Returns:
      None
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'MakeDb'
    log['ALIGNER'] = 'IMGT'
    log['ALIGNER_OUTPUT'] = aligner_output
    log['SEQ_FILE'] = os.path.basename(seq_file) if seq_file else ''
    log['NO_PARSE'] = no_parse
    log['PARTIAL'] = partial
    log['SCORES'] = parse_scores
    log['REGIONS'] = parse_regions
    log['JUNCTION'] = parse_junction
    printLog(log)

    start_time = time()
    printMessage('Loading sequence files', start_time=start_time, width=25)
    # Extract IMGT files
    temp_dir, imgt_files = extractIMGT(aligner_output)
    # Count records in IMGT files
    total_count = countDbFile(imgt_files['summary'])
    # Get (parsed) IDs from fasta file submitted to IMGT
    id_dict = getIDforIMGT(seq_file) if seq_file else {}
    printMessage('Done', start_time=start_time, end=True, width=25)

    # Parse IMGT output and write db
    with open(imgt_files['summary'], 'r') as summary_handle, \
            open(imgt_files['gapped'], 'r') as gapped_handle, \
            open(imgt_files['ntseq'], 'r') as ntseq_handle, \
            open(imgt_files['junction'], 'r') as junction_handle:
        parse_iter = IMGTReader(summary_handle, gapped_handle, ntseq_handle, junction_handle,
                                parse_scores=parse_scores, parse_regions=parse_regions,
                                parse_junction=parse_junction)
        file_prefix = getFilePrefix(aligner_output, out_args)
        writeDb(parse_iter, parse_iter.fields, file_prefix, total_count, id_dict=id_dict,
                no_parse=no_parse, partial=partial, out_args=out_args)

    # Cleanup temp directory
    temp_dir.cleanup()

    return None


# TODO:  may be able to merge with other mains
def parseIgBLAST(aligner_output, seq_file, repo, no_parse=True, partial=False,
                 parse_regions=False, parse_scores=False, parse_igblast_cdr3=False,
                 out_args=default_out_args):
    """
    Main for IgBLAST aligned sample sequences.

    Arguments:
      aligner_output : IgBLAST output file to process.
      seq_file : fasta file input to IgBlast (from which to get sequence).
      repo : folder with germline repertoire files.
      no_parse : if ID is to be parsed for pRESTO output with default delimiters.
      partial : If True put incomplete alignments in the pass file.
      parse_regions : if True add FWR and CDR fields to output file.
      parse_scores : if True add alignment score fields to output file.
      parse_igblast_cdr3 : if True parse CDR3 sequences generated by IgBLAST
      out_args : common output argument dictionary from parseCommonArgs.

    Returns:
      None
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'MakeDB'
    log['ALIGNER'] = 'IgBlast'
    log['ALIGNER_OUTPUT'] = os.path.basename(aligner_output)
    log['SEQ_FILE'] = os.path.basename(seq_file)
    log['NO_PARSE'] = no_parse
    log['PARTIAL'] = partial
    log['SCORES'] = parse_scores
    log['REGIONS'] = parse_regions
    printLog(log)

    start_time = time()
    printMessage('Loading sequence files', start_time=start_time, width=25)
    # Count records in sequence file
    total_count = countSeqFile(seq_file)
    # Get input sequence dictionary
    seq_dict = getSeqDict(seq_file)
    # Create germline repo dictionary
    repo_dict = readRepo(repo)
    printMessage('Done', start_time=start_time, end=True, width=25)

    # Parse and write output
    with open(aligner_output, 'r') as f:
        parse_iter = IgBLASTReader(f, seq_dict, repo_dict,
                                   parse_scores=parse_scores, parse_regions=parse_regions,
                                   parse_igblast_cdr3=parse_igblast_cdr3)
        file_prefix = getFilePrefix(aligner_output, out_args)
        writeDb(parse_iter, parse_iter.fields, file_prefix, total_count,
                no_parse=no_parse, partial=partial, out_args=out_args)

    return None


# TODO:  may be able to merge with other mains
def parseIHMM(aligner_output, seq_file, repo, no_parse=True, partial=False,
              parse_scores=False, parse_regions=False, out_args=default_out_args):
    """
    Main for iHMMuneAlign aligned sample sequences.

    Arguments:
      aligner_output : iHMMune-Align output file to process.
      seq_file : fasta file input to iHMMuneAlign (from which to get sequence).
      repo : folder with germline repertoire files.
      no_parse : if ID is to be parsed for pRESTO output with default delimiters.
      partial : If True put incomplete alignments in the pass file.
      parse_scores : if True parse alignment scores.
      parse_regions : if True add FWR and CDR region fields.
      out_args : common output argument dictionary from parseCommonArgs.

    Returns:
      None
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'MakeDB'
    log['ALIGNER'] = 'iHMMune-Align'
    log['ALIGNER_OUTPUT'] = os.path.basename(aligner_output)
    log['SEQ_FILE'] = os.path.basename(seq_file)
    log['NO_PARSE'] = no_parse
    log['PARTIAL'] = partial
    log['SCORES'] = parse_scores
    log['REGIONS'] = parse_regions
    printLog(log)

    start_time = time()
    printMessage('Loading sequence files', start_time=start_time, width=25)
    # Count records in sequence file
    total_count = countSeqFile(seq_file)
    # Get input sequence dictionary
    seq_dict = getSeqDict(seq_file)
    # Create germline repo dictionary
    repo_dict = readRepo(repo)
    printMessage('Done', start_time=start_time, end=True, width=25)

    # Parse and write output
    with open(aligner_output, 'r') as f:
        parse_iter = IHMMuneReader(f, seq_dict, repo_dict,
                                   parse_scores=parse_scores, parse_regions=parse_regions)
        file_prefix = getFilePrefix(aligner_output, out_args)
        writeDb(parse_iter, parse_iter.fields, file_prefix, total_count,
                no_parse=no_parse, partial=partial, out_args=out_args)

    return None


def getArgParser():
    """
    Defines the ArgumentParser.

    Returns: 
      argparse.ArgumentParser
    """
    fields = dedent(
             '''
              output files:
                  db-pass
                      database of alignment records with functionality information,
                      V and J calls, and a junction region.
                  db-fail
                      database with records that fail due to no functionality information
                      (did not pass IMGT), no V call, no J call, or no junction region.

              universal output fields:
                  SEQUENCE_ID, SEQUENCE_INPUT, SEQUENCE_VDJ, SEQUENCE_IMGT,
                  FUNCTIONAL, IN_FRAME, STOP, MUTATED_INVARIANT, INDELS,
                  V_CALL, D_CALL, J_CALL,
                  V_SEQ_START, V_SEQ_LENGTH,
                  D_SEQ_START, D_SEQ_LENGTH, D_GERM_START, D_GERM_LENGTH,
                  J_SEQ_START, J_SEQ_LENGTH, J_GERM_START, J_GERM_LENGTH,
                  JUNCTION_LENGTH, JUNCTION, NP1_LENGTH, NP2_LENGTH,
                  FWR1_IMGT, FWR2_IMGT, FWR3_IMGT, FWR4_IMGT,
                  CDR1_IMGT, CDR2_IMGT, CDR3_IMGT

              imgt specific output fields:
                  V_GERM_START_IMGT, V_GERM_LENGTH_IMGT,
                  N1_LENGTH, N2_LENGTH, P3V_LENGTH, P5D_LENGTH, P3D_LENGTH, P5J_LENGTH,
                  D_FRAME, V_SCORE, V_IDENTITY, J_SCORE, J_IDENTITY,

              igblast specific output fields:
                  V_GERM_START_VDJ, V_GERM_LENGTH_VDJ,
                  V_EVALUE, V_SCORE, V_IDENTITY, V_BTOP,
                  J_EVALUE, J_SCORE, J_IDENTITY, J_BTOP.
                  CDR3_IGBLAST_NT, CDR3_IGBLAST_AA

              ihmm specific output fields:
                  V_GERM_START_VDJ, V_GERM_LENGTH_VDJ,
                  HMM_SCORE
              ''')
                
    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            formatter_class=CommonHelpFormatter)
    parser.add_argument('--version', action='version',
                        version='%(prog)s:' + ' %s-%s' %(__version__, __date__))
    subparsers = parser.add_subparsers(title='subcommands', dest='command',
                                       help='Aligner used', metavar='')
    # TODO:  This is a temporary fix for Python issue 9253
    subparsers.required = True

    # Parent parser    
    parser_parent = getCommonArgParser(seq_in=False, seq_out=False, log=False)

    # IgBlast Aligner
    parser_igblast = subparsers.add_parser('igblast', parents=[parser_parent],
                                           formatter_class=CommonHelpFormatter,
                                           help='Process IgBLAST output.',
                                           description='Process IgBLAST output.')
    parser_igblast.add_argument('-i', nargs='+', action='store', dest='aligner_outputs',
                                required=True,
                                help='''IgBLAST output files in format 7 with query sequence
                                     (IgBLAST argument \'-outfmt "7 std qseq sseq btop"\').''')
    parser_igblast.add_argument('-r', nargs='+', action='store', dest='repo', required=True,
                                help='''List of folders and/or fasta files containing
                                     IMGT-gapped germline sequences corresponding to the
                                     set of germlines used in the IgBLAST alignment.''')
    parser_igblast.add_argument('-s', action='store', nargs='+', dest='seq_files',
                                required=True,
                                help='''List of input FASTA files (with .fasta, .fna or .fa
                                     extension), containing sequences.''')
    parser_igblast.add_argument('--noparse', action='store_true', dest='no_parse',
                                help='''Specify to prevent input sequence headers from being parsed
                                    to add new columns to database. Parsing of sequence headers requires
                                    headers to be in the pRESTO annotation format, so this should be specified
                                    when sequence headers are incompatible with the pRESTO annotation scheme.
                                    Note, unrecognized header formats will default to this behavior.''')
    parser_igblast.add_argument('--partial', action='store_true', dest='partial',
                                help='''If specified, include incomplete V(D)J alignments in
                                     the pass file instead of the fail file.''')
    parser_igblast.add_argument('--scores', action='store_true', dest='parse_scores',
                                help='''Specify if alignment score metrics should be
                                     included in the output. Adds the V_SCORE, V_IDENTITY,
                                     V_EVALUE, V_BTOP, J_SCORE, J_IDENTITY,
                                     J_BTOP, and J_EVALUE columns.''')
    parser_igblast.add_argument('--regions', action='store_true', dest='parse_regions',
                                help='''Specify if IMGT FWR and CDRs should be
                                     included in the output. Adds the FWR1_IMGT, FWR2_IMGT,
                                     FWR3_IMGT, FWR4_IMGT, CDR1_IMGT, CDR2_IMGT, and
                                     CDR3_IMGT columns.''')
    parser_igblast.add_argument('--cdr3', action='store_true',
                                dest='parse_igblast_cdr3', 
                                help='''Specify if the CDR3 sequences generated by IgBLAST 
                                     should be included in the output. Adds the columns
                                     CDR3_IGBLAST_NT and CDR3_IGBLAST_AA. Requires IgBLAST
                                     version 1.5 or greater.''')
    parser_igblast.set_defaults(func=parseIgBLAST)

    # IMGT aligner
    parser_imgt = subparsers.add_parser('imgt', parents=[parser_parent],
                                        formatter_class=CommonHelpFormatter,
                                        help='''Process IMGT/HighV-Quest output
                                             (does not work with V-QUEST).''',
                                        description='''Process IMGT/HighV-Quest output
                                             (does not work with V-QUEST).''')
    parser_imgt.add_argument('-i', nargs='+', action='store', dest='aligner_outputs',
                             help='''Either zipped IMGT output files (.zip or .txz) or a
                                  folder containing unzipped IMGT output files (which must
                                  include 1_Summary, 2_IMGT-gapped, 3_Nt-sequences,
                                  and 6_Junction).''')
    parser_imgt.add_argument('-s', nargs='*', action='store', dest='seq_files',
                             required=False,
                             help='''List of input FASTA files (with .fasta, .fna or .fa
                                  extension) containing sequences.''')
    parser_imgt.add_argument('--noparse', action='store_true', dest='no_parse', 
                             help='''Specify to prevent input sequence headers from being parsed
                                  to add new columns to database. Parsing of sequence headers requires
                                  headers to be in the pRESTO annotation format, so this should be specified
                                  when sequence headers are incompatible with the pRESTO annotation scheme.
                                  Note, unrecognized header formats will default to this behavior.''')
    parser_imgt.add_argument('--partial', action='store_true', dest='partial',
                             help='''If specified, include incomplete V(D)J alignments in
                                  the pass file instead of the fail file.''')
    parser_imgt.add_argument('--scores', action='store_true', dest='parse_scores',
                             help='''Specify if alignment score metrics should be
                                  included in the output. Adds the V_SCORE, V_IDENTITY,
                                  J_SCORE and J_IDENTITY.''')
    parser_imgt.add_argument('--regions', action='store_true', dest='parse_regions',
                             help='''Specify if IMGT FWRs and CDRs should be
                                  included in the output. Adds the FWR1_IMGT, FWR2_IMGT,
                                  FWR3_IMGT, FWR4_IMGT, CDR1_IMGT, CDR2_IMGT, and
                                  CDR3_IMGT columns.''')
    parser_imgt.add_argument('--junction', action='store_true', dest='parse_junction',
                             help='''Specify if detailed junction fields should be
                                  included in the output. Adds the columns 
                                  N1_LENGTH, N2_LENGTH, P3V_LENGTH, P5D_LENGTH, P3D_LENGTH,
                                  P5J_LENGTH, D_FRAME.''')
    parser_imgt.set_defaults(func=parseIMGT)

    # iHMMuneAlign Aligner
    parser_ihmm = subparsers.add_parser('ihmm', parents=[parser_parent],
                                        formatter_class=CommonHelpFormatter,
                                        help='Process iHMMune-Align output.',
                                        description='Process iHMMune-Align output.')
    parser_ihmm.add_argument('-i', nargs='+', action='store', dest='aligner_outputs',
                             required=True,
                             help='''iHMMune-Align output file.''')
    parser_ihmm.add_argument('-r', nargs='+', action='store', dest='repo', required=True,
                             help='''List of folders and/or FASTA files containing
                                  IMGT-gapped germline sequences corresponding to the
                                  set of germlines used in the IgBLAST alignment.''')
    parser_ihmm.add_argument('-s', action='store', nargs='+', dest='seq_files',
                             required=True,
                             help='''List of input FASTA files (with .fasta, .fna or .fa
                                  extension) containing sequences.''')
    parser_ihmm.add_argument('--noparse', action='store_true', dest='no_parse',
                             help='''Specify to prevent input sequence headers from being parsed
                                  to add new columns to database. Parsing of sequence headers requires
                                  headers to be in the pRESTO annotation format, so this should be specified
                                  when sequence headers are incompatible with the pRESTO annotation scheme.
                                  Note, unrecognized header formats will default to this behavior.''')
    parser_ihmm.add_argument('--partial', action='store_true', dest='partial',
                             help='''If specified, include incomplete V(D)J alignments in
                                  the pass file instead of the fail file.''')
    parser_ihmm.add_argument('--scores', action='store_true', dest='parse_scores',
                             help='''Specify if alignment score metrics should be
                                  included in the output. Adds the path score of the
                                  iHMMune-Align hidden Markov model to HMM_SCORE.''')
    parser_ihmm.add_argument('--regions', action='store_true', dest='parse_regions',
                             help='''Specify if IMGT FWRs and CDRs should be
                                  included in the output. Adds the FWR1_IMGT, FWR2_IMGT,
                                  FWR3_IMGT, FWR4_IMGT, CDR1_IMGT, CDR2_IMGT, and
                                  CDR3_IMGT columns.''')
    parser_ihmm.set_defaults(func=parseIHMM)

    return parser
    
    
if __name__ == "__main__":
    """
    Parses command line arguments and calls main
    """
    parser = getArgParser()
    checkArgs(parser)
    args = parser.parse_args()
    args_dict = parseCommonArgs(args, in_arg='aligner_outputs')

    # Set no ID parsing if sequence files are not provided
    if 'seq_files' in args_dict and not args_dict['seq_files']:
        args_dict['no_parse'] = True

    # Delete
    if 'seq_files' in args_dict: del args_dict['seq_files']
    if 'aligner_outputs' in args_dict: del args_dict['aligner_outputs']
    if 'command' in args_dict: del args_dict['command']
    if 'func' in args_dict: del args_dict['func']           
    
    if args.command == 'imgt':
        for i in range(len(args.__dict__['aligner_outputs'])):
            args_dict['aligner_output'] = args.__dict__['aligner_outputs'][i]
            args_dict['seq_file'] = args.__dict__['seq_files'][i] \
                                    if args.__dict__['seq_files'] else None
            args.func(**args_dict)
    elif args.command == 'igblast' or args.command == 'ihmm':
        for i in range(len(args.__dict__['aligner_outputs'])):
            args_dict['aligner_output'] =  args.__dict__['aligner_outputs'][i]
            args_dict['seq_file'] = args.__dict__['seq_files'][i]
            args.func(**args_dict)
