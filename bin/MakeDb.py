#!/usr/bin/env python3
"""
Create tab-delimited database file to store sequence alignment information
"""
# Info
__author__ = 'Namita Gupta, Jason Anthony Vander Heiden'
from changeo import __version__, __date__

# Imports
import csv
import os
import re
from argparse import ArgumentParser
from collections import OrderedDict
from shutil import rmtree
from textwrap import dedent
from time import time
from Bio import SeqIO
from Bio.Alphabet import IUPAC

# Presto and changeo imports
from presto.Defaults import default_out_args
from presto.Annotation import parseAnnotation
from presto.IO import countSeqFile, printLog, printProgress
from changeo.Commandline import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from changeo.IO import countDbFile, extractIMGT, getDbWriter, getRepo
from changeo.Parsers import IgBLASTReader, IMGTReader, gapV, getRegions, inferJunction
from changeo.Receptor import IgRecord, parseAllele, v_allele_regex, d_allele_regex, \
                             j_allele_regex


# Default parameters
default_delimiter = ('\t', ',', '-')


def getInputSeq(seq_file):
    """
    Fetch input sequences for IgBlast or iHMMune-align queries

    Arguments:
    seq_file = a fasta file of sequences input to IgBlast or iHMMune-align

    Returns:
    a dictionary of {ID:Seq}
    """

    seq_dict = SeqIO.index(seq_file, "fasta", IUPAC.ambiguous_dna)

    # Create a seq_dict ID translation using IDs truncate up to space or 50 chars
    seqs = {}
    for seq in seq_dict.values():
        seqs.update({seq.description:str(seq.seq)})

    return seqs

    
# TODO:  should be more readable
def readIHMM(ihmm_output, seq_dict, repo_dict,
             score_fields=False, region_fields=False):
    """
    Reads iHMMuneAlign output

    Arguments:
    ihmm_output = iHMMuneAlign output file
    seq_dict = a dictionary of {ID:Seq} from input fasta file
    repo_dict = dictionary of IMGT gapped germline sequences
    score_fields = if True parse alignment scores
    region_fields = if True add FWR and CDR region fields

    Returns:
    a generator of dictionaries containing alignment data
    """
    # iHMMuneAlign columns
    # Courtesy of Katherine Jackson
    #
    #  1: Identifier - sequence identifer from FASTA input file
    #  2: IGHV - IGHV gene match from the IGHV repertoire, if multiple genes had equally
    #            good alignments both will be listed, if indels were found this will be
    #            listed, in case of multiple IGHV all further data is reported with
    #            respect to the first listed gene
    #  3: IGHD - IGHD gene match, if no IGHD could be found or the IGHD that was found
    #            failed to meet confidence criteria this will be 'NO_DGENE_ALIGNMENT'
    #  4: IGHJ - IGHJ gene match, only a single best matching IGHJ is reported, if indels
    #            are found then 'indel' will be listed
    #  5: V-REGION - portion of input sequence that matches to the germline IGHV, were
    #                nucleotide are missing at start or end the sequence is padded back
    #                to full length with '.' (the exonuclease loss from the end of the
    #                gene will therefore be equal to the number of '.' characters at the
    #                5` end), mismatches between germline and rearranged are in uppercase,
    #                matches are in lowercase
    #  6: N1-REGION - sequence between V- and D-REGIONs
    #  7: D-REGION - portion of input sequence that matches to the germline IGHD
    #                (model doesn't currently permit indels in the IGHD), where IGHD is
    #                reported as 'NO_DGENE_ALIGNMENT' this field contains all nucleotides
    #                between the V- and J-REGIONs
    #  8: N2-REGION - sequence between D- and J-REGIONs
    #  9: J-REGION - portion of the input sequence that matches germline IGHJ, padded
    #                5` and 3` to length of germline match
    # 10: V mutation count - count of mismatches in the V-REGION
    # 11: D mutation count - count of mismatches in the D-REGION
    # 12: J mutation count - count of mismatches in the J-REGION
    # 13: count of ambigious nts - count of 'n' or 'x' nucleotides in the input sequence
    # 14: IGHJ in-frame - 'true' is IGHJ is in-frame and 'false' if IGHJ is out-of-frame,
    #                     WARNING indels and germline IGHV database sequences that are
    #                     not RF1 can cause this to report inaccurately
    # 15: IGHV start offset - offset for start of alignment between input sequence and
    #                         germline IGHV
    #                         NOTE: appears to be base 1 indexing.
    # 16: stop codons - count of stop codons in the sequence, WARNING indels and germline
    #                   IGHV database sequence that are not RF can cause this to be inaccurate
    # 17: IGHD probability - probability that N-nucleotide addition could have created the
    #                        D-REGION sequence
    # 18: HMM path score - path score from HMM
    # 19: reverse complement - 0 for no reverse complement, 1 if alignment was to reverse
    #                          complement NOTE currently this version only functions with
    #                          input in coding orientation
    # 20: mutations in common region - count of mutations in common region, which is a
    #                                  portion of the IGHV that is highly conserved,
    #                                  mutations in this region are used to set various
    #                                  probabilities in the HMM
    # 21: ambigious nts in common region - count of 'n' or 'x' nucleotides in the
    #                                      common region
    # 22: IGHV start offset  - offset for start of alignment between input sequence and
    #                          germline IGHV
    #                          NOTE: appears to be base 0 indexing.
    #                          NOTE: don't know if this differs from 15; it doesn't appear to.
    # 23: IGHV gene length - length of IGHV gene
    # 24: A score - A score probability is calculated from the common region mutations
    #               and is used for HMM calculations relating to expected mutation
    #               probability at different positions in the rearrangement
    record_headers = ['SEQUENCE_ID',
                      'V_CALL',
                      'D_CALL',
                      'J_CALL',
                      'V_SEQ',
                      'NP1_SEQ',
                      'D_SEQ',
                      'NP2_SEQ',
                      'J_SEQ',
                      'V_MUT',
                      'D_MUT',
                      'J_MUT',
                      'NX_COUNT',
                      'J_INFRAME',
                      'V_SEQ_START',
                      'STOP_COUNT',
                      'D_PROB',
                      'HMM_SCORE',
                      'RC',
                      'COMMON_MUT',
                      'COMMON_NX_COUNT',
                      'V_SEQ_START2',
                      'V_SEQ_LENGTH',
                      'A_SCORE']

    record_handle = open(ihmm_output, 'rU')
    csv.register_dialect('semicol', delimiter=';', quotechar='"')
    records = csv.DictReader(record_handle, fieldnames=record_headers, dialect='semicol')

    # loop through records
    for row in records:
        # Skip empty lines
        if not row:
            continue

        # Process row if not empty
        db = {'SEQUENCE_ID': row['SEQUENCE_ID'],
              'SEQUENCE_INPUT': seq_dict[row['SEQUENCE_ID']]}

        if not row['V_CALL'] or row['V_CALL'].startswith('NA - ') or \
                row['V_CALL'].startswith('State path'):
            db['FUNCTIONAL'] = None
            db['V_CALL'] = None
            yield IgRecord(db)
            continue

        # Check stop codons
        db['FUNCTIONAL'] = 'F' if int(row['STOP_COUNT']) > 0 else 'T'

        # Check whether J is in-frame
        if row['J_INFRAME'] != 'true':
            db['IN_FRAME'] = 'F'
            db['FUNCTIONAL'] = 'F'

        # Check for indels
        if re.search('\[indels\]', row['V_CALL'] + row['D_CALL'] + row['J_CALL']):
            db['INDELS'] = 'T'
        else:
            db['INDELS'] = 'F'

        # Find V-REGION
        v_call = parseAllele(row['V_CALL'], v_allele_regex, action='list')
        vkey = (v_call[0],)
        if vkey in repo_dict:
            db['V_CALL'] = ','.join(v_call) if v_call is not None else 'None'
            sample_vseq = (row['V_SEQ']).strip('.')
            db['V_SEQ_START'] = int(row['V_SEQ_START'])
            db['V_SEQ_LENGTH'] = len(sample_vseq)
        else:
            yield IgRecord(db)
            continue

        # Find D-REGION
        if not row['D_CALL'] or row['D_CALL'] == 'NO_DGENE_ALIGNMENT':
            dgene, sample_dseq = 'NA', ''
        else:
            d_call = parseAllele(row['D_CALL'], d_allele_regex, action='list')
            dkey = (d_call[0],)
            if dkey in repo_dict:
                db['D_CALL'] = ','.join(d_call) if d_call is not None else 'None'
                sample_dseq = (row['D_SEQ']).strip('.')
                db['D_GERM_START'] = len(row['D_SEQ']) - len((row['D_SEQ']).lstrip('.'))
                db['D_GERM_LENGTH'] = len(sample_dseq)
            else:
                yield IgRecord(db)
                continue

        # Find J-REGION
        if not row['J_CALL'] or row['J_CALL'] == 'NO_JGENE_ALIGNMENT':
            db['FUNCTIONAL'] = 'F'
            db['J_CALL'] = ''
            yield IgRecord(db)
            continue
        else:
            j_call = parseAllele(row['J_CALL'], j_allele_regex, action='list')
            jkey = (j_call[0],)
            if jkey in repo_dict:
                db['J_CALL'] = ','.join(j_call) if j_call is not None else 'None'
                sample_jseq = (row['J_SEQ']).strip('.')
                db['J_GERM_START'] = len(row['J_SEQ']) - len((row['J_SEQ']).lstrip('.'))
                db['J_GERM_LENGTH'] = len(sample_jseq)
            else:
                yield IgRecord(db)
                continue

        # Assemble sample sequence
        db['SEQUENCE_VDJ'] = ''.join([sample_vseq,
                                      row['NP1_SEQ'],
                                      sample_dseq,
                                      row['NP2_SEQ'],
                                      sample_jseq])

        # Germline positions
        db['V_GERM_START_VDJ'] = 1
        db['V_GERM_LENGTH_VDJ'] = len(sample_vseq)

        # N/P lengths
        db['NP1_LENGTH'] = len(row['NP1_SEQ'])
        db['NP2_LENGTH'] = len(row['NP2_SEQ'])

        # D positions
        db['D_SEQ_START'] = sum(int(i) for i in [db['V_SEQ_START'],
                                                 db['V_SEQ_LENGTH'],
                                                 db['NP1_LENGTH']] if i)
        db['D_SEQ_LENGTH'] = len(sample_dseq)

        # J positions
        db['J_SEQ_START'] = sum(int(i) for i in [db['V_SEQ_START'],
                                                 db['V_SEQ_LENGTH'],
                                                 db['NP1_LENGTH'],
                                                 db['D_SEQ_LENGTH'],
                                                 db['NP2_LENGTH']] if i)
        db['J_SEQ_LENGTH'] = len(sample_jseq)

        # Add IMGT gapped sequence
        if v_call:
            db.update(gapV(db, repo_dict))

        # Extract junction regions
        if v_call and j_call:
            db.update(inferJunction(db, repo_dict))

         # Overall alignment score
        if score_fields:
            try: db['HMM_SCORE'] = float(row['HMM_SCORE'])
            except (TypeError, ValueError): db['HMM_SCORE'] = ''

        # FWR and CDR regions
        if region_fields:
            db.update(getRegions(db))

        yield IgRecord(db)


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


def writeDb(db, file_prefix, total_count, id_dict={}, no_parse=True, partial=False,
            score_fields=False, region_fields=False, junction_fields=False,
            out_args=default_out_args):
    """
    Writes tab-delimited database file in output directory
    
    Arguments:
    db = a iterator of IgRecord objects containing alignment data
    file_prefix = directory and prefix for CLIP tab-delim file
    total_count = number of records (for progress bar)
    id_dict = a dictionary of {IMGT ID: full seq description}
    no_parse = if ID is to be parsed for pRESTO output with default delimiters
    partial = If True put incomplete alignments in the pass file
    score_fields = if True add alignment score fields to output file
    region_fields = if True add FWR and CDR region fields to output file
    junction_fields = if True add D FRAME junction field to output file
    out_args = common output argument dictionary from parseCommonArgs

    Returns:
    None
    """
    # Function to check for valid records strictly
    def _pass_strict(rec):
        valid = [rec.v_call != 'None',
                 rec.j_call != 'None',
                 rec.functional is not None,
                 rec.seq_vdj,
                 rec.junction]
        return all(valid)

    # Function to check for valid records loosely
    def _pass_gentle(rec):
        valid = [rec.v_call != 'None',
                 rec.d_call != 'None',
                 rec.j_call != 'None']
        return any(valid)

    # Set pass criteria
    _pass = _pass_gentle if partial else _pass_strict

    # Define output file names
    pass_file = '%s_db-pass.tab' % file_prefix
    fail_file = '%s_db-fail.tab' % file_prefix

    # Define core fields
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
                      'V_GERM_START_VDJ',
                      'V_GERM_LENGTH_VDJ',
                      'V_GERM_START_IMGT',
                      'V_GERM_LENGTH_IMGT',
                      'NP1_LENGTH',
                      'D_SEQ_START',
                      'D_SEQ_LENGTH',
                      'D_GERM_START',
                      'D_GERM_LENGTH',
                      'NP2_LENGTH',
                      'J_SEQ_START',
                      'J_SEQ_LENGTH',
                      'J_GERM_START',
                      'J_GERM_LENGTH',
                      'JUNCTION_LENGTH',
                      'JUNCTION']

    # Define optional scoring fields
    if score_fields:
        ordered_fields.extend(['V_SCORE',
                               'V_IDENTITY',
                               'V_EVALUE',
                               'V_BTOP',
                               'J_SCORE',
                               'J_IDENTITY',
                               'J_EVALUE',
                               'J_BTOP',
                               'HMM_SCORE'])

    # Define optional region fields
    if region_fields:
        ordered_fields.extend(['FWR1_IMGT', 'FWR2_IMGT', 'FWR3_IMGT', 'FWR4_IMGT',
                               'CDR1_IMGT', 'CDR2_IMGT', 'CDR3_IMGT'])

    # Define optional junction fields
    if junction_fields:
        ordered_fields.extend(['N1_LENGTH', 'N2_LENGTH', 
                               'P3V_LENGTH', 'P5D_LENGTH', 'P3D_LENGTH', 'P5J_LENGTH',
                               'D_FRAME'])

    # Not currently implemented
    # if ihmm_germ:
    #     ordered_fields.extend(['GERMLINE_IHMM', 'GERMLINE_IHMM_D_MASK'])

    # Initiate handles, writers and counters
    pass_handle = None
    fail_handle = None
    pass_writer = None
    fail_writer = None
    start_time = time()
    rec_count = pass_count = fail_count = 0

    # Validate and write output
    for record in db:
        printProgress(rec_count, total_count, 0.05, start_time)
        rec_count += 1

        # Replace sequence description with full string, if required
        if id_dict and record.id in id_dict:
            record.id = id_dict[record.id]

        # Parse sequence description into new columns
        if not no_parse:
            record.annotations = parseAnnotation(record.id, delimiter=out_args['delimiter'])
            record.id = record.annotations['ID']
            del record.annotations['ID']

            # TODO:  This is not the best approach. should pass in output fields.
            # If first record, use parsed description to define extra columns
            if pass_writer is None and fail_writer is None:
                ordered_fields.extend(list(record.annotations.keys()))

        # Count pass or fail and write to appropriate file
        if _pass(record):
            # Open pass file
            if pass_writer is None:
                pass_handle = open(pass_file, 'wt')
                pass_writer = getDbWriter(pass_handle, add_fields=ordered_fields)

            # Write row to pass file
            pass_count += 1
            pass_writer.writerow(record.toDict())
        else:
            # Open failed file
            if out_args['failed'] and fail_writer is None:
                fail_handle = open(fail_file, 'wt')
                fail_writer = getDbWriter(fail_handle, add_fields=ordered_fields)

            # Write row to fail file if specified
            fail_count += 1
            if fail_writer is not None:
                fail_writer.writerow(record.toDict())

    # Print log
    printProgress(rec_count, total_count, 0.05, start_time)

    log = OrderedDict()
    log['OUTPUT'] = pass_file
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    log['END'] = 'MakeDb'
    printLog(log)
    
    if pass_handle is not None: pass_handle.close()
    if fail_handle is not None: fail_handle.close()


# TODO:  may be able to merge with other mains
def parseIgBLAST(igblast_output, seq_file, repo, no_parse=True, partial=False,
                 score_fields=False, region_fields=False, out_args=default_out_args):
    """
    Main for IgBLAST aligned sample sequences

    Arguments:
      igblast_output : IgBLAST output file to process
      seq_file : fasta file input to IgBlast (from which to get sequence)
      repo : folder with germline repertoire files
      no_parse : if ID is to be parsed for pRESTO output with default delimiters
      partial : If True put incomplete alignments in the pass file
      score_fields : if True add alignment score fields to output file
      region_fields : if True add FWR and CDR region fields to output file
      out_args : common output argument dictionary from parseCommonArgs

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
    log['PARTIAL'] = partial
    log['SCORE_FIELDS'] = score_fields
    log['REGION_FIELDS'] = region_fields
    printLog(log)

    # Get input sequence dictionary
    seq_dict = getInputSeq(seq_file)

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

    # Create germline repo dictionary
    repo_dict = getRepo(repo)

    # Parse and write output
    with open(igblast_output, 'r') as f:
        parse_iter = IgBLASTReader(f, seq_dict, repo_dict, score_fields=score_fields,
                                   region_fields=region_fields)
        writeDb(parse_iter, file_prefix, total_count, no_parse=no_parse,
                score_fields=score_fields, region_fields=region_fields, out_args=out_args)

    return None


# TODO:  may be able to merge with other mains
def parseIMGT(imgt_output, seq_file=None, no_parse=True, partial=False,
              score_fields=False, region_fields=False, junction_fields=False,
              out_args=default_out_args):
    """
    Main for IMGT aligned sample sequences

    Arguments:
    imgt_output = zipped file or unzipped folder output by IMGT
    seq_file = FASTA file input to IMGT (from which to get seqID)
    no_parse = if ID is to be parsed for pRESTO output with default delimiters
    partial = If True put incomplete alignments in the pass file
    score_fields = if True add alignment score fields to output file
    region_fields = if True add FWR and CDR region fields to output file
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
    log['PARTIAL'] = partial
    log['SCORE_FIELDS'] = score_fields
    log['REGION_FIELDS'] = region_fields
    log['JUNCTION_FIELDS'] = junction_fields
    printLog(log)
        
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

    # Extract IMGT files
    temp_dir, imgt_files = extractIMGT(imgt_output)
    total_count = countDbFile(imgt_files['summary'])

    # Get (parsed) IDs from fasta file submitted to IMGT
    id_dict = getIDforIMGT(seq_file) if seq_file else {}
    
    # Parse IMGT output and write db
    with open(imgt_files['summary'], 'r') as summary_handle, \
            open(imgt_files['gapped'], 'r') as gapped_handle, \
            open(imgt_files['ntseq'], 'r') as ntseq_handle, \
            open(imgt_files['junction'], 'r') as junction_handle:
        parse_iter = IMGTReader(summary_handle, gapped_handle, ntseq_handle, junction_handle,
                                score_fields=score_fields, region_fields=region_fields,
                                junction_fields=junction_fields)
        writeDb(parse_iter, file_prefix, total_count, id_dict=id_dict, no_parse=no_parse,
                score_fields=score_fields, region_fields=region_fields, junction_fields=junction_fields,
                out_args=out_args)

    # Cleanup temp directory
    temp_dir.cleanup()

    return None


# TODO:  may be able to merge with other mains
def parseIHMM(ihmm_output, seq_file, repo, no_parse=True, partial=False,
              score_fields=False, region_fields=False, out_args=default_out_args):
    """
    Main for iHMMuneAlign aligned sample sequences

    Arguments:
    ihmm_output = iHMMuneAlign output file to process
    seq_file = fasta file input to iHMMuneAlign (from which to get sequence)
    repo = folder with germline repertoire files
    no_parse = if ID is to be parsed for pRESTO output with default delimiters
    partial = If True put incomplete alignments in the pass file
    score_fields = if True parse alignment scores
    region_fields = if True add FWR and CDR region fields
    out_args = common output argument dictionary from parseCommonArgs

    Returns:
    None
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'MakeDB'
    log['ALIGNER'] = 'iHMMuneAlign'
    log['ALIGN_RESULTS'] = os.path.basename(ihmm_output)
    log['SEQ_FILE'] = os.path.basename(seq_file)
    log['NO_PARSE'] = no_parse
    log['PARTIAL'] = partial
    log['SCORE_FIELDS'] = score_fields
    log['REGION_FIELDS'] = region_fields
    printLog(log)

    # Get input sequence dictionary
    seq_dict = getInputSeq(seq_file)

    # Formalize out_dir and file-prefix
    if not out_args['out_dir']:
        out_dir = os.path.split(ihmm_output)[0]
    else:
        out_dir = os.path.abspath(out_args['out_dir'])
        if not os.path.exists(out_dir):  os.mkdir(out_dir)
    if out_args['out_name']:
        file_prefix = out_args['out_name']
    else:
        file_prefix = os.path.basename(os.path.splitext(ihmm_output)[0])
    file_prefix = os.path.join(out_dir, file_prefix)

    total_count = countSeqFile(seq_file)

    # Create
    repo_dict = getRepo(repo)
    ihmm_dict = readIHMM(ihmm_output, seq_dict, repo_dict,
                         score_fields=score_fields, region_fields=region_fields)
    writeDb(ihmm_dict, file_prefix, total_count,
            score_fields=score_fields, region_fields=region_fields,
            no_parse=no_parse, out_args=out_args)

    return None


def getArgParser():
    """
    Defines the ArgumentParser

    Arguments: 
    None
                      
    Returns: 
    an ArgumentParser object
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
                  J_EVALUE, J_SCORE, J_IDENTITY, J_BTOP

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
    parser_igblast.add_argument('-i', nargs='+', action='store', dest='aligner_files',
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
                                help='''Specify if input IDs should not be parsed to add
                                    new columns to database.''')
    parser_igblast.add_argument('--partial', action='store_true', dest='partial',
                                help='''If specified, include incomplete V(D)J alignments in
                                     the pass file instead of the fail file.''')
    parser_igblast.add_argument('--scores', action='store_true', dest='score_fields',
                                help='''Specify if alignment score metrics should be
                                     included in the output. Adds the V_SCORE, V_IDENTITY,
                                     V_EVALUE, V_BTOP, J_SCORE, J_IDENTITY,
                                     J_BTOP, and J_EVALUE columns.''')
    parser_igblast.add_argument('--regions', action='store_true', dest='region_fields',
                                help='''Specify if IMGT framework and CDR regions should be
                                     included in the output. Adds the FWR1_IMGT, FWR2_IMGT,
                                     FWR3_IMGT, FWR4_IMGT, CDR1_IMGT, CDR2_IMGT, and
                                     CDR3_IMGT columns.''')
    parser_igblast.set_defaults(func=parseIgBLAST)

    # IMGT aligner
    parser_imgt = subparsers.add_parser('imgt', parents=[parser_parent],
                                        formatter_class=CommonHelpFormatter,
                                        help='''Process IMGT/HighV-Quest output
                                             (does not work with V-QUEST).''',
                                        description='''Process IMGT/HighV-Quest output
                                             (does not work with V-QUEST).''')
    parser_imgt.add_argument('-i', nargs='+', action='store', dest='aligner_files',
                             help='''Either zipped IMGT output files (.zip or .txz) or a
                                  folder containing unzipped IMGT output files (which must
                                  include 1_Summary, 2_IMGT-gapped, 3_Nt-sequences,
                                  and 6_Junction).''')
    parser_imgt.add_argument('-s', nargs='*', action='store', dest='seq_files',
                             required=False,
                             help='''List of input FASTA files (with .fasta, .fna or .fa
                                  extension) containing sequences.''')
    parser_imgt.add_argument('--noparse', action='store_true', dest='no_parse', 
                             help='''Specify if input IDs should not be parsed to add new
                                  columns to database.''')
    parser_imgt.add_argument('--partial', action='store_true', dest='partial',
                             help='''If specified, include incomplete V(D)J alignments in
                                  the pass file instead of the fail file.''')
    parser_imgt.add_argument('--scores', action='store_true', dest='score_fields',
                             help='''Specify if alignment score metrics should be
                                  included in the output. Adds the V_SCORE, V_IDENTITY,
                                  J_SCORE and J_IDENTITY. Note, this will also add
                                  the columns V_EVALUE, V_BTOP, J_EVALUE and J_BTOP,
                                  but they will be empty for IMGT output.''')
    parser_imgt.add_argument('--regions', action='store_true', dest='region_fields',
                             help='''Specify if IMGT framework and CDR regions should be
                                  included in the output. Adds the FWR1_IMGT, FWR2_IMGT,
                                  FWR3_IMGT, FWR4_IMGT, CDR1_IMGT, CDR2_IMGT, and
                                  CDR3_IMGT columns.''')
    parser_imgt.add_argument('--junction', action='store_true', dest='junction_fields',
                             help='''Specify if junction fields should be
                                  included in the output. Adds the columns 
                                  N1_LENGTH, N2_LENGTH, P3V_LENGTH, P5D_LENGTH, P3D_LENGTH,
                                  P5J_LENGTH, D_FRAME.''')
    parser_imgt.set_defaults(func=parseIMGT)

    # iHMMuneAlign Aligner
    parser_ihmm = subparsers.add_parser('ihmm', parents=[parser_parent],
                                        formatter_class=CommonHelpFormatter,
                                        help='Process iHMMuneAlign output.',
                                        description='Process iHMMuneAlign output.')
    parser_ihmm.add_argument('-i', nargs='+', action='store', dest='aligner_files',
                             required=True,
                             help='''iHMMuneAlign output file.''')
    parser_ihmm.add_argument('-r', nargs='+', action='store', dest='repo', required=True,
                             help='''List of folders and/or FASTA files containing
                                  IMGT-gapped germline sequences corresponding to the
                                  set of germlines used in the IgBLAST alignment.''')
    parser_ihmm.add_argument('-s', action='store', nargs='+', dest='seq_files',
                             required=True,
                             help='''List of input FASTA files (with .fasta, .fna or .fa
                                  extension) containing sequences.''')
    parser_ihmm.add_argument('--noparse', action='store_true', dest='no_parse',
                             help='''Specify if input IDs should not be parsed to add
                                  new columns to database.''')
    parser_ihmm.add_argument('--partial', action='store_true', dest='partial',
                             help='''If specified, include incomplete V(D)J alignments in
                                  the pass file instead of the fail file.''')
    parser_ihmm.add_argument('--scores', action='store_true', dest='score_fields',
                             help='''Specify if alignment score metrics should be
                                  included in the output. Adds the path score of the
                                  iHMMuneAlign hidden Markov model to HMM_SCORE.''')
    parser_ihmm.add_argument('--regions', action='store_true', dest='region_fields',
                             help='''Specify if IMGT framework and CDR regions should be
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
    elif args.command == 'ihmm':
        for i in range(len(args.__dict__['aligner_files'])):
            args_dict['ihmm_output'] =  args.__dict__['aligner_files'][i]
            args_dict['seq_file'] = args.__dict__['seq_files'][i]
            args.func(**args_dict)
