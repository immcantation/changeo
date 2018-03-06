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
from itertools import groupby
from textwrap import dedent
from time import time

# Presto and change imports
from presto.Defaults import default_out_args
from presto.IO import getOutputHandle, printLog, printMessage, printProgress
from changeo.Defaults import default_v_field, default_d_field, default_j_field, default_clone_field, \
                             default_format
from changeo.Commandline import CommonHelpFormatter, checkArgs, getCommonArgParser, parseCommonArgs
from changeo.IO import countDbFile, getDbFields, readGermlines
from changeo.Parsers import AIRRReader, AIRRWriter, AIRRSchema, ChangeoReader, ChangeoWriter, ChangeoSchema

# Defaults
default_seq_field = 'SEQUENCE_IMGT'
default_germ_types = 'dmask'

# TODO: this is not generalized for non-IMGT gapped sequences!
def getVGermline(receptor, references, v_field=default_v_field):
    """
    Extract V allele and germline sequence

    Arguments:
      receptor : Receptor object
      references : dictionary of germline sequences
      v_field : field containing the V allele assignment

    Returns:
      tuple : V allele name, V segment germline sequence
    """
    # Extract V allele call
    vgene = receptor.getVAllele(action='first', field=v_field)

    # Build V segment germline sequence
    if vgene is None:
        try:  vlen = int(receptor.v_germ_length_imgt)
        except (TypeError, ValueError):  vlen = 0
        germ_vseq = 'N' * vlen
    elif vgene in references:
        # Define V germline positions
        try:  vstart = int(receptor.v_germ_start_imgt) - 1
        except (TypeError, ValueError):  vstart = 0
        try:  vlen = int(receptor.v_germ_length_imgt)
        except (TypeError, ValueError):  vlen = 0
        # Define V germline sequence
        vseq = references[vgene]
        vpad = vlen - len(vseq[vstart:])
        if vpad < 0: vpad = 0
        germ_vseq = vseq[vstart:(vstart + vlen)] + ('N' * vpad)
    else:
        germ_vseq = None

    return vgene, germ_vseq


def getDGermline(receptor, references, d_field=default_d_field):
    """
    Extract D allele and germline sequence

    Arguments:
      receptor : Receptor object
      references : dictionary of germline sequences
      d_field : field containing the D allele assignment

    Returns:
      tuple : D allele name, D segment germline sequence
    """
    # Extract D allele call
    dgene = receptor.getDAllele(action='first', field=d_field)

    # Build D segment germline sequence
    if dgene is None:
        germ_dseq = ''
    elif dgene in references:
        # Define D germline positions
        try:  dstart = int(receptor.d_germ_start) - 1
        except (TypeError, ValueError):  dstart = 0
        try:  dlen = int(receptor.d_germ_length)
        except (TypeError, ValueError):  dlen = 0
        # Define D germline sequence
        dseq = references[dgene]
        germ_dseq = dseq[dstart:(dstart + dlen)]

    return dgene, germ_dseq


def getJGermline(receptor, references, j_field=default_j_field):
    """
    Extract J allele and germline sequence

    Arguments:
      receptor : Receptor object
      references : dictionary of germline sequences
      j_field : field containing the J allele assignment

    Returns:
      tuple : J allele name, J segment germline sequence
    """
    # Extract J allele call
    jgene = receptor.getJAllele(action='first', field=j_field)

    # Build J segment germline sequence
    if jgene is None:
        try:  jlen = int(receptor.j_germ_length)
        except (TypeError, ValueError):  jlen = 0
        germ_jseq = 'N' * jlen
    elif jgene in references:
        jseq = references[jgene]
        # Define J germline positions
        try:  jstart = int(receptor.j_germ_start) - 1
        except (TypeError, ValueError):  jstart = 0
        try:  jlen = int(receptor.j_germ_length)
        except (TypeError, ValueError):  jlen = 0
        # Define J germline sequence
        jpad = jlen - len(jseq[jstart:])
        if jpad < 0: jpad = 0
        germ_jseq = jseq[jstart:(jstart + jlen)] + ('N' * jpad)

    return jgene, germ_jseq


def stitchVDJ(receptor, v_seq, d_seq, j_seq):
    """
    Assemble full length germline sequence

    Arguments:
      receptor : Receptor object
      v_seq : V segment sequence as a string
      d_seq : D segment sequence as a string
      j_seq : J segment sequence as a string

    Returns:
      str : full germline sequence
    """
    # Assemble pieces starting with V segment
    sequence = v_seq

    # Add Ns for first N/P region
    try:  np1_len = int(receptor.np1_length)
    except (TypeError, ValueError):  np1_len = 0
    sequence += 'N' * np1_len

    # Add D segment
    sequence += d_seq

    # Add Ns for second N/P region
    try:  np2_len = int(receptor.np2_length)
    except (TypeError, ValueError):  np2_len = 0
    sequence += 'N' * np2_len

    # Add J segment
    sequence += j_seq

    return sequence


def stitchRegions(receptor, v_seq, d_seq, j_seq):
    """
    Assemble full length region encoding

    Arguments:
      receptor : Receptor object
      v_seq : V segment germline sequence as a string
      d_seq : D segment germline sequence as a string
      j_seq : J segment germline sequence as a string

    Returns:
      str : string defining germline regions
    """
    # Set mode for region definitions
    full_junction = True if getattr(receptor, 'n1_length', None) is not None else False

    # Assemble pieces starting with V segment
    regions = 'V' * len(v_seq)

    # NP nucleotide additions after V
    if not full_junction:
        # PNP nucleotide additions after V
        try:  np1_len = int(receptor.np1_length)
        except (TypeError, ValueError):  np1_len = 0
        regions += 'N' * np1_len
    else:
        # P nucleotide additions before N1
        try:  p3v_len = int(receptor.p3v_length)
        except (TypeError, ValueError):  p3v_len = 0
        # N1 nucleotide additions
        try:  n1_len = int(receptor.n1_length)
        except (TypeError, ValueError):  n1_len = 0
        # P nucleotide additions before D
        try:  p5d_len = int(receptor.p5d_length)
        except (TypeError, ValueError):  p5d_len = 0

        # Update regions
        regions += 'P' * p3v_len
        regions += 'N' * n1_len
        regions += 'P' * p5d_len

    # Add D segment
    regions += 'D' * len(d_seq)

    # NP nucleotide additions before J
    if not full_junction:
        # NP nucleotide additions
        try:  np2_len = int(receptor.np2_length)
        except (TypeError, ValueError):  np2_len = 0
        regions += 'N' * np2_len
    else:
        # P nucleotide additions after D
        try: p3d_len = int(receptor.p3d_length)
        except (TypeError, ValueError): p3d_len = 0
        # N2 nucleotide additions
        try:  n2_len = int(receptor.n2_length)
        except (TypeError, ValueError): n2_len = 0
        # P nucleotide additions before J
        try:  p5j_len = int(receptor.p5j_length)
        except (TypeError, ValueError):  p5j_len = 0

        # Update regions
        regions += 'P' * p3d_len
        regions += 'N' * n2_len
        regions += 'P' * p5j_len

    # Add J segment
    regions += 'J' * len(j_seq)

    return regions


def buildGermline(receptor, references, seq_field=default_seq_field, v_field=default_v_field,
                  d_field=default_d_field, j_field=default_j_field):
    """
    Join gapped germline sequences aligned with sample sequences

    Arguments:
      receptor : Receptor object
      references : dictionary of IMGT gapped germline sequences
      seq_field : field in which to look for sequence
      v_field : field in which to look for V call
      d_field : field in which to look for V call
      j_field : field in which to look for V call

    Returns:
      tuple : log dictionary, dictionary of {germline_type: germline_sequence}, dictionary of {segment: gene call}
    """
    # Return objects
    log = OrderedDict()
    germlines = {'full': '', 'dmask': '', 'vonly': '', 'regions': ''}

    # Build V segment germline sequence
    vgene, germ_vseq = getVGermline(receptor, references, v_field=v_field)
    log['V_CALL'] = vgene
    if germ_vseq is None:
        log['ERROR'] = 'Allele %s in not in the provided germline database.' % vgene
        return log, None, None

    # Build D segment germline sequence
    dgene, germ_dseq = getDGermline(receptor, references, d_field=d_field)
    log['D_CALL'] = dgene
    if germ_dseq is None:
        log['ERROR'] = 'Allele %s in not in the provided germline database.' % vgene
        return log, None, None

    # Build J segment germline sequence
    jgene, germ_jseq = getJGermline(receptor, references, j_field=j_field)
    log['J_CALL'] = jgene
    if germ_jseq is None:
        log['ERROR'] = 'Allele %s in not in the provided germline database.' % vgene
        return log, None, None

    # Stitch complete germlines
    germ_seq = stitchVDJ(receptor, germ_vseq, germ_dseq, germ_jseq)
    regions = stitchRegions(receptor, germ_vseq, germ_dseq, germ_jseq)

    # Update log
    log['SEQUENCE'] = receptor.getField(seq_field)
    log['GERMLINE'] = germ_seq
    log['REGIONS'] = regions

    # Check that input and germline sequence match
    if len(receptor.getField(seq_field)) == 0:
        log['ERROR'] = 'Sequence is missing from the %s field' % seq_field
        return log, None, None

    if len(germ_seq) != len(receptor.getField(seq_field)):
        log['ERROR'] = 'Germline sequence is %d nucleotides longer than input sequence' % \
                              (len(germlines['full']) - len(receptor.getField(seq_field)))
        return log, None, None

    # Define return germlines object
    germ_dmask = germ_seq[:len(germ_vseq)] + \
                 'N' * (len(germ_seq) - len(germ_vseq) - len(germ_jseq)) + \
                 germ_seq[-len(germ_jseq):]
    germlines = {'full': germ_seq, 'dmask': germ_dmask, 'vonly': germ_vseq, 'regions': regions}
    for k, v in germlines.items():  germlines[k] = v.upper()

    # Define return genes object
    genes = {'v': log['V_CALL'],
             'd': log['D_CALL'],
             'j': log['J_CALL']}

    return log, germlines, genes


# TODO: Should do 'first' method for ambiguous V/J groups. And explicit allele extraction.
def buildClonalGermline(receptors, references, seq_field=default_seq_field,
                        v_field=default_v_field, d_field=default_d_field, j_field=default_j_field):
    """
    Determine consensus clone sequence and create germline for clone

    Arguments:
      receptors : list of Receptor objects
      references : dictionary of IMGT gapped germline sequences
      seq_field : field in which to look for sequence
      v_field : field in which to look for V call
      d_field : field in which to look for D call
      j_field : field in which to look for J call

    Returns:
      tuple : log dictionary, dictionary of {germline_type: germline_sequence},
              dictionary of consensus {segment: gene call}
    """
    # Log
    log = OrderedDict()

    # Create dictionaries to count observed V/J calls
    v_dict = OrderedDict()
    j_dict = OrderedDict()

    # Find longest sequence in clone
    max_length = 0
    for rec in receptors:
        v = rec.getVAllele(action='first', field=v_field)
        v_dict[v] = v_dict.get(v, 0) + 1
        j = rec.getJAllele(action='first', field=j_field)
        j_dict[j] = j_dict.get(j, 0) + 1
        seq_len = len(rec.getField(seq_field))
        if seq_len > max_length:
            max_length = seq_len

    # Consensus V and J having most observations
    v_cons = [k for k in list(v_dict.keys()) if v_dict[k] == max(v_dict.values())]
    j_cons = [k for k in list(j_dict.keys()) if j_dict[k] == max(j_dict.values())]

    # Consensus sequence(s) with consensus V/J calls and longest sequence
    cons = [x for x in receptors if x.getVAllele(action='first', field=v_field) in v_cons and \
                                    x.getJAllele(action='first', field=j_field) in j_cons and \
                                    len(x.getField(seq_field)) == max_length]
    # Consensus sequence(s) with consensus V/J calls but not the longest sequence
    if not cons:
        cons = [x for x in receptors if x.getVAllele(action='first', field=v_field) in v_cons and \
                                        x.getJAllele(action='first', field=j_field) in j_cons]

    # Return without germline if no sequence has both consensus V and J call
    if not cons:
        log['V_CALL'] = ','.join(v_cons)
        log['J_CALL'] = ','.join(j_cons)
        log['ERROR'] = 'No sequence found with both consensus V and J calls.'
        return log, None, None

    # Select consensus Receptor, resolving ties by alphabetical ordering of sequence id.
    cons = sorted(cons, key=lambda x: x.sequence_id)[0]

    # Pad end of consensus sequence with gaps to make it the max length
    gap_length = max_length - len(cons.getField(seq_field))
    if gap_length > 0:
        cons.j_germ_length = int(cons.j_germ_length or 0) + gap_length
        cons.setField(seq_field, cons.getSeq(seq_field) + ('N' * gap_length))

    # Update lengths padded to longest sequence in clone
    for rec in receptors:
        x = max_length - len(rec.getField(seq_field))
        rec.j_germ_length = int(rec.j_germ_length or 0) + x
        rec.setField(seq_field, rec.getSeq(seq_field) + ('N' * x))

    # Stitch consensus germline
    cons_log, germlines, genes = buildGermline(cons, references, seq_field=seq_field, v_field=v_field,
                                               d_field=d_field, j_field=j_field)

    # Update log
    log['CONSENSUS'] = cons.sequence_id
    log.update(cons_log)

    # Return log
    return log, germlines, genes


def createGermlines(db_file, repo, seq_field=default_seq_field, v_field=default_v_field,
                    d_field=default_d_field, j_field=default_j_field,
                    cloned=False, clone_field=default_clone_field, germ_types=default_germ_types,
                    format=default_format, out_args=default_out_args):
    """
    Write germline sequences to tab-delimited database file

    Arguments:
      db_file : input tab-delimited database file
      repo : folder with germline repertoire files
      seq_field : field in which to look for sequence
      v_field : field in which to look for V call
      d_field : field in which to look for D call
      j_field : field in which to look for J call
      cloned : if True build germlines by clone, otherwise build individual germlines
      clone_field : field containing clone identifiers; ignored if cloned=False.
      germ_types : list of germline sequence types to be output from the set of 'full', 'dmask', 'vonly', 'regions'
      format : input and output format
      out_args : arguments for output preferences

    Returns:
      None
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'CreateGermlines'
    log['DB_FILE'] = os.path.basename(db_file)
    log['GERM_TYPES'] = ','.join(germ_types)
    log['SEQ_FIELD'] = seq_field
    log['V_FIELD'] = v_field
    log['D_FIELD'] = d_field
    log['J_FIELD'] = j_field
    log['CLONED'] = cloned
    if cloned:  log['CLONE_FIELD'] = clone_field
    printLog(log)

    # Format options
    if format == 'changeo':
        reader = ChangeoReader
        writer = ChangeoWriter
        schema = ChangeoSchema
        germline_fields = OrderedDict()
        seq_type = seq_field.split('_')[-1]
        if 'full' in germ_types:  germline_fields['full'] = 'GERMLINE_' + seq_type
        if 'dmask' in germ_types:  germline_fields['dmask'] = 'GERMLINE_' + seq_type + '_D_MASK'
        if 'vonly' in germ_types:  germline_fields['vonly'] = 'GERMLINE_' + seq_type + '_V_REGION'
        if 'regions' in germ_types:  germline_fields['regions'] = 'GERMLINE_REGIONS'
        if cloned:
            germline_fields['v'] = 'GERMLINE_V_CALL'
            germline_fields['d'] = 'GERMLINE_D_CALL'
            germline_fields['j'] = 'GERMLINE_J_CALL'
        out_fields = getDbFields(db_file, add=list(germline_fields.values()), reader=reader)
    elif format == 'airr':
        reader = AIRRReader
        writer = AIRRWriter
        schema = AIRRSchema
        germline_fields = OrderedDict()
        # TODO: this won't work for AIRR necessarily
        seq_type = seq_field.split('_')[-1]
        if 'full' in germ_types:  germline_fields['full'] = 'germline_' + seq_type
        if 'dmask' in germ_types:  germline_fields['dmask'] = 'germline_' + seq_type + '_d_mask'
        if 'vonly' in germ_types:  germline_fields['vonly'] = 'germline_' + seq_type + '_v_region'
        if 'regions' in germ_types:  germline_fields['regions'] = 'germline_regions'
        if cloned:
            germline_fields['v'] = 'germline_v_call'
            germline_fields['d'] = 'germline_d_call'
            germline_fields['j'] = 'germline_j_call'
        out_fields = getDbFields(db_file, add=list(germline_fields.values()), reader=reader)
    else:
        sys.exit('Error:  Invalid format %s' % format)

    # Get repertoire and open Db reader
    references = readGermlines(repo)
    db_handle = open(db_file, 'rt')
    reader = reader(db_handle)

    # Check for existence of fields
    for f in [v_field, d_field, j_field, seq_field]:
        if f not in reader.fields:
            sys.exit('Error: %s field does not exist in input database file.' % f)
    #  Translate to Receptor attribute names
    v_field = schema.asReceptor(v_field)
    d_field = schema.asReceptor(d_field)
    j_field = schema.asReceptor(j_field)
    seq_field = schema.asReceptor(seq_field)
    clone_field = schema.asReceptor(clone_field)

    # Define log handle
    if out_args['log_file'] is None:  log_handle = None
    else:  log_handle = open(out_args['log_file'], 'w')

    # Create output file handle and Db writer
    pass_handle = getOutputHandle(db_file,
                                  out_label='germ-pass',
                                  out_dir=out_args['out_dir'],
                                  out_name=out_args['out_name'],
                                  out_type='tsv')
    pass_writer = writer(pass_handle, fields=out_fields)

    if out_args['failed']:
        fail_handle = getOutputHandle(db_file,
                                      out_label='germ-fail',
                                      out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'],
                                      out_type='tsv')
        fail_writer = writer(fail_handle, fields=out_fields)
    else:
        fail_handle = None
        fail_writer = None

    # Define Receptor iterator
    if cloned:
        start_time = time()
        printMessage('Sorting by clone', start_time=start_time, width=20)
        sorted_records = sorted(reader, key=lambda x: x.getField(clone_field))
        printMessage('Done', start_time=start_time, end=True, width=20)
        receptor_iter = groupby(sorted_records, lambda x: x.getField(clone_field))
        #receptor_iter = groupby(reader, lambda x: x.getField(clone_field))
    else:
        receptor_iter = ((x.sequence_id, [x]) for x in reader)

    # Initialize time and total count for progress bar
    start_time = time()
    total_count = countDbFile(db_file)
    rec_count = pass_count = fail_count = 0
    # Iterate over rows
    for key, records in receptor_iter:
        # Print progress
        printProgress(rec_count, total_count, 0.05, start_time)

        # Define iteration variables
        records = list(records)
        rec_log = OrderedDict([('ID', key)])
        rec_count += len(records)

        # Build germline for records
        if len(records) == 1:
            germ_log, germlines, genes = buildGermline(records[0], references, seq_field=seq_field, v_field=v_field,
                                                       d_field=d_field, j_field=j_field)
        else:
            germ_log, germlines, genes = buildClonalGermline(records, references, seq_field=seq_field, v_field=v_field,
                                                             d_field=d_field, j_field=j_field)
        rec_log.update(germ_log)

        # Write row to pass or fail file
        if germlines is not None:
            pass_count += len(records)

            # Add germlines to Receptor record
            annotations = {}
            if 'full' in germ_types:  annotations[germline_fields['full']] = germlines['full']
            if 'dmask' in germ_types:  annotations[germline_fields['dmask']] = germlines['dmask']
            if 'vonly' in germ_types:  annotations[germline_fields['vonly']] = germlines['vonly']
            if 'regions' in germ_types:  annotations[germline_fields['regions']] = germlines['regions']
            if cloned:
                annotations[germline_fields['v']] = genes['v']
                annotations[germline_fields['d']] = genes['d']
                annotations[germline_fields['j']] = genes['j']

            # Write records
            for r in records:
                r.setDict(annotations)
                pass_writer.writeReceptor(r)
        else:
            fail_count += len(records)
            if fail_writer is not None:
                for r in records: fail_writer.writeReceptor(r)

        # Write log
        printLog(rec_log, handle=log_handle)

    # Print log
    printProgress(rec_count, total_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    log['END'] = 'CreateGermlines'
    printLog(log)

    # Close file handles
    db_handle.close()
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
    parser_parent = getCommonArgParser(format=True)

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
                             germline with D segment masked, or germline for V segment only.''')
    parser.add_argument('--cloned', action='store_true', dest='cloned',
                        help='''Specify to create only one germline per clone. Note, if allele 
                             calls are ambiguous within a clonal group, this will place the germline call 
                             used for the entire clone within the
                             GERMLINE_V_CALL, GERMLINE_D_CALL and GERMLINE_J_CALL fields.''')
    parser.add_argument('--sf', action='store', dest='seq_field', default=None,
                        help='Field containing the alinged sequence.')
    parser.add_argument('--vf', action='store', dest='v_field', default=None,
                        help='Field containing the germline V segment call.')
    parser.add_argument('--df', action='store', dest='d_field', default=None,
                        help='Field containing the germline D segment call.')
    parser.add_argument('--jf', action='store', dest='j_field', default=None,
                        help='Field containing the germline J segment call.')

    return parser


if __name__ == '__main__':
    """
    Parses command line arguments and calls main
    """

    # Parse command line arguments
    parser = getArgParser()
    checkArgs(parser)
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)
    del args_dict['db_files']

    # Set default fields if not specified.
    default_fields = {'seq_field': default_seq_field,
                      'v_field': default_v_field,
                      'd_field': default_d_field,
                      'j_field': default_j_field}

    # Default Change-O fields
    if args_dict['format'] == 'changeo':
        for f in default_fields:
            if args_dict[f] is None:  args_dict[f] = default_fields[f]
            else: args_dict[f] = args_dict[f].upper()

    # Default AIRR fields
    if args_dict['format'] == 'airr':
        for f in default_fields:
            if args_dict[f] is None:  args_dict[f] = ChangeoSchema.asAIRR(default_fields[f])
            else: args_dict[f] = args_dict[f].lower()

    # Call main for each input file
    for f in args.__dict__['db_files']:
        args_dict['db_file'] = f
        createGermlines(**args_dict)
