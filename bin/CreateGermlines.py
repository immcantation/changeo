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
from textwrap import dedent
from time import time

# Presto and change imports
from presto.Defaults import default_out_args
from presto.IO import getOutputHandle, printLog, printProgress
from changeo.Defaults import default_format, default_v_field, default_d_field, default_j_field
from changeo.Commandline import CommonHelpFormatter, checkArgs, getCommonArgParser, parseCommonArgs
from changeo.IO import countDbFile, getDbFields, readRepo
from changeo.Receptor import allele_regex, parseAllele, Receptor
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
    vgene = parseAllele(receptor.getField(v_field), allele_regex, 'first')

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
    dgene = parseAllele(receptor.getField(d_field), allele_regex, 'first')

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
    jgene = parseAllele(receptor.getField(j_field), allele_regex, 'first')

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
                  d_field=default_d_field, j_field=default_j_field,
                  germ_types=default_germ_types):
    """
    Join gapped germline sequences aligned with sample sequences

    Arguments:
      receptor : Receptor object
      references : dictionary of IMGT gapped germline sequences
      seq_field : field in which to look for sequence
      v_field : field in which to look for V call
      d_field : field in which to look for V call
      j_field : field in which to look for V call
      germ_types : types of germline sequences to be output
                   (full germline, D-region masked, only V-region germline)

    Returns:
      tuple : log dictionary, dictionary of {germline_type: germline_sequence}
    """
    germlines = {'full': '', 'dmask': '', 'vonly': '', 'regions': ''}

    # Define log
    log = OrderedDict()
    log['ID'] = receptor.sequence_id

    # Build V segment germline sequence
    vgene, germ_vseq = getVGermline(receptor, references, v_field=v_field)
    log['V_CALL'] = vgene
    if germ_vseq is None:
        log['ERROR'] = 'Allele %s in not in the provided germline database.' % vgene
        return log, germlines

    # Build D segment germline sequence
    dgene, germ_dseq = getDGermline(receptor, references, d_field=d_field)
    log['D_CALL'] = dgene
    if germ_dseq is None:
        log['ERROR'] = 'Allele %s in not in the provided germline database.' % vgene
        return log, germlines

    # Build J segment germline sequence
    jgene, germ_jseq = getJGermline(receptor, references, j_field=j_field)
    log['J_CALL'] = jgene
    if germ_jseq is None:
        log['ERROR'] = 'Allele %s in not in the provided germline database.' % vgene
        return log, germlines

    # Stitch complete germlines
    germ_seq = stitchVDJ(receptor, germ_vseq, germ_dseq, germ_jseq)
    regions = stitchRegions(receptor, germ_vseq, germ_dseq, germ_jseq)

    # Define return germlines
    germlines['full'] = germ_seq
    germlines['regions'] = regions

    if 'dmask' in germ_types:
        germlines['dmask'] = germ_seq[:len(germ_vseq)] + \
                             'N' * (len(germ_seq) - len(germ_vseq) - len(germ_jseq)) + \
                             germ_seq[-len(germ_jseq):]
    if 'vonly' in germ_types:
        germlines['vonly'] = germ_vseq

    # Convert to uppercase
    for k, v in germlines.items():  germlines[k] = v.upper()

    # Update log
    log['SEQUENCE'] = receptor.getField(seq_field)
    log['GERMLINE'] = germ_seq
    log['REGIONS'] = regions

    # Check that input and germline sequence match
    if len(receptor.getField(seq_field)) == 0:
        log['ERROR'] = 'Sequence is missing from the %s field' % seq_field
    elif len(germlines['full']) != len(receptor.getField(seq_field)):
        log['ERROR'] = 'Germline sequence is %d nucleotides longer than input sequence' % \
                              (len(germlines['full']) - len(receptor.getField(seq_field)))

    return log, germlines


def joinGermline(align, references, seq_field=default_seq_field, v_field=default_v_field,
                 d_field=default_d_field, j_field=default_j_field,
                 germ_types=default_germ_types):
    """
    Join gapped germline sequences aligned with sample sequences
    
    Arguments:
    align = iterable yielding dictionaries of sample sequence data
    references = dictionary of IMGT gapped germline sequences
    seq_field = field in which to look for sequence
    v_field = field in which to look for V call
    d_field = field in which to look for V call
    j_field = field in which to look for V call
    germ_types = types of germline sequences to be output
                 (full germline, D-region masked, only V-region germline)

    Returns:
    dictionary of germline_type: germline_sequence
    """
    j_field = 'J_CALL'
    germlines = {'full': '', 'dmask': '', 'vonly': '', 'regions': ''}

    # Define log
    result_log = OrderedDict()
    result_log['ID'] = align['SEQUENCE_ID']

    # Set mode for region definitions
    try:
        int(align['P3V_LENGTH'])
        int(align['N1_LENGTH'])
        int(align['P5D_LENGTH'])
        int(align['P3D_LENGTH'])
        int(align['N2_LENGTH'])
        int(align['P5J_LENGTH'])
    except:
        regions_style = 'IgBLAST'
    else:
        regions_style = 'IMGT'

    # Find germline V-region gene
    vgene = parseAllele(align[v_field], allele_regex, 'first')

    # Build V-region germline
    if vgene is not None:
        result_log['V_CALL'] = vgene
        if vgene in references:
            vseq = references[vgene]
            # Germline start
            try: vstart = int(align['V_GERM_START_IMGT']) - 1
            except (TypeError, ValueError): vstart = 0
            # Germline length
            try: vlen = int(align['V_GERM_LENGTH_IMGT'])
            except (TypeError, ValueError): vlen = 0
            # TODO:  not sure what this line is doing here. it no make no sense.
            vpad = vlen - len(vseq[vstart:])
            if vpad < 0: vpad = 0
            germ_vseq = vseq[vstart:(vstart + vlen)] + ('N' * vpad)
        else:
            result_log['ERROR'] = 'Germline %s not in repertoire' % vgene
            return result_log, germlines
    else:
        result_log['V_CALL'] = None
        try: vlen = int(align['V_GERM_LENGTH_IMGT'])
        except (TypeError, ValueError): vlen = 0
        germ_vseq = 'N' * vlen

    # Find germline D-region gene
    dgene = parseAllele(align[d_field], allele_regex, 'first')

    # Build D-region germline
    if dgene is not None:
        result_log['D_CALL'] = dgene
        if dgene in references:
            dseq = references[dgene]
            # Germline start
            try: dstart = int(align['D_GERM_START']) - 1
            except (TypeError, ValueError): dstart = 0
            # Germline length
            try: dlen = int(align['D_GERM_LENGTH'])
            except (TypeError, ValueError): dlen = 0
            germ_dseq = dseq[dstart:(dstart + dlen)]
        else:
            result_log['ERROR'] = 'Germline %s not in repertoire' % dgene
            return result_log, germlines
    else:
        result_log['D_CALL'] = None
        germ_dseq = ''

    # Find germline J-region gene
    jgene = parseAllele(align[j_field], allele_regex, 'first')

    # Build D-region germline
    if jgene is not None:
        result_log['J_CALL'] = jgene
        if jgene in references:
            jseq = references[jgene]
            # Germline start
            try: jstart = int(align['J_GERM_START']) - 1
            except (TypeError, ValueError): jstart = 0
            # Germline length
            try: jlen = int(align['J_GERM_LENGTH'])
            except (TypeError, ValueError): jlen = 0
            # TODO:  not sure what this line is doing either
            jpad = jlen - len(jseq[jstart:])
            if jpad < 0: jpad = 0
            germ_jseq = jseq[jstart:(jstart + jlen)] + ('N' * jpad)
        else:
            result_log['ERROR'] = 'Germline %s not in repertoire' % jgene
            return result_log, germlines
    else:
        result_log['J_CALL'] = None
        try: jlen = int(align['J_GERM_LENGTH'])
        except (TypeError, ValueError): jlen = 0
        germ_jseq = 'N' * jlen

    # Assemble pieces starting with V-region
    germ_seq = germ_vseq
    regions = 'V' * len(germ_vseq)

    try:
        np1_len = int(align['NP1_LENGTH'])
    except (TypeError, ValueError):
        np1_len = 0

    # NP nucleotide additions after V
    if regions_style == 'IMGT':
        # P nucleotide additions
        try:
            p3v_len = int(align['P3V_LENGTH'])
        except (TypeError, ValueError):
            p3v_len = 0
        if p3v_len < 0:
            result_log['ERROR'] = 'P3V_LENGTH is negative'
            return result_log, germlines

        regions += 'P' * p3v_len

        # N1 nucleotide additions
        try:
            n1_len = int(align['N1_LENGTH'])
        except (TypeError, ValueError):
            n1_len = 0
        if n1_len < 0:
            result_log['ERROR'] = 'N1_LENGTH is negative'
            return result_log, germlines

        regions += 'N' * n1_len

        # P nucleotide additions before D
        try: p5d_len = int(align['P5D_LENGTH'])
        except (TypeError, ValueError): p5d_len = 0
        if p5d_len < 0:
            result_log['ERROR'] = 'P5D_LENGTH is negative'
            return result_log, germlines

        regions += 'P' * p5d_len
    else:
        # IgBLAST style
        # PNP nucleotide additions after V
        if np1_len < 0:
            result_log['ERROR'] = 'NP1_LENGTH is negative'
            return result_log, germlines

        regions += 'N' * np1_len

    germ_seq += 'N' * np1_len

    # Add D-region
    germ_seq += germ_dseq
    regions += 'D' * len(germ_dseq)

    #print 'VD>', germ_seq, '\nVD>', regions

    try:
        np2_len = int(align['NP2_LENGTH'])
    except (TypeError, ValueError):
        np2_len = 0

    # NP nucleotide additions before J
    if regions_style == 'IMGT':
        # P nucleotide additions
        try:
            p3d_len = int(align['P3D_LENGTH'])
        except (TypeError, ValueError):
            p3d_len = 0
        if p3d_len < 0:
            result_log['ERROR'] = 'P3D_LENGTH is negative'
            return result_log, germlines

        regions += 'P' * p3d_len

        # N2 nucleotide additions
        try:
            n2_len = int(align['N2_LENGTH'])
        except (TypeError, ValueError):
            n2_len = 0
        if n2_len < 0:
            result_log['ERROR'] = 'N2_LENGTH is negative'
            return result_log, germlines

        regions += 'N' * n2_len

        # P nucleotide additions
        try:
            p5j_len = int(align['P5J_LENGTH'])
        except (TypeError, ValueError):
            p5j_len = 0
        if p5j_len < 0:
            result_log['ERROR'] = 'P5J_LENGTH is negative'
            return result_log, germlines

        regions += 'P' * p5j_len
    else:
        # IgBLAST style
        # NP nucleotide additions
        if np2_len < 0:
            result_log['ERROR'] = 'NP2_LENGTH is negative'
            return result_log, germlines

        regions += 'N' * np2_len

    germ_seq += 'N' * np2_len

    # Add J-region
    germ_seq += germ_jseq
    regions += 'J' * len(germ_jseq)

    #print('\nREGIONS>',regions,'\n')

    # Define return germlines
    germlines['full'] = germ_seq
    germlines['regions'] = regions

    if 'dmask' in germ_types:
        germlines['dmask'] = germ_seq[:len(germ_vseq)] + \
                             'N' * (len(germ_seq) - len(germ_vseq) - len(germ_jseq)) + \
                             germ_seq[-len(germ_jseq):]
    if 'vonly' in germ_types:
        germlines['vonly'] = germ_vseq

    # Check that input and germline sequence match
    if len(align[seq_field]) == 0:
        result_log['ERROR'] = 'Sequence is missing from %s column' % seq_field
    elif len(germlines['full']) != len(align[seq_field]):
        result_log['ERROR'] = 'Germline sequence is %d nucleotides longer than input sequence' % \
                              (len(germlines['full']) - len(align[seq_field]))

    # Convert to uppercase
    for k, v in germlines.items():  germlines[k] = v.upper()

    return result_log, germlines


def createGermlines(db_file, repo, seq_field=default_seq_field, v_field=default_v_field,
                    d_field=default_d_field, j_field=default_j_field,
                    germ_types=default_germ_types, cloned=False,
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
      germ_types : list of germline sequence types to be output from the set of 'full', 'dmask', 'vonly', 'regions'
      cloned : if True build germlines by clone, otherwise build individual germlines
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
        out_fields = getDbFields(db_file, add=list(germline_fields.values()), reader=reader)
    else:
        sys.exit('Error:  Invalid format %s' % format)

    # Get repertoire and open Db reader
    references = readRepo(repo)
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

    # Initialize time and total count for progress bar
    start_time = time()
    rec_count = countDbFile(db_file)
    clone_count = pass_count = fail_count = 0
    clone_id = None
    receptor_list = []
    # Iterate over rows
    for i, rec in enumerate(reader):
        # Print progress
        printProgress(i, rec_count, 0.05, start_time)

        # TODO: if cloned: (1) check clone id, (2) update receptor list, (3) make clonal germline, (4) write all records in clone
        # TODO: add check for sorted clones earlier. when counting clones.
        # TODO: whole cloned/not-cloned loops could be two processing/writer function. with log output embedded.
        # TODO: or just a writer than takes a Receptor list, germlines, and log
        result_log, germlines = buildGermline(rec, references, seq_field=seq_field, v_field=v_field,
                                              d_field=d_field, j_field=j_field, germ_types=germ_types)

        # Add germlines to Receptor record
        # TODO: this should probably pass through the schemas instead so the output fields are correct
        annotations = {}
        if 'full' in germ_types:  annotations[germline_fields['full']] = germlines['full']
        if 'dmask' in germ_types:  annotations[germline_fields['dmask']] = germlines['dmask']
        if 'vonly' in germ_types:  annotations[germline_fields['vonly']] = germlines['vonly']
        if 'regions' in germ_types:  annotations[germline_fields['regions']] = germlines['regions']
        rec.setDict(annotations)

        # Write row to pass or fail file
        if 'ERROR' in result_log:
            fail_count += 1
            if fail_writer is not None: fail_writer.writeReceptor(rec)
        else:
            pass_count += 1
            pass_writer.writeReceptor(rec)
        printLog(result_log, handle=log_handle)

    # Print log
    printProgress(i + 1, rec_count, 0.05, start_time)
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


def buildClonalGermline(receptors, references, seq_field=default_seq_field,
                        v_field=default_v_field, d_field=default_d_field, j_field=default_j_field,
                        germ_types=default_germ_types):
    """
    Determine consensus clone sequence and create germline for clone

    Arguments:
      receptors : list of Receptor objects
      references : dictionary of IMGT gapped germline sequences
      seq_field : field in which to look for sequence
      v_field : field in which to look for V call
      d_field : field in which to look for D call
      j_field : field in which to look for J call
      germ_types : list of germline sequence types to be output from the set of 'full', 'dmask', 'vonly', 'regions'
      out_args = arguments for output preferences

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
        v = rec.getField(v_field)
        v_dict[v] = v_dict.get(v, 0) + 1
        j = rec.getField(j_field)
        j_dict[j] = j_dict.get(j, 0) + 1
        seq_len = len(rec.getField(seq_field))
        if seq_len > max_length:
            max_length = seq_len

    # Consensus V and J having most observations
    v_cons = [k for k in list(v_dict.keys()) if v_dict[k] == max(v_dict.values())]
    j_cons = [k for k in list(j_dict.keys()) if j_dict[k] == max(j_dict.values())]

    # Consensus sequence(s) with consensus V/J calls and longest sequence
    cons = [x for x in receptors if x.getField(v_field) in v_cons and \
                                    x.getField(j_field) in j_cons and \
                                    len(x.getField(seq_field)) == max_length]
    # Consensus sequence(s) with consensus V/J calls but not the longest sequence
    if not cons:
        cons = [x for x in receptors if x.getField(v_field) in v_cons and \
                                        x.getField(j_field) in j_cons]

    # Return without germline if no sequence has both consensus V and J call
    if not cons:
        log['V_CALL'] = ','.join(v_cons)
        log['J_CALL'] = ','.join(j_cons)
        log['ERROR'] = 'No sequence found with both consensus V and J calls.'
        return log, None, None

    # Assign consensus Receptor
    cons = cons[0]
    genes = {'v': cons.getField(v_field),
             'd': cons.getField(d_field),
             'j': cons.getField(j_field)}

    # Pad end of consensus sequence with gaps to make it the max length
    gap_length = max_length - cons.getField(seq_field)
    if gap_length > 0:
        cons.j_germ_length = int(cons.j_germ_length or 0) + gap_length
        cons.setField(seq_field, cons.getSeq(seq_field) + '-' * gap_length)

    # Update lengths padded to longest sequence in clone
    for rec in receptors:
        x = max_length - rec.getField(seq_field)
        rec.j_germ_length = int(rec.j_germ_length or 0) + x
        rec.setField(seq_field, rec.getSeq(seq_field) + '-' * x)

    # Stitch consensus germline
    join_log, germlines = buildGermline(cons, references, seq_field=seq_field, v_field=v_field,
                                        d_field=d_field, j_field=j_field, germ_types=germ_types)

    # Update log
    log['CONSENSUS'] = cons.sequence_id
    log['V_CALL'] = genes['v']
    log['D_CALL'] = genes['d']
    log['J_CALL'] = genes['j']
    log.update(join_log)

    # Return log
    return log, germlines, genes


def makeCloneGermline(clone, clone_dict, references, germ_types, v_field,
                      seq_field, counts, writers, out_args):
    """
    Determine consensus clone sequence and create germline for clone

    Arguments:
    clone = clone ID
    clone_dict = iterable yielding dictionaries of sequence data from clone
    references = dictionary of IMGT gapped germline sequences
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
    for val in clone_dict.values():
        v = val[v_field]
        v_dict[v] = v_dict.get(v, 0) + 1
        j = val[j_field]
        j_dict[j] = j_dict.get(j, 0) + 1
        if len(val[seq_field]) > max_length:
            max_length = len(val[seq_field])

    # Consensus V and J having most observations
    v_cons = [k for k in list(v_dict.keys()) if v_dict[k] == max(v_dict.values())]
    j_cons = [k for k in list(j_dict.keys()) if j_dict[k] == max(j_dict.values())]

    # Consensus sequence(s) with consensus V/J calls and longest sequence
    cons = [val for val in list(clone_dict.values()) \
            if val.get(v_field, '') in v_cons and \
            val.get(j_field, '') in j_cons and \
            len(val[seq_field]) == max_length]

    # Sequence(s) with consensus V/J are not longest
    if not cons:
        # Sequence(s) with consensus V/J (not longest)
        cons = [val for val in list(clone_dict.values()) \
                if val.get(v_field, '') in v_cons and val.get(j_field, '') in j_cons]

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
            cons[seq_field] += '.' * (max_length - len(cons[seq_field]))
            result_log, germlines = joinGermline(cons, references, seq_field=seq_field, v_field=v_field,
                                                 germ_types=germ_types)
            result_log['ID'] = clone
            result_log['CONSENSUS'] = cons['SEQUENCE_ID']
    else:
        cons = cons[0]
        result_log, germlines = joinGermline(cons, references, seq_field=seq_field, v_field=v_field,
                                             germ_types=germ_types)
        result_log['ID'] = clone
        result_log['CONSENSUS'] = cons['SEQUENCE_ID']

    # Write sequences of clone
    for val in clone_dict.values():
        if 'ERROR' not in result_log:
            # Update lengths padded to longest sequence in clone
            val['J_GERM_LENGTH'] = str(int(val['J_GERM_LENGTH'] or 0) + max_length - len(val[seq_field]))
            val[seq_field] += '.' * (max_length - len(val[seq_field]))

            # Add column(s) to tab-delimited database file
            if 'full' in germ_types: val['GERMLINE_' + seq_type] = germlines['full']
            if 'dmask' in germ_types: val['GERMLINE_' + seq_type + '_D_MASK'] = germlines['dmask']
            if 'vonly' in germ_types: val['GERMLINE_' + seq_type + '_V_REGION'] = germlines['vonly']
            if 'regions' in germ_types: val['GERMLINE_REGIONS'] = germlines['regions']

            # Add field
            val['GERMLINE_V_CALL'] = result_log['V_CALL']
            val['GERMLINE_D_CALL'] = result_log['D_CALL']
            val['GERMLINE_J_CALL'] = result_log['J_CALL']
            
            result_log['SEQUENCE'] = cons[seq_field]
            result_log['GERMLINE'] = germlines['full']
            result_log['REGIONS'] = germlines['regions']

            # Write to pass file
            counts['pass'] += 1
            writers['pass'].writeDict(val)
        else:
            # Write to fail file
            counts['fail'] += 1
            if writers['fail'] is not None:
                writers['fail'].writeDict(val)
    # Return log
    return result_log


def assembleCloneGermline(db_file, repo, seq_field=default_seq_field, v_field=default_v_field,
                          d_field=default_d_field, j_field=default_j_field,
                          germ_types=default_germ_types,
                          format=default_format, out_args=default_out_args):
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
    log['GERM_TYPES'] = germ_types if isinstance(germ_types, str) else ','.join(germ_types)
    log['CLONED'] = 'True'
    log['V_FIELD'] = v_field
    log['SEQ_FIELD'] = seq_field
    printLog(log)

    # Get repertoire and open Db reader
    references = readRepo(repo)
    db_handle = open(db_file, 'rt')
    reader = ChangeoReader(db_handle, receptor=False)

    # Exit if V call field does not exist in reader
    if v_field not in reader.fields:
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
    if 'regions' in germ_types: add_fields += ['GERMLINE_REGIONS']

    add_fields += ['GERMLINE_V_CALL']
    add_fields += ['GERMLINE_D_CALL']
    add_fields += ['GERMLINE_J_CALL']
    
    # Create output file handle and Db writer
    writers = {}
    pass_handle = getOutputHandle(db_file,
                                  out_label='germ-pass',
                                  out_dir=out_args['out_dir'],
                                  out_name=out_args['out_name'],
                                  out_type='tsv')
    out_fields = getDbFields(db_file, add=add_fields)
    writers['pass'] = ChangeoWriter(pass_handle, fields=out_fields)

    if out_args['failed']:
        fail_handle = getOutputHandle(db_file,
                                      out_label='germ-fail',
                                      out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'],
                                      out_type='tsv')
        writers['fail'] = ChangeoWriter(fail_handle, fields=out_fields)
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
    for i, row in enumerate(reader):
        # Print progress
        printProgress(i, rec_count, 0.05, start_time)

        # Clone isn't over yet
        if row.get('CLONE', '') == clone:
            clone_dict[i] = row
        # Clone just finished
        elif clone_dict:
            clone_count += 1
            result_log = makeCloneGermline(clone, clone_dict, references, germ_types,
                                           v_field, seq_field, counts, writers, out_args)
            printLog(result_log, handle=log_handle)
            # Now deal with current row (first of next clone)
            clone = row['CLONE']
            clone_dict = OrderedDict([(i, row)])
        # Last case is only for first row of file
        else:
            clone = row['CLONE']
            clone_dict = OrderedDict([(i, row)])

    clone_count += 1
    result_log = makeCloneGermline(clone, clone_dict, references, germ_types, v_field,
                                   seq_field, counts, writers, out_args)
    printLog(result_log, handle=log_handle)

    # Print log
    printProgress(i + 1, rec_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['CLONES'] = clone_count
    log['RECORDS'] = rec_count
    log['PASS'] = counts['pass']
    log['FAIL'] = counts['fail']
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
                             germline with D-region masked, or germline for V region only.''')
    parser.add_argument('--cloned', action='store_true', dest='cloned',
                        help='''Specify to create only one germline per clone. Assumes input file is
                             sorted by clone column, and will not yield correct results if the data
                             is unsorted. Note, if allele calls are ambiguous within a clonal group,
                             this will place the germline call used for the entire clone within the
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


if __name__ == "__main__":
    """
    Parses command line arguments and calls main
    """

    # Parse command line arguments
    parser = getArgParser()
    checkArgs(parser)
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)
    del args_dict['db_files']
    del args_dict['cloned']

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


    for f in args.__dict__['db_files']:
        args_dict['db_file'] = f
        if args.__dict__['cloned']:
            assembleCloneGermline(**args_dict)
        else:
            createGermlines(**args_dict)
