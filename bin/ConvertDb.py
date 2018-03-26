#!/usr/bin/env python3
"""
Parses tab delimited database files
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from changeo import __version__, __date__

# Imports
import csv
import os
import re
import sys
from argparse import ArgumentParser
from collections import OrderedDict
from itertools import chain

from textwrap import dedent
from time import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

# Presto and changeo imports
from presto.Annotation import flattenAnnotation
from presto.IO import getOutputHandle, printLog, printProgress
from changeo.Defaults import default_csv_size, default_format, default_out_args
from changeo.Commandline import CommonHelpFormatter, checkArgs, getCommonArgParser, parseCommonArgs
from changeo.IO import countDbFile, TSVReader
from changeo.Parsers import AIRRReader, AIRRWriter, ChangeoReader, ChangeoWriter
from changeo.Receptor import c_gene_regex, parseAllele, AIRRSchema, ChangeoSchema

# System settings
csv.field_size_limit(default_csv_size)

# Defaults
default_id_field = 'SEQUENCE_ID'
default_seq_field = 'SEQUENCE_IMGT'
default_germ_field = 'GERMLINE_IMGT_D_MASK'
default_db_xref = 'IMGT/GENE-DB'
default_molecule='mRNA'
default_product='immunoglobulin heavy chain'

# TODO:  convert SQL-ish operations to modify_func() as per ParseHeaders

def buildSeqRecord(db_record, id_field, seq_field, meta_fields=None):
    """
    Parses a database record into a SeqRecord

    Arguments: 
    db_record = a dictionary containing a database record
    id_field = the field containing identifiers
    seq_field = the field containing sequences
    meta_fields = a list of fields to add to sequence annotations

    Returns: 
    a SeqRecord
    """
    # Return None if ID or sequence fields are empty
    if not db_record[id_field] or not db_record[seq_field]:
        return None
    
    # Create description string
    desc_dict = OrderedDict([('ID', db_record[id_field])])
    if meta_fields is not None:
        desc_dict.update([(f, db_record[f]) for f in meta_fields if f in db_record]) 
    desc_str = flattenAnnotation(desc_dict)
    
    # Create SeqRecord
    seq_record = SeqRecord(Seq(db_record[seq_field], IUPAC.ambiguous_dna),
                           id=desc_str, name=desc_str, description='')
        
    return seq_record


def convertDbAIRR(db_file, out_file=None, out_args=default_out_args):
    """
    Converts a Change-O formatted file into an AIRR formatted file

    Arguments:
      db_file : the database file name.
      out_file : output file name. Automatically generated from the input file if None.
      out_args : common output argument dictionary from parseCommonArgs.

    Returns:
     str : output file name
    """
    log = OrderedDict()
    log['START'] = 'ConvertDb'
    log['COMMAND'] = 'airr'
    log['FILE'] = os.path.basename(db_file)
    printLog(log)

    # Open input
    db_handle = open(db_file, 'rt')
    db_iter = ChangeoReader(db_handle, receptor=True)

    # Set fields
    fields_delete = ['v_seq_length', 'v_germ_length_vdj',
                     'd_seq_length', 'd_germ_length',
                     'j_seq_length', 'j_germ_length']
    fields_add = ['v_sequence_start', 'v_germline_start',
                  'd_sequence_start', 'd_germline_start',
                  'j_sequence_start', 'j_germline_start',
                  'v_sequence_end', 'v_germline_end',
                  'd_sequence_end', 'd_germline_end',
                  'j_sequence_end', 'j_germline_end']
    fields = [ChangeoSchema.asAIRR(x) for x in db_iter.fields]
    fields = [x for x in fields if x not in fields_delete]
    fields.extend([x for x in fields_add if x not in fields])

    # Open output writer
    if out_file is not None:
        pass_handle = open(out_file, 'w')
    else:
        pass_handle = getOutputHandle(db_file, out_label='airr', out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'], out_type='tsv')
    pass_writer = AIRRWriter(pass_handle, fields=fields)

    # Count records
    result_count = countDbFile(db_file)

    # Iterate over records
    start_time = time()
    rec_count = pass_count = fail_count = 0
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time)
        rec_count += 1
        # Write records
        pass_writer.writeReceptor(rec)

    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['END'] = 'ConvertDb'
    printLog(log)

    # Close file handles
    pass_handle.close()
    db_handle.close()

    return pass_handle.name


def convertDbChangeo(db_file, out_file=None, out_args=default_out_args):
    """
    Converts a AIRR formatted file into an Change-O formatted file

    Arguments:
      db_file = the database file name.
      out_file : output file name. Automatically generated from the input file if None.
      out_args = common output argument dictionary from parseCommonArgs.

    Returns:
      str : output file name.
    """
    log = OrderedDict()
    log['START'] = 'ConvertDb'
    log['COMMAND'] = 'changeo'
    log['FILE'] = os.path.basename(db_file)
    printLog(log)

    # Open input
    db_handle = open(db_file, 'rt')
    db_iter = AIRRReader(db_handle, receptor=True)

    # Set fields
    fields_delete = ['V_SEQ_END', 'V_GERM_END_VDJ',
                     'D_SEQ_END', 'D_GERM_END',
                     'J_SEQ_END', 'J_GERM_END']
    fields_add = ['V_SEQ_START', 'V_GERM_START_VDJ',
                  'D_SEQ_START', 'D_GERM_START',
                  'J_SEQ_START', 'J_GERM_START',
                  'V_SEQ_LENGTH', 'V_GERM_LENGTH_VDJ',
                  'D_SEQ_LENGTH', 'D_GERM_LENGTH',
                  'J_SEQ_LENGTH', 'J_GERM_LENGTH']
    fields = [AIRRSchema.asChangeo(x) for x in db_iter.fields]
    fields = [x for x in fields if x not in fields_delete]
    fields.extend([x for x in fields_add if x not in fields])

    # Open output writer
    if out_file is not None:
        pass_handle = open(out_file, 'w')
    else:
        pass_handle = getOutputHandle(db_file, out_label='changeo', out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'], out_type='tab')
    pass_writer = ChangeoWriter(pass_handle, fields=fields)

    # Count records
    result_count = countDbFile(db_file)

    # Iterate over records
    start_time = time()
    rec_count = pass_count = fail_count = 0
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time)
        rec_count += 1
        # Write records
        pass_writer.writeReceptor(rec)

    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['END'] = 'ConvertDb'
    printLog(log)

    # Close file handles
    pass_handle.close()
    db_handle.close()

    return pass_handle.name


# TODO:  SHOULD ALLOW FOR UNSORTED CLUSTER COLUMN
# TODO:  SHOULD ALLOW FOR GROUPING FIELDS
def convertDbBaseline(db_file, id_field=default_id_field, seq_field=default_seq_field,
                      germ_field=default_germ_field, cluster_field=None,
                      meta_fields=None, out_file=None, out_args=default_out_args):
    """
    Builds fasta files from database records

    Arguments: 
      db_file : the database file name.
      id_field : the field containing identifiers.
      seq_field : the field containing sample sequences.
      germ_field : the field containing germline sequences.
      cluster_field : the field containing clonal groupings;
                    if None write the germline for each record.
      meta_fields : a list of fields to add to sequence annotations.
      out_file : output file name. Automatically generated from the input file if None.
      out_args : common output argument dictionary from parseCommonArgs.
                    
    Returns: 
     str : output file name
    """
    log = OrderedDict()
    log['START'] = 'ConvertDb'
    log['COMMAND'] = 'fasta'
    log['FILE'] = os.path.basename(db_file)
    log['ID_FIELD'] = id_field
    log['SEQ_FIELD'] = seq_field
    log['GERM_FIELD'] = germ_field
    log['CLUSTER_FIELD'] = cluster_field
    if meta_fields is not None:  log['META_FIELDS'] = ','.join(meta_fields)
    printLog(log)
    
    # Open input
    db_handle = open(db_file, 'rt')
    db_iter = TSVReader(db_handle)
    result_count = countDbFile(db_file)

    # Open output
    if out_file is not None:
        pass_handle = open(out_file, 'w')
    else:
        pass_handle = getOutputHandle(db_file, out_label='sequences', out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'], out_type='clip')
    # Iterate over records
    start_time = time()
    rec_count = germ_count = pass_count = fail_count = 0
    cluster_last = None
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time)
        rec_count += 1
        
        # Update cluster ID
        cluster = rec.get(cluster_field, None)
        
        # Get germline SeqRecord when needed
        if cluster_field is None:
            germ = buildSeqRecord(rec, id_field, germ_field, meta_fields)
            germ.id = '>' + germ.id
        elif cluster != cluster_last:
            germ = buildSeqRecord(rec, cluster_field, germ_field)
            germ.id = '>' + germ.id            
        else:
            germ = None

        # Get read SeqRecord
        seq = buildSeqRecord(rec, id_field, seq_field, meta_fields)
        
        # Write germline
        if germ is not None:
            germ_count += 1
            SeqIO.write(germ, pass_handle, 'fasta')
        
        # Write sequences
        if seq is not None:
            pass_count += 1
            SeqIO.write(seq, pass_handle, 'fasta')
        else:
            fail_count += 1
        
        # Set last cluster ID
        cluster_last = cluster
        
    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['GERMLINES'] = germ_count
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    log['END'] = 'ConvertDb'
    printLog(log)

    # Close file handles
    pass_handle.close()
    db_handle.close()

    return pass_handle.name


def convertDbFasta(db_file, id_field=default_id_field, seq_field=default_seq_field,
                 meta_fields=None, out_file=None, out_args=default_out_args):
    """
    Builds fasta files from database records

    Arguments: 
      db_file : the database file name.
      id_field : the field containing identifiers.
      seq_field : the field containing sequences.
      meta_fields : a list of fields to add to sequence annotations.
      out_file : output file name. Automatically generated from the input file if None.
      out_args : common output argument dictionary from parseCommonArgs.
                    
    Returns: 
      str : output file name.
    """
    log = OrderedDict()
    log['START'] = 'ConvertDb'
    log['COMMAND'] = 'fasta'
    log['FILE'] = os.path.basename(db_file)
    log['ID_FIELD'] = id_field
    log['SEQ_FIELD'] = seq_field
    if meta_fields is not None:  log['META_FIELDS'] = ','.join(meta_fields)
    printLog(log)
    
    # Open input
    out_type = 'fasta'
    db_handle = open(db_file, 'rt')
    db_iter = TSVReader(db_handle)
    result_count = countDbFile(db_file)

    # Open output
    if out_file is not None:
        pass_handle = open(out_file, 'w')
    else:
        pass_handle = getOutputHandle(db_file, out_label='sequences', out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'], out_type=out_type)

    # Iterate over records
    start_time = time()
    rec_count = pass_count = fail_count = 0
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time)
        rec_count += 1

        # Get SeqRecord
        seq = buildSeqRecord(rec, id_field, seq_field, meta_fields)

        # Write sequences
        if seq is not None:
            pass_count += 1
            SeqIO.write(seq, pass_handle, out_type)
        else:
            fail_count += 1
        
    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    log['END'] = 'ConvertDb'
    printLog(log)

    # Close file handles
    pass_handle.close()
    db_handle.close()

    return pass_handle.name


def makeGenbankFeatures(record, start=None, end=None, product=default_product,
                        inference=None, db_xref=default_db_xref,
                        cregion_field=None, keep_stop=False):
    """
    Creates a feature table for GenBank submissions

    Arguments:
      record : Receptor record.
      start : start position of the modified seqeuence in the input sequence. Used for feature position offsets.
      end : end position of the modified seqeuence in the input sequence. Used for feature position offsets.
      product : Product (protein) name.
      inference : Reference alignment tool.
      db_xref : Reference database name.
      cregion_field : column containing the C region gene call.
      keep_stop : if True retain records with junctions having stop codons.

    Returns:
      dict : dictionary defining GenBank features where the key is a tuple
             (start, end, feature key) and values are a list of
             tuples contain (qualifier key, qualifier value).
    """
    # .tbl file format
    #   Line 1, Column 1: Start location of feature
    #   Line 1, Column 2: Stop location of feature
    #   Line 1, Column 3: Feature key
    #   Line 2, Column 4: Qualifier key
    #   Line 2, Column 5: Qualifier value

    # Get genes and check for valid record
    v_gene = record.getVGene()
    d_gene = record.getDGene()
    j_gene = record.getJGene()
    if v_gene is None or j_gene is None:
        return None

    # Set position offsets if required
    start_trim = 0 if start is None else start
    end_trim = 0 if end is None else len(record.sequence_input) - end

    # Define return object
    result = OrderedDict()

    # C_region
    #     gene
    #     db_xref
    #     inference
    c_region_start = record.j_seq_end + 1 - start_trim
    c_region_length = len(record.sequence_input[(c_region_start + start_trim - 1):]) - end_trim
    if c_region_length > 0:
        if cregion_field is None:
            c_region = []
        else:
            c_call = record.getField(cregion_field)
            c_gene = parseAllele(c_call, c_gene_regex, action='first')
            c_region = [('gene', c_gene),
                        ('db_xref', '%s:%s' % (db_xref, c_gene))]

        # Assign C_region feature
        result[(c_region_start, '>%i' % (c_region_start + c_region_length - 1), 'C_region')] = c_region

        # Preserve J segment end position
        j_end = record.j_seq_end
    else:
        # Trim J segment end position
        j_end = record.j_seq_end + c_region_length

    # V_region
    variable_start = max(record.v_seq_start - start_trim, 1)
    variable_end = j_end - start_trim
    variable_region = []
    result[(variable_start, variable_end, 'V_region')] = variable_region

    # Product feature
    result[(variable_start, variable_end, 'misc_feature')] = [('note', '%s variable region' % product)]

    # V_segment
    #     gene (gene name)
    #     allele (allele only, without gene name, don't use if ambiguous)
    #     db_xref (database link)
    #     inference (reference alignment tool)
    v_segment = [('gene', v_gene),
                 ('allele', record.getVAlleleNumber()),
                 ('db_xref', '%s:%s' % (db_xref, v_gene))]
    if inference is not None:
        v_segment.append(('inference', 'COORDINATES:alignment:%s' % inference))
    result[(variable_start, record.v_seq_end - start_trim, 'V_segment')] = v_segment

    # D_segment
    #     gene
    #     allele
    #     db_xref
    #     inference
    if d_gene:
        d_segment = [('gene', d_gene),
                     ('allele', record.getDAlleleNumber()),
                     ('db_xref', '%s:%s' % (db_xref, d_gene))]
        if inference is not None:
            d_segment.append(('inference', 'COORDINATES:alignment:%s' % inference))
        result[(record.d_seq_start - start_trim, record.d_seq_end - start_trim, 'D_segment')] = d_segment

    # J_segment
    #     gene
    #     allele
    #     db_xref
    #     inference
    j_segment = [('gene', j_gene),
                 ('allele', record.getVAlleleNumber()),
                 ('db_xref', '%s:%s' % (db_xref, j_gene))]
    if inference is not None:
        j_segment.append(('inference', 'COORDINATES:alignment:%s' % inference))
    result[(record.j_seq_start - start_trim, j_end - start_trim, 'J_segment')] = j_segment

    # CDS
    #     codon_start (must indicate codon offset)
    #     function = JUNCTION
    #     inference
    if record.junction_start is not None and record.junction_end is not None:
        # Define junction boundaries
        junction_start = record.junction_start - start_trim
        junction_end = record.junction_end - start_trim

        # CDS record
        cds_start = '<%i' % junction_start
        cds_end = '>%i' % junction_end
        cds_record = [('function', 'JUNCTION')]
        if inference is not None:
            cds_record.append(('inference', 'COORDINATES:protein motif:%s' % inference))

        # Check for valid translation
        junction_seq = record.sequence_input[(junction_start - 1):junction_end]
        if len(junction_seq) % 3 > 0:  junction_seq = junction_seq + 'N' * (3 - len(junction_seq) % 3)
        junction_aa = junction_seq.translate()

        # Return invalid record upon junction stop codon
        if '*' in junction_aa and not keep_stop:
            return None
        elif '*' in junction_aa:
            cds_record.append(('note', '%s junction region' % product))
            result[(cds_start, cds_end, 'misc_feature')] = cds_record
        else:
            cds_record.append(('product', '%s junction region' % product))
            cds_record.append(('codon_start', 1))
            result[(cds_start, cds_end, 'CDS')] = cds_record

    return result


def makeGenbankSequence(record, name=None, label=None, organism=None, sex=None, isolate=None,
                        tissue=None, cell=None, molecule=default_molecule):
    """
    Creates a sequence for GenBank submissions

    Arguments:
      record : Receptor record.
      name : sequence identifier for the output sequence. If None,
             use the original sequence identifier.
      label : a string to use as a label for the ID. if None do not add a field label.
      organism : scientific name of the organism.
      sex : sex.
      isolate : sample identifier.
      tissue : tissue type.
      cell : cell type.
      molecule : source molecule (eg, "mRNA", "genomic DNA")

    Returns:
      dict : dictionary with {'record': SeqRecord,
                              'start': start position in raw sequence,
                              'end': end position in raw sequence}
    """
    # Replace gaps with N
    seq = str(record.sequence_input)
    seq = seq.replace('-', 'N').replace('.', 'N')

    # Strip leading and trailing Ns
    head_match = re.search('^N+', seq)
    tail_match = re.search('N+$', seq)
    seq_start = head_match.end() if head_match else 0
    seq_end = tail_match.start() if tail_match else len(seq)

    # Define ID
    if name is None:
        name = record.sequence_id.split(' ')[0]
    if label is not None:
        name = '%s=%s' % (label, name)
    if organism is not None:
        name = '%s [organism=%s]' % (name, organism)
    if sex is not None:
        name = '%s [sex=%s]' % (name, sex)
    if isolate is not None:
        name = '%s [isolate=%s]' % (name, isolate)
    if tissue is not None:
        name = '%s [tissue-type=%s]' % (name, tissue)
    if cell is not None:
        name = '%s [cell-type=%s]' % (name, cell)
    name = '%s [moltype=%s] [keyword=TLS; Targeted Locus Study; AIRR; MiAIRR:1.0]' % (name, molecule)

    # Return SeqRecord and positions
    record = SeqRecord(Seq(seq[seq_start:seq_end], IUPAC.ambiguous_dna), id=name,
                       name=name, description='')
    result = {'record': record, 'start': seq_start, 'end': seq_end}

    return result


def convertDbGenbank(db_file, inference=None, db_xref=None, organism=None, sex=None,
                     isolate=None, tissue=None, cell=None, molecule=default_molecule,
                     product=default_product, cregion_field=None, label=None,
                     keep_id=False, keep_stop=False,
                     format=default_format, out_file=None, out_args=default_out_args):
    """
    Builds a GenBank submission tbl file from records

    Arguments:
      db_file : the database file name.
      inference : reference alignment tool.
      db_xref : reference database link.
      organism : scientific name of the organism.
      sex : sex.
      isolate : sample identifier.
      tissue : tissue type.
      cell : cell type.
      molecule : source molecule (eg, "mRNA", "genomic DNA")
      product : Product (protein) name.
      cregion_field : column containing the C region gene call.
      label : a string to use as a label for the ID. if None do not add a field label.
      keep_id : if True use the original sequence ID for the output IDs
      keep_stop : if True retain records with junctions having stop codons.
      format : input and output format.
      out_file : output file name without extension. Automatically generated from the input file if None.
      out_args : common output argument dictionary from parseCommonArgs.

    Returns:
      tuple : the output (feature table, fasta) file names.
    """
    log = OrderedDict()
    log['START'] = 'ConvertDb'
    log['COMMAND'] = 'genbank'
    log['FILE'] = os.path.basename(db_file)
    printLog(log)

    # Format options
    if format == 'changeo':
        reader = ChangeoReader
    elif format == 'airr':
        reader = AIRRReader
    else:
        sys.exit('Error:  Invalid format %s' % format)

    # Open input
    db_handle = open(db_file, 'rt')
    db_iter = reader(db_handle)

    # Open output
    if out_file is not None:
        out_name, __ = os.path.splitext(out_file)
        fsa_handle = open('%s.fsa' % out_name, 'w')
        tbl_handle = open('%s.tbl' % out_name, 'w')
    else:
        fsa_handle = getOutputHandle(db_file, out_label='genbank', out_dir=out_args['out_dir'],
                                     out_name=out_args['out_name'], out_type='fsa')
        tbl_handle = getOutputHandle(db_file, out_label='genbank', out_dir=out_args['out_dir'],
                                     out_name=out_args['out_name'], out_type='tbl')

    # Count records
    result_count = countDbFile(db_file)

    # Define writer
    writer = csv.writer(tbl_handle, delimiter='\t', quoting=csv.QUOTE_NONE)

    # Iterate over records
    start_time = time()
    rec_count = pass_count = fail_count = 0
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time)
        rec_count += 1

        # Extract table dictionary
        name = None if keep_id else rec_count
        seq = makeGenbankSequence(rec, name=name, label=label, organism=organism, sex=sex, isolate=isolate,
                                  tissue=tissue, cell=cell, molecule=molecule)
        tbl = makeGenbankFeatures(rec, start=seq['start'], end=seq['end'], product=product,
                                  db_xref=db_xref, inference=inference, cregion_field=cregion_field,
                                  keep_stop=keep_stop)

        if tbl is not None:
            pass_count +=1
            # Write table
            writer.writerow(['>Features', seq['record'].id])
            for feature, qualifiers in tbl.items():
                writer.writerow(feature)
                if qualifiers:
                    for x in qualifiers:
                        writer.writerow(list(chain(['', '', ''], x)))

            # Write sequence
            SeqIO.write(seq['record'], fsa_handle, 'fasta')
        else:
            fail_count += 1

    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT_TBL'] = os.path.basename(tbl_handle.name)
    log['OUTPUT_FSA'] = os.path.basename(fsa_handle.name)
    log['RECORDS'] = rec_count
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    log['END'] = 'ConvertDb'
    printLog(log)

    # Close file handles
    tbl_handle.close()
    fsa_handle.close()
    db_handle.close()

    return (tbl_handle.name, fsa_handle.name)


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
                 airr
                     AIRR formatted database files.
                 changeo
                     Change-O formatted database files.
                 sequences
                     FASTA formatted sequences output from the subcommands fasta and clip.
                 genbank
                     feature tables and fasta files containing MiAIRR compliant input for tbl2asn.

             required fields:
                 SEQUENCE_ID, SEQUENCE_INPUT, JUNCTION, V_CALL, D_CALL, J_CALL, 
                 V_SEQ_START, V_SEQ_LENGTH, D_SEQ_START, D_SEQ_LENGTH, J_SEQ_START, J_SEQ_LENGTH
                 
             optional fields:
                 SEQUENCE_IMGT, GERMLINE_IMGT, GERMLINE_IMGT_D_MASK, CLONE, CREGION
                
             output fields:
                 None
             ''')
    
    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            formatter_class=CommonHelpFormatter, add_help=False)
    group_help = parser.add_argument_group('help')
    group_help.add_argument('--version', action='version',
                            version='%(prog)s:' + ' %s-%s' %(__version__, __date__))
    group_help.add_argument('-h', '--help', action='help', help='show this help message and exit')
    subparsers = parser.add_subparsers(title='subcommands', dest='command', metavar='',
                                       help='Database operation')
    # TODO:  This is a temporary fix for Python issue 9253
    subparsers.required = True

    # Define parent parsers
    default_parent = getCommonArgParser(failed=False, log=False, format=False)
    genbank_parent = getCommonArgParser(failed=False, log=False)

    # Subparser to convert changeo to AIRR files
    parser_airr = subparsers.add_parser('airr', parents=[default_parent],
                                       formatter_class=CommonHelpFormatter, add_help=False,
                                       help='Converts a Change-O database to an AIRR database.',
                                       description='Converts a Change-O database to an AIRR database.')
    parser_airr.set_defaults(func=convertDbAIRR)

    # Subparser to convert AIRR to changeo files
    parser_airr = subparsers.add_parser('changeo', parents=[default_parent],
                                       formatter_class=CommonHelpFormatter, add_help=False,
                                       help='Converts an AIRR database to a Change-O database.',
                                       description='Converts an AIRR database to a Change-O database.')
    parser_airr.set_defaults(func=convertDbChangeo)

    # Subparser to convert database entries to sequence file
    parser_fasta = subparsers.add_parser('fasta', parents=[default_parent],
                                       formatter_class=CommonHelpFormatter, add_help=False,
                                       help='Creates a fasta file from database records.',
                                       description='Creates a fasta file from database records.')
    group_fasta = parser_fasta.add_argument_group('conversion arguments')
    group_fasta.add_argument('--if', action='store', dest='id_field',
                              default=default_id_field,
                              help='The name of the field containing identifiers')
    group_fasta.add_argument('--sf', action='store', dest='seq_field',
                              default=default_seq_field,
                              help='The name of the field containing sequences')
    group_fasta.add_argument('--mf', nargs='+', action='store', dest='meta_fields',
                              help='List of annotation fields to add to the sequence description')
    parser_fasta.set_defaults(func=convertDbFasta)
    
    # Subparser to convert database entries to clip-fasta file
    parser_baseln = subparsers.add_parser('baseline', parents=[default_parent],
                                          formatter_class=CommonHelpFormatter, add_help=False,
                                          description='Creates a BASELINe fasta file from database records.',
                                          help='''Creates a specially formatted fasta file
                                               from database records for input into the BASELINe
                                               website. The format groups clonally related sequences
                                               sequentially, with the germline sequence preceding
                                               each clone and denoted by headers starting with ">>".''')
    group_baseln = parser_baseln.add_argument_group('conversion arguments')
    group_baseln.add_argument('--if', action='store', dest='id_field',
                               default=default_id_field,
                               help='The name of the field containing identifiers')
    group_baseln.add_argument('--sf', action='store', dest='seq_field',
                               default=default_seq_field,
                               help='The name of the field containing reads')
    group_baseln.add_argument('--gf', action='store', dest='germ_field',
                               default=default_germ_field,
                               help='The name of the field containing germline sequences')
    group_baseln.add_argument('--cf', action='store', dest='cluster_field', default=None,
                               help='The name of the field containing containing sorted clone IDs')
    group_baseln.add_argument('--mf', nargs='+', action='store', dest='meta_fields',
                               help='List of annotation fields to add to the sequence description')
    parser_baseln.set_defaults(func=convertDbBaseline)

    # Subparser to convert database entries to a GenBank fasta and feature table file
    parser_gb = subparsers.add_parser('genbank', parents=[genbank_parent],
                                       formatter_class=CommonHelpFormatter, add_help=False,
                                       help='Creates a fasta and feature table file for GenBank submissions.',
                                       description='Creates a fasta and feature table file for GenBank submissions.')
    group_gb = parser_gb.add_argument_group('conversion arguments')
    group_gb.add_argument('--product', action='store', dest='product', default=default_product,
                            help='''The product name, such as "immunoglobulin heavy chain".''')
    group_gb.add_argument('--inf', action='store', dest='inference', default=None,
                            help='''Name and version of the inference tool used for reference alignment in the 
                                 form tool:version.''')
    group_gb.add_argument('--db', action='store', dest='db_xref', default=default_db_xref,
                            help='Name of the reference database used for alignment.')
    group_gb.add_argument('--mol', action='store', dest='molecule', default=default_molecule,
                            help='''The source molecule type. Usually one of "mRNA" or "genomic DNA".''')
    group_gb.add_argument('--organism', action='store', dest='organism', default=None,
                            help='The scientific name of the organism.')
    group_gb.add_argument('--sex', action='store', dest='sex', default=None,
                            help='''If specified, adds the given sex annotation 
                                 to the fasta headers.''')
    group_gb.add_argument('--isolate', action='store', dest='isolate', default=None,
                            help='''If specified, adds the given isolate annotation 
                                 (sample label) to the fasta headers.''')
    group_gb.add_argument('--tissue', action='store', dest='tissue', default=None,
                            help='''If specified, adds the given tissue-type annotation 
                                 to the fasta headers.''')
    group_gb.add_argument('--cell', action='store', dest='cell', default=None,
                            help='''If specified, adds the given cell-type annotation 
                                 to the fasta headers.''')
    group_gb.add_argument('--cregion', action='store', dest='cregion_field', default=None,
                            help='''Field containing the C region call. If unspecified, the C region gene 
                                 call will be excluded from the feature table.''')
    group_gb.add_argument('--label', action='store', dest='label', default=None,
                            help='''If specified, add a field name to the sequence identifier. 
                                 Sequence identifiers will be output in the form <label>=<id>.''')
    group_gb.add_argument('--id', action='store_true', dest='keep_id',
                            help='''If specified, use the existing sequence identifier for the output identifier. 
                                 By default, only the row number will be used as the identifier to avoid
                                 the 50 character limit.''')
    group_gb.add_argument('--stop', action='store_true', dest='keep_stop',
                            help='''If specified, retain records in the output with stop codons in the junction region.
                                 In such records the CDS will be removed and replaced with a similar misc_feature in 
                                 the feature table.''')
    parser_gb.set_defaults(func=convertDbGenbank)

    return parser


if __name__ == '__main__':
    """
    Parses command line arguments and calls main function
    """
    # Parse arguments
    parser = getArgParser()
    checkArgs(parser)
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)

    # Check argument pairs
    if args.command == 'add' and len(args_dict['fields']) != len(args_dict['values']):
        parser.error('You must specify exactly one value (-u) per field (-f)')
    elif args.command == 'rename' and len(args_dict['fields']) != len(args_dict['names']):
        parser.error('You must specify exactly one new name (-k) per field (-f)')
    elif args.command == 'update' and len(args_dict['values']) != len(args_dict['updates']):
        parser.error('You must specify exactly one value (-u) per replacement (-t)')

    # Clean arguments dictionary
    del args_dict['command']
    del args_dict['func']
    del args_dict['db_files']
    if 'out_files' in args_dict: del args_dict['out_files']

    # Call main function for each input file
    for i, f in enumerate(args.__dict__['db_files']):
        args_dict['db_file'] = f
        args_dict['out_file'] = args.__dict__['out_files'][i] \
            if args.__dict__['out_files'] else None
        args.func(**args_dict)

