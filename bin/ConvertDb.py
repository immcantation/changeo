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
import shutil
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
from presto.IO import printLog, printMessage, printProgress, printError, printWarning
from changeo.Alignment import gapV
from changeo.Applications import default_tbl2asn_exec, runASN
from changeo.Defaults import default_id_field, default_seq_field, default_germ_field, \
                             default_csv_size, default_format, default_out_args
from changeo.Commandline import CommonHelpFormatter, checkArgs, getCommonArgParser, parseCommonArgs
from changeo.Gene import c_gene_regex, parseAllele, buildGermline
from changeo.IO import countDbFile, getFormatOperators, getOutputHandle, AIRRReader, AIRRWriter, \
                       ChangeoReader, ChangeoWriter, TSVReader, ReceptorData, readGermlines, \
                       checkFields, yamlDict
from changeo.Receptor import AIRRSchema, ChangeoSchema

# System settings
csv.field_size_limit(default_csv_size)

# Defaults
default_db_xref = 'IMGT/GENE-DB'
default_molecule = 'mRNA'
default_product = 'immunoglobulin heavy chain'
default_allele_delim = '*'


def buildSeqRecord(db_record, id_field, seq_field, meta_fields=None):
    """
    Parses a database record into a SeqRecord

    Arguments: 
      db_record : a dictionary containing a database record.
      id_field : the field containing identifiers.
      seq_field : the field containing sequences.
      meta_fields : a list of fields to add to sequence annotations.

    Returns: 
      Bio.SeqRecord.SeqRecord: record.
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


def correctIMGTFields(receptor, references):
    """
    Add IMGT-gaps to IMGT fields in a Receptor object

    Arguments:
      receptor (changeo.Receptor.Receptor): Receptor object to modify.
      references (dict): dictionary of IMGT-gapped references sequences.

    Returns:
      changeo.Receptor.Receptor: modified Receptor with IMGT-gapped fields.
    """
    # Initialize update object
    imgt_dict = {'sequence_imgt': None,
                 'v_germ_start_imgt': None,
                 'v_germ_length_imgt': None,
                 'germline_imgt': None}

    try:
        if not all([receptor.sequence_imgt,
                    receptor.v_germ_start_imgt,
                    receptor.v_germ_length_imgt,
                    receptor.v_call]):
            raise AttributeError
    except AttributeError:
        return None

    # Update IMGT fields
    try:
        gapped = gapV(receptor.sequence_imgt,
                      receptor.v_germ_start_imgt,
                      receptor.v_germ_length_imgt,
                      receptor.v_call,
                      references)
    except KeyError as e:
        printWarning(e)
        return None

    # Verify IMGT-gapped sequence and junction concur
    try:
        check = (receptor.junction == gapped['sequence_imgt'][309:(309 + receptor.junction_length)])
    except TypeError:
        check = False
    if not check:
        return None

    # Rebuild germline sequence
    __, germlines, __ = buildGermline(receptor, references)
    if germlines is None:
        return None
    else:
        gapped['germline_imgt'] = germlines['full']

    # Update return object
    imgt_dict.update(gapped)

    return imgt_dict


def insertGaps(db_file, references=None, format=default_format,
               out_file=None, out_args=default_out_args):
    """
    Inserts IMGT numbering into V fields

    Arguments:
      db_file : the database file name.
      references : folder with germline repertoire files. If None, do not updated alignment columns wtih IMGT gaps.
      format : input format.
      out_file : output file name. Automatically generated from the input file if None.
      out_args : common output argument dictionary from parseCommonArgs.

    Returns:
     str : output file name
    """
    log = OrderedDict()
    log['START'] = 'ConvertDb'
    log['COMMAND'] = 'imgt'
    log['FILE'] = os.path.basename(db_file)
    printLog(log)

    # Define format operators
    try:
        reader, writer, schema = getFormatOperators(format)
    except ValueError:
        printError('Invalid format %s.' % format)

    # Open input
    db_handle = open(db_file, 'rt')
    db_iter = reader(db_handle)

    # Check for required columns
    try:
        required = ['sequence_imgt', 'v_germ_start_imgt']
        checkFields(required, db_iter.fields, schema=schema)
    except LookupError as e:
        printError(e)

    # Load references
    reference_dict = readGermlines(references)

    # Check for IMGT-gaps in germlines
    if all('...' not in x for x in reference_dict.values()):
        printWarning('Germline reference sequences do not appear to contain IMGT-numbering spacers. Results may be incorrect.')

    # Open output writer
    if out_file is not None:
        pass_handle = open(out_file, 'w')
    else:
        pass_handle = getOutputHandle(db_file, out_label='gap', out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'], out_type=schema.out_type)
    pass_writer = writer(pass_handle, fields=db_iter.fields)

    # Count records
    result_count = countDbFile(db_file)

    # Iterate over records
    start_time = time()
    rec_count = pass_count = 0
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time=start_time)
        rec_count += 1
        # Update IMGT fields
        imgt_dict = correctIMGTFields(rec, reference_dict)
        # Write records
        if imgt_dict is not None:
            pass_count += 1
            rec.setDict(imgt_dict, parse=False)
            pass_writer.writeReceptor(rec)

    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time=start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['PASS'] = pass_count
    log['FAIL'] = rec_count - pass_count
    log['END'] = 'ConvertDb'
    printLog(log)

    # Close file handles
    pass_handle.close()
    db_handle.close()

    return pass_handle.name


def convertToAIRR(db_file, format=default_format,
                  out_file=None, out_args=default_out_args):
    """
    Converts a Change-O formatted file into an AIRR formatted file

    Arguments:
      db_file : the database file name.
      format : input format.
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

    # Define format operators
    try:
        reader, __, schema = getFormatOperators(format)
    except ValueError:
        printError('Invalid format %s.' % format)

    # Open input
    db_handle = open(db_file, 'rt')
    db_iter = reader(db_handle)

    # Set output fields replacing length with end fields
    in_fields = [schema.toReceptor(f) for f in db_iter.fields]
    out_fields = [ReceptorData.length_fields[f][1] if f in ReceptorData.length_fields else f \
                  for f in in_fields]
    out_fields = [AIRRSchema.fromReceptor(f) for f in out_fields]

    # Open output writer
    if out_file is not None:
        pass_handle = open(out_file, 'w')
    else:
        pass_handle = getOutputHandle(db_file, out_label='airr', out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'], out_type=AIRRSchema.out_type)
    pass_writer = AIRRWriter(pass_handle, fields=out_fields)

    # Count records
    result_count = countDbFile(db_file)

    # Iterate over records
    start_time = time()
    rec_count = 0
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time=start_time)
        rec_count += 1
        # Write records
        pass_writer.writeReceptor(rec)

    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time=start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['RECORDS'] = rec_count
    log['END'] = 'ConvertDb'
    printLog(log)

    # Close file handles
    pass_handle.close()
    db_handle.close()

    return pass_handle.name


def convertToChangeo(db_file, out_file=None, out_args=default_out_args):
    """
    Converts an AIRR formatted file into an Change-O formatted file

    Arguments:
      db_file: the database file name.
      out_file : output file name. Automatically generated from the input file if None.
      out_args : common output argument dictionary from parseCommonArgs.

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
    db_iter = AIRRReader(db_handle)

    # Set output fields replacing length with end fields
    in_fields = [AIRRSchema.toReceptor(f) for f in db_iter.fields]
    out_fields = [ReceptorData.end_fields[f][1] if f in ReceptorData.end_fields else f \
                  for f in in_fields]
    out_fields = [ChangeoSchema.fromReceptor(f) for f in out_fields]

    # Open output writer
    if out_file is not None:
        pass_handle = open(out_file, 'w')
    else:
        pass_handle = getOutputHandle(db_file, out_label='changeo', out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'], out_type=ChangeoSchema.out_type)
    pass_writer = ChangeoWriter(pass_handle, fields=out_fields)

    # Count records
    result_count = countDbFile(db_file)

    # Iterate over records
    start_time = time()
    rec_count = 0
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time=start_time)
        rec_count += 1
        # Write records
        pass_writer.writeReceptor(rec)

    # Print counts
    printProgress(rec_count, result_count, 0.05, start_time=start_time)
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
def convertToBaseline(db_file, id_field=default_id_field, seq_field=default_seq_field,
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
    rec_count, germ_count, pass_count, fail_count = 0, 0, 0, 0
    cluster_last = None
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time=start_time)
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
    printProgress(rec_count, result_count, 0.05, start_time=start_time)
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


def convertToFasta(db_file, id_field=default_id_field, seq_field=default_seq_field,
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
    rec_count, pass_count, fail_count = 0, 0, 0
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time=start_time)
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
    printProgress(rec_count, result_count, 0.05, start_time=start_time)
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
                        c_field=None, allow_stop=False, asis_calls=False,
                        allele_delim=default_allele_delim):
    """
    Creates a feature table for GenBank submissions

    Arguments:
      record : Receptor record.
      start : start position of the modified seqeuence in the input sequence. Used for feature position offsets.
      end : end position of the modified seqeuence in the input sequence. Used for feature position offsets.
      product : Product (protein) name.
      inference : Reference alignment tool.
      db_xref : Reference database name.
      c_field : column containing the C region gene call.
      allow_stop : if True retain records with junctions having stop codons.
      asis_calls : if True do not parse gene calls for IMGT nomenclature.
      allele_delim : delimiter separating the gene name from the allele number when asis_calls=True.

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

    # Get genes and alleles
    c_gene = None
    if not asis_calls:
        # V gene
        v_gene = record.getVGene()
        v_allele = record.getVAlleleNumber()
        # D gene
        d_gene = record.getDGene()
        d_allele = record.getDAlleleNumber()
        # J gene
        j_gene = record.getJGene()
        j_allele = record.getJAlleleNumber()
        # C region
        if c_field is not None:
            c_gene = parseAllele(record.getField(c_field), c_gene_regex, action='first')
    else:
        # V gene
        v_split = iter(record.v_call.rsplit(allele_delim, maxsplit=1))
        v_gene = next(v_split, None)
        v_allele = next(v_split, None)
        # D gene
        d_split = iter(record.d_call.rsplit(allele_delim, maxsplit=1))
        d_gene = next(d_split, None)
        d_allele = next(d_split, None)
        # J gene
        j_split = iter(record.j_call.rsplit(allele_delim, maxsplit=1))
        j_gene = next(j_split, None)
        j_allele = next(j_split, None)
        # C region
        if c_field is not None:
            c_gene = record.getField(c_field)

    # Fail if V or J is missing
    if v_gene is None or j_gene is None:
        return None

    # Set position offsets if required
    start_trim = 0 if start is None else start
    end_trim = 0 if end is None else len(record.sequence_input) - end
    source_len = len(record.sequence_input) - end_trim

    # Define return object
    result = OrderedDict()

    # C_region
    #   gene
    #   db_xref
    #   inference
    c_region_start = record.j_seq_end + 1 - start_trim
    c_region_length = len(record.sequence_input[(c_region_start + start_trim - 1):]) - end_trim
    if c_region_length > 0:
        if c_gene is not None:
            c_region = [('gene', c_gene),
                        ('db_xref', '%s:%s' % (db_xref, c_gene))]
        else:
            c_region = []

        # Assign C_region feature
        c_region_end = c_region_start + c_region_length - 1
        result[(c_region_start, '>%i' % c_region_end, 'C_region')] = c_region

        # Preserve J segment end position
        j_end = record.j_seq_end

        # Check for range error
        if c_region_end > source_len:
            return None
    else:
        # Trim J segment end position
        j_end = record.j_seq_end + c_region_length

    # V_region
    variable_start = max(record.v_seq_start - start_trim, 1)
    variable_end = j_end - start_trim
    result[(variable_start, variable_end, 'V_region')] = []

    # Check for range error
    if variable_end > source_len:
        return None

    # Product feature
    result[(variable_start, variable_end, 'misc_feature')] = [('note', '%s variable region' % product)]

    # V_segment
    #   gene (gene name)
    #   allele (allele only, without gene name, don't use if ambiguous)
    #   db_xref (database link)
    #   inference (reference alignment tool)
    v_segment = [('gene', v_gene)]
    if v_allele is not None:
        v_segment.append(('allele', v_allele))
    if db_xref is not None:
        v_segment.append(('db_xref', '%s:%s' % (db_xref, v_gene)))
    if inference is not None:
        v_segment.append(('inference', 'COORDINATES:alignment:%s' % inference))
    result[(variable_start, record.v_seq_end - start_trim, 'V_segment')] = v_segment

    # D_segment
    #   gene
    #   allele
    #   db_xref
    #   inference
    if d_gene:
        d_segment = [('gene', d_gene)]
        if d_allele is not None:
            d_segment.append(('allele', d_allele))
        if db_xref is not None:
            d_segment.append(('db_xref', '%s:%s' % (db_xref, d_gene)))
        if inference is not None:
            d_segment.append(('inference', 'COORDINATES:alignment:%s' % inference))
        result[(record.d_seq_start - start_trim, record.d_seq_end - start_trim, 'D_segment')] = d_segment

    # J_segment
    #   gene
    #   allele
    #   db_xref
    #   inference
    j_segment = [('gene', j_gene)]
    if j_allele is not None:
        j_segment.append(('allele', j_allele))
    if db_xref is not None:
        j_segment.append(('db_xref', '%s:%s' % (db_xref, j_gene)))
    if inference is not None:
        j_segment.append(('inference', 'COORDINATES:alignment:%s' % inference))
    result[(record.j_seq_start - start_trim, j_end - start_trim, 'J_segment')] = j_segment

    # CDS
    #   codon_start (must indicate codon offset)
    #   function = JUNCTION
    #   inference
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
        junction_aa = Seq(junction_seq).translate()

        # Return invalid record upon junction stop codon
        if '*' in junction_aa and not allow_stop:
            return None
        elif '*' in junction_aa:
            cds_record.append(('note', '%s junction region' % product))
            result[(cds_start, cds_end, 'misc_feature')] = cds_record
        else:
            cds_record.append(('product', '%s junction region' % product))
            cds_record.append(('codon_start', 1))
            result[(cds_start, cds_end, 'CDS')] = cds_record

    return result


def makeGenbankSequence(record, name=None, label=None, count_field=None, index_field=None,
                        molecule=default_molecule, features=None):
    """
    Creates a sequence for GenBank submissions

    Arguments:
      record : Receptor record.
      name : sequence identifier for the output sequence. If None,
             use the original sequence identifier.
      label : a string to use as a label for the ID. if None do not add a field label.
      count_field : field name to populate the AIRR_READ_COUNT note.
      index_field : field name to populate the AIRR_CELL_INDEX note.
      molecule : source molecule (eg, "mRNA", "genomic DNA")
      features : dictionary of sample features (BioSample attributes) to add to the description of each record.

    Returns:
      dict: dictionary with {'record': SeqRecord,
                             'start': start position in raw sequence,
                             'end': end position in raw sequence}
    """
    # Replace gaps with N
    seq = record.sequence_input
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
    if features is not None:
        sample_desc = ' '.join(['[%s=%s]' % (k, v) for k, v in features.items()])
        name = '%s %s' % (name, sample_desc)
    name = '%s [moltype=%s] [keyword=TLS; Targeted Locus Study; AIRR; MiAIRR:1.0]' % (name, molecule)

    # Notes
    note_dict = OrderedDict()
    if count_field is not None:
        note_dict['AIRR_READ_COUNT'] = record.getField(count_field)
    if index_field is not None:
        note_dict['AIRR_CELL_INDEX'] = record.getField(index_field)
    if note_dict:
        note = '; '.join(['%s:%s' % (k, v) for k, v in note_dict.items()])
        name = '%s [note=%s]' % (name, note)

    # Return SeqRecord and positions
    record = SeqRecord(Seq(seq[seq_start:seq_end], IUPAC.ambiguous_dna), id=name,
                       name=name, description='')
    result = {'record': record, 'start': seq_start, 'end': seq_end}

    return result


def convertToGenbank(db_file, inference=None, db_xref=None, molecule=default_molecule,
                     product=default_product, features=None, c_field=None, label=None,
                     count_field=None, index_field=None, allow_stop=False,
                     asis_id=False, asis_calls=False, allele_delim=default_allele_delim,
                     build_asn=False, asn_template=None, tbl2asn_exec=default_tbl2asn_exec,
                     format=default_format, out_file=None,
                     out_args=default_out_args):
    """
    Builds GenBank submission fasta and table files

    Arguments:
      db_file : the database file name.
      inference : reference alignment tool.
      db_xref : reference database link.
      molecule : source molecule (eg, "mRNA", "genomic DNA")
      product : Product (protein) name.
      features : dictionary of sample features (BioSample attributes) to add to the description of each record.
      c_field : column containing the C region gene call.
      label : a string to use as a label for the ID. if None do not add a field label.
      count_field : field name to populate the AIRR_READ_COUNT note.
      index_field : field name to populate the AIRR_CELL_INDEX note.
      allow_stop : if True retain records with junctions having stop codons.
      asis_id : if True use the original sequence ID for the output IDs.
      asis_calls : if True do not parse gene calls for IMGT nomenclature.
      allele_delim : delimiter separating the gene name from the allele number when asis_calls=True.
      build_asn : if True run tbl2asn on the generated .tbl and .fsa files.
      asn_template : template file (.sbt) to pass to tbl2asn.
      tbl2asn_exec : name of or path to the tbl2asn executable.
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

    # Define format operators
    try:
        reader, __, schema = getFormatOperators(format)
    except ValueError:
        printError('Invalid format %s.' % format)

    # Open input
    db_handle = open(db_file, 'rt')
    db_iter = reader(db_handle)

    # Check for required columns
    try:
        required = ['sequence_input',
                    'v_call', 'd_call', 'j_call',
                    'v_seq_start', 'd_seq_start', 'j_seq_start',
                    'np1_length', 'np2_length']
        checkFields(required, db_iter.fields, schema=schema)
    except LookupError as e:
        printError(e)

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
    rec_count, pass_count, fail_count = 0, 0, 0
    for rec in db_iter:
        # Print progress for previous iteration
        printProgress(rec_count, result_count, 0.05, start_time=start_time)
        rec_count += 1

        # Extract table dictionary
        name = None if asis_id else rec_count
        seq = makeGenbankSequence(rec, name=name, label=label, count_field=count_field, index_field=index_field,
                                  molecule=molecule, features=features)
        tbl = makeGenbankFeatures(rec, start=seq['start'], end=seq['end'], product=product,
                                  db_xref=db_xref, inference=inference, c_field=c_field,
                                  allow_stop=allow_stop, asis_calls=asis_calls, allele_delim=allele_delim)

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

    # Final progress bar
    printProgress(rec_count, result_count, 0.05, start_time=start_time)

    # Run tbl2asn
    if build_asn:
        start_time = time()
        printMessage('Running tbl2asn', start_time=start_time, width=25)
        result = runASN(fsa_handle.name, template=asn_template, exec=tbl2asn_exec)
        printMessage('Done', start_time=start_time, end=True, width=25)

    # Print ending console log
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
                 V_SEQ_START, V_SEQ_LENGTH, D_SEQ_START, D_SEQ_LENGTH, J_SEQ_START, J_SEQ_LENGTH,
                 NP1_LENGTH, NP2_LENGTH
                 SEQUENCE_IMGT, V_GERM_START_IMGT, V_GERM_LENGTH_IMGT
                 
             optional fields:
                 GERMLINE_IMGT, GERMLINE_IMGT_D_MASK, CLONE, C_CALL
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
    format_parent = getCommonArgParser(failed=False, log=False)

    # Subparser to convert changeo to AIRR files
    parser_airr = subparsers.add_parser('airr', parents=[default_parent],
                                        formatter_class=CommonHelpFormatter, add_help=False,
                                        help='Converts input to an AIRR TSV file.',
                                        description='Converts input to an AIRR TSV file.')
    parser_airr.set_defaults(func=convertToAIRR)

    # Subparser to convert AIRR to changeo files
    parser_changeo = subparsers.add_parser('changeo', parents=[default_parent],
                                       formatter_class=CommonHelpFormatter, add_help=False,
                                       help='Converts input into a Change-O TSV file.',
                                       description='Converts input into a Change-O TSV file.')
    parser_changeo.set_defaults(func=convertToChangeo)

    # Subparser to insert IMGT-gaps
    # desc_gap = dedent('''
    #                   Inserts IMGT numbering spacers into the observed sequence
    #                   (SEQUENCE_IMGT, sequence_alignment) and rebuilds the germline sequence
    #                   (GERMLINE_IMGT, germline_alignment) if present. Also adjusts the values
    #                   in the V germline coordinate fields (V_GERM_START_IMGT, V_GERM_LENGTH_IMGT;
    #                   v_germline_end, v_germline_start), which are required.
    #                   ''')
    # parser_gap = subparsers.add_parser('gap', parents=[format_parent],
    #                                     formatter_class=CommonHelpFormatter, add_help=False,
    #                                     help='Inserts IMGT numbering spacers into the V region.',
    #                                     description=desc_gap)
    # group_gap = parser_gap.add_argument_group('conversion arguments')
    # group_gap.add_argument('-r', nargs='+', action='store', dest='references', required=False,
    #                         help='''List of folders and/or fasta files containing
    #                                 IMGT-gapped germline sequences corresponding to the
    #                                 set of germlines used for the alignment.''')
    # parser_gap.set_defaults(func=insertGaps)

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
    parser_fasta.set_defaults(func=convertToFasta)
    
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
    parser_baseln.set_defaults(func=convertToBaseline)

    # Subparser to convert database entries to a GenBank fasta and feature table file
    parser_gb = subparsers.add_parser('genbank', parents=[format_parent],
                                       formatter_class=CommonHelpFormatter, add_help=False,
                                       help='Creates files for GenBank/TLS submissions.',
                                       description='Creates files for GenBank/TLS submissions.')
    # Genbank source information arguments
    group_gb_src = parser_gb.add_argument_group('source information arguments')
    group_gb_src.add_argument('--product', action='store', dest='product', default=default_product,
                          help='''The product name, such as "immunoglobulin heavy chain".''')
    group_gb_src.add_argument('--inf', action='store', dest='inference', default=None,
                          help='''Name and version of the inference tool used for reference alignment in the 
                               form tool:version.''')
    group_gb_src.add_argument('--db', action='store', dest='db_xref', default=default_db_xref,
                          help='Name of the reference database used for alignment.')
    group_gb_src.add_argument('--mol', action='store', dest='molecule', default=default_molecule,
                          help='''The source molecule type. Usually one of "mRNA" or "genomic DNA".''')
    # Genbank sample information arguments
    group_gb_sam = parser_gb.add_argument_group('sample information arguments')
    group_gb_sam.add_argument('--organism', action='store', dest='organism', default=None,
                          help='The scientific name of the organism.')
    group_gb_sam.add_argument('--sex', action='store', dest='sex', default=None,
                          help='''If specified, adds the given sex annotation 
                               to the fasta headers.''')
    group_gb_sam.add_argument('--isolate', action='store', dest='isolate', default=None,
                          help='''If specified, adds the given isolate annotation 
                               (sample label) to the fasta headers.''')
    group_gb_sam.add_argument('--tissue', action='store', dest='tissue', default=None,
                          help='''If specified, adds the given tissue-type annotation 
                               to the fasta headers.''')
    group_gb_sam.add_argument('--cell-type', action='store', dest='cell_type', default=None,
                          help='''If specified, adds the given cell-type annotation 
                               to the fasta headers.''')
    group_gb_sam.add_argument('-y', action='store', dest='yaml_config', default=None,
                          help='''A yaml file specifying sample features (BioSample attributes) 
                               in the form \'variable: value\'. If specified, any features provided in the 
                               yaml file will override those provided at the commandline. Note,
                               this config file applies to sample features only and
                               cannot be used for required source features such as 
                               the --product or --mol argument.''')
    # General genbank conversion arguments
    group_gb_cvt = parser_gb.add_argument_group('conversion arguments')
    group_gb_cvt.add_argument('--label', action='store', dest='label', default=None,
                              help='''If specified, add a field name to the sequence identifier. 
                                Sequence identifiers will be output in the form <label>=<id>.''')
    group_gb_cvt.add_argument('--cf', action='store', dest='c_field', default=None,
                              help='''Field containing the C region call. If unspecified, the C region gene 
                                call will be excluded from the feature table.''')
    group_gb_cvt.add_argument('--nf', action='store', dest='count_field', default=None,
                              help='''If specified, use the provided column to add the AIRR_READ_COUNT 
                                note to the feature table.''')
    group_gb_cvt.add_argument('--if', action='store', dest='index_field', default=None,
                              help='''If specified, use the provided column to add the AIRR_CELL_INDEX 
                                note to the feature table.''')
    group_gb_cvt.add_argument('--allow-stop', action='store_true', dest='allow_stop',
                              help='''If specified, retain records in the output with stop codons in the junction region.
                                In such records the CDS will be removed and replaced with a similar misc_feature in 
                                the feature table.''')
    group_gb_cvt.add_argument('--asis-id', action='store_true', dest='asis_id',
                              help='''If specified, use the existing sequence identifier for the output identifier. 
                                By default, only the row number will be used as the identifier to avoid
                                the 50 character limit.''')
    group_gb_cvt.add_argument('--asis-calls', action='store_true', dest='asis_calls',
                              help='''Specify to prevent alleles from being parsed using the IMGT nomenclature.
                                Note, this requires the gene assignments to be exact matches to valid 
                                records in the references database specified by the --db argument.''')
    group_gb_cvt.add_argument('--allele-delim', action='store', dest='allele_delim', default=default_allele_delim,
                              help='''The delimiter to use for splitting the gene name from the allele number.
                                Note, this only applies when specifying --asis-calls. By default,
                                this argument will be ignored and allele numbers extracted under the
                                expectation of IMGT nomenclature consistency.''')
    group_gb_cvt.add_argument('--asn', action='store_true', dest='build_asn',
                          help='''If specified, run tbl2asn to generate the .sqn submission file after making 
                               the .fsa and .tbl files.''')
    group_gb_cvt.add_argument('--sbt', action='store', dest='asn_template', default=None,
                          help='''If provided along with --asn, use the specified file for the template file
                               argument to tbl2asn.''')
    group_gb_cvt.add_argument('--exec', action='store', dest='tbl2asn_exec', default=default_tbl2asn_exec,
                          help='The name or location of the tbl2asn executable.')
    parser_gb.set_defaults(func=convertToGenbank)

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

    # Genbank mode
    if args.command == 'genbank':
        # Create sample feature dictionary
        features = OrderedDict()
        if args_dict['organism'] is not None:
            features['organism'] = args_dict['organism']
        if args_dict['sex'] is not None:
            features['sex'] = args_dict['sex']
        if args_dict['isolate'] is not None:
            features['isolate'] = args_dict['isolate']
        if args_dict['tissue'] is not None:
            features['tissue_type'] = args_dict['tissue']
        if args_dict['cell_type'] is not None:
            features['cell_type'] = args_dict['cell_type']
        # Add yaml values
        if args_dict['yaml_config'] is not None:
            yaml_config = args_dict['yaml_config']
            if not os.path.exists(yaml_config):
                parser.error('%s does not exist' % yaml_config)
            else:
                features.update(yamlDict(yaml_config))
        # Update arguments
        args_dict['features'] = features
        # Clean arguments dictionary
        del args_dict['organism']
        del args_dict['sex']
        del args_dict['isolate']
        del args_dict['tissue']
        del args_dict['cell_type']
        del args_dict['yaml_config']

        # Check tbl2asn execution arguments
        if args_dict['build_asn']:
            template_file = args_dict['asn_template']
            tbl2asn_exec = args_dict['tbl2asn_exec']
            if not shutil.which(tbl2asn_exec):
                parser.error('%s does not exist or is not executable.' % tbl2asn_exec)
            if template_file is not None and not os.path.exists(template_file):
                parser.error('%s does not exist' % template_file)

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
