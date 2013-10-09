#!/usr/bin/env python
"""
Core functions shared by pRESTO modules
"""

__author__    = 'Jason Anthony Vander Heiden'
__copyright__ = 'Copyright 2013 Kleinstein Lab, Yale University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported'
__version__   = '0.4.0'
__date__      = '2013.9.24'

# Imports
import math, os, re, sys
import multiprocessing as mp
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from itertools import izip, izip_longest, product
from collections import OrderedDict
from time import time, strftime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

# Defaults
default_delimiter = ('|', '=', ',')
default_coord_choices = ['illumina', 'solexa', 'sra', '454', 'presto']
default_action_choices = ['min', 'max', 'sum', 'first', 'last', 'set']
default_coord_type = 'presto'
default_barcode_field = 'BARCODE'
default_primer_field = 'PRIMER'
default_missing_chars = ['-', '.', 'N']
default_min_freq = 0.7
default_min_qual = 20
default_out_args = {'log_file':None, 
                    'delimiter':default_delimiter,
                    'out_dir':None,
                    'out_name':None,
                    'out_type':None,
                    'clean':False}


def parseAnnotation(record, fields=None, delimiter=default_delimiter):
    """
    Extracts annotations from a FASTA/FASTQ sequence description

    Arguments:
    record = description string to extract annotations from
    fields = a list of fields to subset the return dictionary to
             if None return all fields
    delimiter = a tuple of delimiters for (fields, values, value lists) 
    
    Returns: 
    a OrderedDict of {field:value}
    """
    annotation = record.split(delimiter[0])
    field_dict = OrderedDict([('ID', annotation.pop(0))])
    for a in annotation:
        vals = a.split(delimiter[1])
        field_dict[vals[0].upper()] = vals[1] 

    # Subset return dictionary to requested fields
    if fields is not None:
        if not isinstance(fields, list):  fields = [fields]
        for f in set(field_dict).difference(fields):  del field_dict[f]
    
    return field_dict


def flattenAnnotation(ann_dict, delimiter=default_delimiter):
    """
    Convered annotations from a dictionary to a FASTA/FASTQ sequence description

    Arguments: 
    ann_dict = a dictionary of field/value pairs
    delimiter = a tuple of delimiters for (fields, values, value lists) 
    
    Returns: 
    a formatted sequence description string
    """
    annotation = ann_dict.get('ID', 'NONE')
    for k, v in ann_dict.iteritems():
        if k.upper() == 'ID':  continue
        if isinstance(v, list):
            v = delimiter[2].join([str(x) for x in v])
        annotation += '%s%s%s%s' % (delimiter[0], k.upper(), delimiter[1], v) 

    return annotation


def mergeAnnotation(ann_dict_1, ann_dict_2, delimiter=default_delimiter):
    """
    Merges non-ID field annotations from one field dictionary into another

    Arguments:
    ann_dict_1 = a dictionary of field/value pairs to add to
    ann_dict_2 = a dictionary of field/value pairs to merge with fields_2
    delimiter = a tuple of delimiters for (fields, values, value lists) 
    
    Returns: 
    a modified fields_1 dictonary of {field:value}
    """
    for k, v in ann_dict_2.iteritems():
        if k.upper() == 'ID':  continue
        if k in ann_dict_1:  
            if isinstance(ann_dict_1[k], list):
                ann_dict_1[k] = delimiter[2].join([str(x) for x in ann_dict_1[k]])
            if isinstance(v, list):
                v = delimiter[2].join([str(x) for x in v])
            ann_dict_1[k] += '%s%s' % (delimiter[2], v)
        else:  
            ann_dict_1[k.upper()] = v

    return ann_dict_1


def renameAnnotation(ann_dict, old_field, new_field, delimiter=default_delimiter):
    """
    Renames an annotation and merges annotations if the new name already exists

    Arguments:
    ann_dict = a dictionary of field/value pairs
    old_field = the old field name
    new_field = the new field name
    delimiter = a tuple of delimiters for (fields, values, value lists) 
    
    Returns: 
    a modified fields dictonary of {field:value}
    """
    if new_field in ann_dict:
        new_ann = ann_dict.copy()
        del new_ann[old_field]
        mergeAnnotation(new_ann, {new_field:ann_dict[old_field]}, delimiter=delimiter)
    else:
        new_ann = OrderedDict([(new_field, v) if k == old_field else (k, v) \
                               for k, v in ann_dict.iteritems()])

    return new_ann


def collapseAnnotation(ann_dict, action, fields=None, delimiter=default_delimiter):
    """
    Convert annotations from a dictionary to a FASTA/FASTQ sequence description

    Arguments: 
    ann_dict = a dictionary of field/value pairs
    action = the collapse action to take;
             one of {min, max, sum, first, last, set}
    fields = a subset of ann_dict to collapse;
             if None collapse all but the ID field
    delimiter = a tuple of delimiters for (fields, values, value lists) 
    
    Returns: 
    a modified field dictionary of {field:value}
    """
    # Define collapse action            
    if action == 'set':
        def collapse(value):  return list(set(value))
    elif action == 'first':
        def collapse(value):  return value[0]
    elif action == 'last':
        def collapse(value):  return value[-1]
    elif action == 'min':
        def collapse(value):  return '%.2g' % min([float(x or 0) for x in value])
    elif action == 'max':
        def collapse(value):  return '%.2g' % max([float(x or 0) for x in value])
    elif action == 'sum':
        def collapse(value):  return '%.2g' % sum([float(x or 0) for x in value])
    else:
        def collapse(value):  return value

    # Collapse fields        
    for k, v in ann_dict.iteritems():
        if k.upper() == 'ID':
            continue
        if fields is None or k in fields:
            # Convert field to list
            if not isinstance(v, list) and isinstance(v, str):
                v = v.split(delimiter[2])
            elif not isinstance(v, list):
                v = [v]
            # Perform collapse and reassign field
            ann_dict[k] = collapse(v)       

    return ann_dict


def getAnnotationValues(seq_iter, field, unique=False, delimiter=default_delimiter):
    """
    Gets the set of unique annotation values in a sequence set

    Arguments: 
    seq_iter = an iterator or list of SeqRecord objects
    field = the annotation field to retrieve values for
    unique = if True return a list of only the unique values;
             if Falser return a list of all values
    delimiter = a tuple of delimiters for (fields, values, value lists) 

    Returns: 
    the list of values for the field
    """
    # Parse annotations from seq_list records
    ann_iter = (parseAnnotation(s.description, delimiter=delimiter) for s in seq_iter)
    values = [a[field] for a in ann_iter]

    return list(set(values)) if unique else values


def annotationConsensus(seq_iter, field, delimiter=default_delimiter):
    """
    Calculate a consensus annotation for a set of sequences
    
    Arguments:
    seq_iter = an iterator or list of SeqRecord objects
    field = the annotation field to take a consensus of
    delimiter = a tuple of delimiters for (annotations, field/values, value lists) 
    
    Returns:
    a dictionary {set:list of annotation values, count:annotation counts,
                  cons:consensus annotation,freq:majority annotation frequency)
    """
    # Define return dictionary
    cons_dict = {'set':None, 'count':None, 'cons':None, 'freq':None}
    
    # Parse annotations from seq_list records
    val_list = getAnnotationValues(seq_iter, field, delimiter=delimiter)
    
    # Define annotation set and counts
    cons_dict['set'] = sorted(set(val_list))
    cons_dict['count'] = [val_list.count(v) for v in cons_dict['set']]

    # Define consensus annotation
    i = cons_dict['count'].index(max(cons_dict['count']))
    cons_dict['cons'] = cons_dict['set'][i]
    cons_dict['freq'] = float(cons_dict['count'][i]) / len(val_list)
                
    return cons_dict


def subsetSeqSet(seq_iter, field, values, delimiter=default_delimiter):
    """
    Subsets a sequence set by annotation value
    
    Arguments:
    seq_iter = an iterator or list of SeqRecord objects
    field = the annotation field to select by
    values = a list of annotation values that define the retained sequences
    delimiter = a tuple of delimiters for (annotations, field/values, value lists) 
    
    Returns:
    a modified list of SeqRecord objects
    """
    # Parse annotations from seq_list records
    ann_list = [parseAnnotation(s.description, delimiter=delimiter) for s in seq_iter]
    
    # Subset seq_list by annotation
    if not isinstance(values, list):  values = [values]
    seq_subset = [seq_iter[i] for i, a in enumerate(ann_list) if a[field] in values]

    return seq_subset


def subsetSeqIndex(seq_dict, field, values, delimiter=default_delimiter):
    """
    Subsets a sequence set by annotation value
    
    Arguments:
    seq_dict = a dictionary index of sequences returned from SeqIO.index()
    field = the annotation field to select keys by
    values = a list of annotation values that define the retained keys
    delimiter = a tuple of delimiters for (annotations, field/values, value lists) 
    
    Returns:
    a list of keys
    """
    # Parse annotations from seq_dict and subset keys
    key_subset = [k for k in seq_dict \
                  if parseAnnotation(seq_dict[k].description, delimiter=delimiter)[field] \
                  in values]

    return key_subset


def compilePrimers(primers):
    """
    Translates IUPAC Ambiguous Nucleotide characters to regular expressions and compiles them

    Arguments: 
    key = a dictionary of sequences to translate

    Returns:
    dictionary of compiled regular expressions
    """
    
    primers_regex = {k: re.compile(re.sub(r'([RYSWKMBDHVN])', translateIUPAC, v)) 
                     for k, v in primers.iteritems()}  
    
    return primers_regex


def readPrimerFile(primer_file):
    """
    Processes primer sequences from file

    Arguments: 
    primer_file = name of file containing primer sequences
        
    Returns: 
    dictionary of {primer id, primer sequence}
    """
    
    # Parse primers from .fasta and .regex files
    ext_name = os.path.splitext(primer_file)[1]
    if ext_name == '.fasta':
        with open(primer_file, 'rU') as primer_handle:
            primer_iter = SeqIO.parse(primer_handle, 'fasta', IUPAC.ambiguous_dna)        
            primers = {p.description: str(p.seq).upper() 
                        for p in primer_iter}
    elif ext_name == '.regex':
        with open(primer_file, 'rU') as primer_handle:
            primer_list = [a.split(':') for a in primer_handle]
            primers = {a[0].strip(): re.sub(r'\[([^\[^\]]+)\]', translateIUPAC, a[1].strip().upper()) 
                        for a in primer_list}
    else:
        sys.exit('ERROR:  %s primer file is not a supported type' % ext_name.upper())
        primers = None
    
    return primers


def readSeqFile(seq_file, index=False):
    """
    Reads FASTA/FASTQ files

    Arguments: 
    seq_file = a FASTA or FASTQ file containing sample sequences
    index = if True return a dictionary from SeqIO.index();
            if False return an iterator from SeqIO.parse()
    Returns: 
    a tuple of (input file type, sequence record object)
    """
    # Read and check file
    try:
        seq_type = os.path.splitext(seq_file)[1].lower().lstrip('.')
        if index:  
            seq_records = SeqIO.index(seq_file, seq_type, IUPAC.ambiguous_dna)
        else:  
            seq_records = SeqIO.parse(seq_file, seq_type, IUPAC.ambiguous_dna)
    except IOError:
        sys.exit('ERROR:  File %s cannot be read' % seq_file)
    except:
        sys.exit('ERROR:  File %s is invalid' % seq_file)
    
    return seq_type, seq_records


def countSeqFile(seq_file):
    """
    Counts the records in FASTA/FASTQ files

    Arguments: 
    seq_file = a FASTA or FASTQ file containing sample sequences

    Returns: 
    the count of records in the sequence file
    """
    # Count records and check file
    try:
        seq_type = os.path.splitext(seq_file)[1].lower().lstrip('.')
        result_count = len(SeqIO.index(seq_file, seq_type))
    except IOError:
        sys.exit('ERROR:  File %s cannot be read' % seq_file)
    except:
        sys.exit('ERROR:  File %s is invalid' % seq_file)
    else:
        if result_count == 0:  sys.exit('ERROR:  File %s is empty' % seq_file)
        
    return result_count

 
def getOutputHandle(in_file, out_label=None, out_dir=None, out_name=None, out_type=None):
    """
    Opens an output file handle
    
    Arguments: 
    in_file = the input filename
    out_label = text to be inserted before the file extension;
                if None do not add a label
    out_type = the file extension of the output file; 
               if None use input file extension
    out_dir = the output directory;
              if None use directory of input file
    out_name = the short file anme to use for the output file; 
               if None use input file short name
    
    Returns:
    a file handle
    """
    # Get in_file components
    dir_name, file_name = os.path.split(in_file)
    short_name, ext_name = os.path.splitext(file_name)
    
    # Define output directory
    if out_dir is None:
        out_dir = dir_name
    else:
        out_dir = os.path.abspath(out_dir)
        if not os.path.exists(out_dir):  os.mkdir(out_dir)
    # Define output file prefix
    if out_name is None:  out_name = short_name
    # Define output file extension
    if out_type is None:  out_type = ext_name.lstrip('.')

    # Define output file name
    if out_label is None:  
        out_file = os.path.join(out_dir, '%s.%s' % (out_name, out_type))
    else:
        out_file = os.path.join(out_dir, '%s_%s.%s' % (out_name, out_label, out_type))
    
    return open(out_file, 'w')


def translateIUPAC(key):
    """
    Translates IUPAC Ambiguous Nucleotide characters to or from character sets

    Arguments: 
    key = a string or re.search object containing the character set to translate

    Returns:
    character translation
    """
    
    # Define valid characters and character translations
    IUPAC_uniq = '-.ACGT'
    IUPAC_ambig = 'BDHKMNRSVWY'
    IUPAC_trans = {'AG':'R', 'CT':'Y', 'CG':'S', 'AT':'W', 'GT':'K', 'AC':'M',
                  'CGT':'B', 'AGT':'D', 'ACT':'H', 'ACG':'V', 'ABCDGHKMRSTVWY':'N',
                  '-.':'.'}

    # Convert passed regular expression match to a string
    if hasattr(key, 'group'): key = key.group(1)
    
    # Sort character in alphabetic order
    key = ''.join(sorted(key))
    
    # Return input character if no translation needed
    if len(key) == 1 and key in IUPAC_uniq:
        return key
    # Return regular expression string for ambiguous single character
    elif len(key) == 1 and key in IUPAC_ambig:
        return ['[' + k + ']' for k, v in IUPAC_trans.iteritems() if v == key][0]
    # Return single ambiguous character for character set 
    elif key in IUPAC_trans:
        return IUPAC_trans[key]
    else:
        return 'N'


def scoreIUPAC(a, b, n_score=None, gap_score=None):
    """
    Returns the score for a pair of IUPAC Ambiguous Nucleotide characters

    Arguments: 
    a = first character
    b = second character
    n_score = score for all matches against an N character;
              if None score according to IUPAC character identity
    gap_score = score for all matches against a [-, .] character;
                if None score according to IUPAC character identity

    Returns:
    score for the character pair
    """
    # Define ambiguous character translations    
    IUPAC_trans = {'AG':'R', 'CT':'Y', 'CG':'S', 'AT':'W', 'GT':'K', 'AC':'M',
                   'CGT':'B', 'AGT':'D', 'ACT':'H', 'ACG':'V', 'ABCDGHKMRSTVWY':'N',
                   '-.':'.'}
    # Create list of tuples of synonymous character pairs
    IUPAC_matches = [p for k, v in IUPAC_trans.iteritems() for p in list(product(k, v))]

    # Check gap condition
    if gap_score is not None and (a in '-.' or b in '-.'):
        return gap_score

    # Check N-value condition
    if n_score is not None and (a == 'N' or b == 'N'):
        return n_score

    # Determine and return score for IUPAC match conditions
    # Symmetric and reflexive
    if a == b:
        return 1
    elif (a, b) in IUPAC_matches:
        return 1
    elif (b, a) in IUPAC_matches:
        return 1
    else:
        return 0


def getScoreDict(n_score=None, gap_score=None):
    """
    Generates a score dictionary

    Arguments: 
    n_score = score for all matches against an N character;
              if None score according to IUPAC character identity
    gap_score = score for all matches against a [-, .] character;
                if None score according to IUPAC character identity
                
    Returns:
    a score dictionary of the form {(char1, char2) : score}
    """
    IUPAC_chars = '-.ACGTRYSWKMBDHVN'
    IUPAC_dict = {k:scoreIUPAC(*k, n_score=n_score, gap_score=gap_score) 
                  for k in product(IUPAC_chars, repeat=2)}

    return IUPAC_dict


def reverseComplement(seq):
    """
    Takes the reverse complement of a sequence

    Arguments: 
    seq = a SeqRecord object, Seq object or string to reverse complement
    
    Returns:
    a object of the same type as the input with the reverse complement sequence 
    """
    
    if isinstance(seq, SeqRecord):
        new_record = seq.reverse_complement(id=True, name=True, description=True, 
                                            features=True, annotations=True, 
                                            letter_annotations=True)
        #new_record.annotations['orientation'] = 'RC'
    elif isinstance(seq, Seq):
        new_record = seq.reverse_complement()
    elif isinstance(seq, str):
        new_record = str(Seq(seq, IUPAC.ambiguous_dna).reverse_complement())
    else:
        sys.exit('ERROR:  invalid record type passed to reverseComplement')
        new_record = None
    
    return new_record


def testSeqEqual(seq1, seq2, ignore_chars=default_missing_chars):
    """
    Determine if two sequences are equal, excluding missing positions
    
    Arguments: 
    seq1 = a SeqRecord object
    seq2 = a SeqRecord object
    ignore_chars = a list of characters to ignore
    
    Returns:
    True if the sequences are equal
    """
    equal = True
    for a, b in izip(seq1, seq2):
        if a != b and a not in ignore_chars and b not in ignore_chars:
            equal = False
            break
    
    return equal
 

def weightSeq(seq):
    """
    Returns a score for a single sequence
    
    Arguments: 
    seq = a SeqRecord object
    
    Returns:
    The sum of the character scores for the sequence
    """
    nuc_score = sum([c in 'ACGTRYSWKMBDHV' for c in seq])
    gap_score = 0
      
    return max(nuc_score + gap_score, 1)


def scoreSeqPair(seq1, seq2, max_error=None, max_weight=None, 
                 score_dict=getScoreDict(n_score=0, gap_score=0)):
    """
    Determine the error rate for a pair of sequences
    
    Arguments: 
    seq1 = a SeqRecord object
    seq2 = a SeqRecord object
    max_error = the maximum error rate; once reached return (0, 0, 1.0)
                if None always return accurate score, weight and error
    max_weight = the maximum weight to use when checking the max_error break condition;
                 if None use the minimum length of seq1,seq2
    score_dict = optional dictionary of alignment scores as {(char1, char2): score}
    
    Returns:
    A tuple of the (score, minimum weight, error rate) for the pair of sequences
    """
    # Determine score
    if max_error is None:
        # Return accurate values when max_error is undefined
        chars = izip(seq1, seq2)
        score = sum([score_dict[(a, b)] for a, b in chars])
        weight = min(weightSeq(seq1), weightSeq(seq2))
        error = 1.0 - float(score) / weight
    else:
        # If max_error defined return when worst case reach
        score = 0
        if not max_weight:  max_weight = min(len(seq1), len(seq2))
        for i, (a, b) in enumerate(izip(seq1, seq2)):
            score += score_dict[(a, b)]
            if (i - float(score)) / max_weight > max_error:
                score, weight, error = 0, 0, 1.0
                break
        else:
            weight = min(weightSeq(seq1), weightSeq(seq2))
            error = 1.0 - float(score) / weight
            
    
    return (score, weight, error)


def calculateDiversity(seq_list, score_dict=getScoreDict(n_score=0, gap_score=0)):
    """
    Determine the average pairwise error rate for a list of sequences

    Arguments: 
    seq_list = a list of SeqRecord objects to score
    score_dict = optional dictionary of alignment scores as {(char1, char2): score}
        
    Returns:
    The average pairwise error rate for the list of sequences
    """
    # Return 0 if less than 2 sequences
    if len(seq_list) <= 1:
        return 0
    
    scores = []
    for i, seq1 in enumerate(seq_list):
        for seq2 in seq_list[(i + 1):]:
            scores.append(scoreSeqPair(seq1, seq2, score_dict=score_dict)[2])
    
    return sum(scores) / len(scores)


def qualityConsensus(seq_list, min_qual=default_min_qual, dependent=False):
    """
    Builds a consensus sequence from a set of sequences

    Arguments: 
    set_seq = a list of SeqRecord objects
    min_qual = the quality cutoff to assign a base
    dependent = if False assume sequences are independent for quality calculation
    
    Returns: 
    a consensus SeqRecord object
    """
    # Return a copy of the input SeqRecord upon singleton
    if len(seq_list) == 1:
        return seq_list[0].upper()
    
    # Define characters to be ignored in consensus building
    ignore_set = set(['.', '-', 'N'])

    # Create sequence and annotation iterators
    # Pad unequal length sequences with character '-' and quality 0 
    seq_str = [str(s.seq) for s in seq_list]
    seq_iter = izip_longest(*seq_str, fillvalue='-')
    ann_list = [s.letter_annotations['phred_quality'] for s in seq_list]
    ann_iter = izip_longest(*ann_list, fillvalue=0)
    
    # Build consensus
    consensus_seq = []
    consensus_qual = []
    for chars, quals in izip(seq_iter, ann_iter):
        char_set = set(chars).difference(ignore_set)
        # Creater per character quality sets, sums and normalization terms
        qual_total = sum(quals)
        qual_set, qual_sum = {}, {}
        for c in char_set:
            qual_set[c] = [q for i, q in enumerate(quals) if chars[i] == c]
            qual_sum[c] = sum(qual_set[c])
        # Calculate per character consensus quality scores
        if dependent:
            qual_cons = {c:int(float(max(qual_set[c])) * qual_sum[c] / qual_total) for c in qual_set}
        else:
            qual_cons = {c:int(float(qual_sum[c]) * qual_sum[c] / qual_total) for c in qual_set}
        # Select character with highest quality over cutoff
        if qual_cons.values():
            qual_cons_max = max(qual_cons.values())
            cons = [(c, min(q, 93)) if q >= min_qual else ('N', min(q, 93)) \
                    for c, q in qual_cons.iteritems() if q == qual_cons_max][0]
        else:
            cons = ('N', 0)
        consensus_seq.append(cons[0])
        consensus_qual.append(cons[1])
    
    # Define return SeqRecord
    record = SeqRecord(Seq(''.join(consensus_seq), IUPAC.ambiguous_dna),
                       id='consensus',
                       name='consensus',
                       description='',
                       letter_annotations={'phred_quality':consensus_qual})

    return record


def frequencyConsensus(seq_list, min_freq=default_min_freq):
    """
    Builds a consensus sequence from a set of sequences

    Arguments: 
    set_seq = a list of SeqRecord objects
    min_freq = the frequency cutoff to assign a base
    
    Returns: 
    a consensus SeqRecord object
    """
    # Return a copy of the input SeqRecord upon singleton
    if len(seq_list) == 1:
        return seq_list[0].upper()

    # Build consensus
    seq_str = [str(s.seq) for s in seq_list]
    consensus_seq = []
    for chars in izip_longest(*seq_str, fillvalue='-'):
        count = float(len(chars))
        freq_dict = {c:chars.count(c) / count for c in set(chars)}
        freq_max = max(freq_dict.values())
        cons = [c if freq_dict[c] >= min_freq else 'N' \
                for c in chars if freq_dict[c] == freq_max][0]
        consensus_seq.append(cons)

    # Define return SeqRecord
    record = SeqRecord(Seq(''.join(consensus_seq), IUPAC.ambiguous_dna), 
                       id='consensus', 
                       name='consensus', 
                       description='')
        
    return record


def indexSeqSets(seq_dict, field=default_barcode_field, delimiter=default_delimiter):
    """
    Identifies sets of sequences with the same ID field

    Arguments: 
    seq_dict = a dictionary index of sequences returned from SeqIO.index()
    field = the annotation field containing set IDs
    delimiter = a tuple of delimiters for (fields, values, value lists) 
    
    Returns: 
    a dictionary of {set name:[record names]}
    """
    set_dict = {}
    for key, rec in seq_dict.iteritems():
        tag = parseAnnotation(rec.description, delimiter=delimiter)[field]
        set_dict.setdefault(tag, []).append(key)

    return set_dict


def indexSeqPairs(seq_dict_1, seq_dict_2, coord_type=default_coord_type, 
                  delimiter=default_delimiter):
    """
    Identifies sequence pairs by coordination information

    Arguments: 
    seq_dict_1 = a dictionary index of sequences returned from SeqIO.index()
    seq_dict_2 = a dictionary index of sequences returned from SeqIO.index()
    coord_type = the sequence header format
    delimiter = a tuple of delimiters for (fields, values, value lists) 
    
    Returns: 
    a dictionary of {coordinate:(seq_dict_1 key, seq_dict_2 key)}
    """
    # Get coordinate dictionaries
    coord_1 = getSeqCoord(seq_dict_1, coord_type, delimiter=delimiter)
    coord_2 = getSeqCoord(seq_dict_2, coord_type, delimiter=delimiter)
    
    # Form sets of keys
    set_1 = set(coord_1)
    set_2 = set(coord_2)
    
    # Find matching entries in key sets
    match_set = set_1.intersection(set_2)
    index_dict = {n:(coord_1[n], coord_2[n]) for n in match_set}
    
    return index_dict


def getUnpairedIndex(seq_dict_1, seq_dict_2, coord_type=default_coord_type, 
                     delimiter=default_delimiter):
    """
    Identifies sequence without a mate pair by coordination information

    Arguments: 
    seq_dict_1 = a dictionary index of sequences returned from SeqIO.index()
    seq_dict_2 = a dictionary index of sequences returned from SeqIO.index()
    coord_type = the sequence header format
    delimiter = a tuple of delimiters for (fields, values, value lists) 
    
    Returns: 
    a tuple of lists containing unique ([seq_dict_1 keys], [seq_dict_2 keys])
    """
    # Get coordinate dictionaries
    coord_1 = getSeqCoord(seq_dict_1, coord_type, delimiter=delimiter)
    coord_2 = getSeqCoord(seq_dict_2, coord_type, delimiter=delimiter)
    
    # Form sets of coordinate keys
    keys_1 = set(coord_1)
    keys_2 = set(coord_2)
    
    # Find unmatched entries in coordinate key sets
    uniq_keys_1 = keys_1.difference(keys_2)
    uniq_keys_2 = keys_2.difference(keys_1)
    
    # Translate unmatched coordinate keys to sequences dictionary keys
    unpaired_1 = [coord_1[n] for n in uniq_keys_1]
    unpaired_2 = [coord_2[n] for n in uniq_keys_2]
    
    return unpaired_1, unpaired_2


def getSeqCoord(seq_dict, coord_type=default_coord_type, delimiter=default_delimiter):
    """
    Create a sequence ID translation using delimiter truncation
    
    Arguments: 
    seq_dict = a sequence dictionary generate by SeqIO.index()
    coord_type = the sequence header format; 
                 one of ['illumina', 'solexa', 'sra', '454', 'presto'];
                 if unrecognized type or None create ID dictionary without translation
    delimiter = a tuple of delimiters for (fields, values, value lists) 
                    
    Returns: 
    a dictionary of {coordinate: full ID}
    """
    # Sequencing identifier formats supported
    # illumina:  @MISEQ:132:000000000-A2F3U:1:1101:14340:1555 2:N:0:ATCACG
    #            @HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1
    # sra:       @SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
    # 454:       @000034_0199_0169 length=437 uaccno=GNDG01201ARRCR
    # pesto      @AATCGGATTTGC|COUNT=2|PRIMER=IGHJ_RT|PRFREQ=1.0
     
    # Define getCoord function by coordinate format
    if coord_type in ('illumina', 'solexa'):
        def getCoord(s):  return s.split()[0].split('#')[0]
    elif coord_type in ('sra', '454'):
        def getCoord(s):  return s.split()[0]
    elif coord_type == 'presto':
        def getCoord(s):  return parseAnnotation(s, delimiter=delimiter)['ID']
    else:
        def getCoord(s):  return s

    # Create a seq_dict ID translation using IDs truncated by getCoord()
    ids = {}
    for seq_id in seq_dict:
        try:
            id_key = getCoord(seq_id)
            ids.update({id_key:seq_id})
        except:
            sys.exit('ERROR:  No ID for sequence %s using coordinate format %s' \
                     % (seq_id, coord_type))
    
    return ids


def feedRecQueue(data_queue, nproc, seq_iter):
    """
    Feeds the data queue with individual sequences for processQueue processes

    Arguments: 
    data_queue = an multiprocessing.Queue to hold data for processing
    nproc = the number of processQueue processes
    seq_iter = an iterator of SeqRecords returned by readSeqFile
    
    Returns: 
    None
    """
    # Iterate over seq_ter and define processQueue arguments 
    for seq in seq_iter:
        # Feed queue
        data_queue.put({'id':seq.id, 'seq':seq})
    
    # Add sentinel object for each processQueue process
    for __ in range(nproc):
        data_queue.put(None)

    return None


def feedSetQueue(data_queue, nproc, seq_dict, index_dict):
    """
    Feeds the data queue with sequence sets for processQueue processes

    Arguments: 
    data_queue = a multiprocessing.Queue to hold data for processing
    nproc = the number of processQueue processes
    seq_dict = a dictionary of SeqRecords returned by readSeqFile
    set_dict = a dictionary returned by indexSeqSets defining set membership
    
    Returns: 
    None
    """
    # Iterate over set_dict and define processQueue arguments 
    for set_id, set_members in index_dict.iteritems():
        # Define sequence set
        seq_list = [seq_dict[k] for k in set_members]
        # Feed queue
        data_queue.put({'id':set_id, 'seq_list':seq_list})
        
    # Add sentinel object for each processQueue process
    for __ in range(nproc):
        data_queue.put(None)

    return None


def feedPairQueue(data_queue, nproc, seq_dict_1, seq_dict_2, index_dict):
    """
    Feeds the data queue with sequence pairs for processQueue processes

    Arguments: 
    data_queue = an multiprocessing.Queue to hold data for processing
    nproc = the number of processQueue processes
    seq_dict_1 = a dictionary of SeqRecords returned by readSeqFile
    seq_dict_2 = a dictionary of SeqRecords returned by readSeqFile
    index_dict = a dictionary returned by indexSeqPairs
         
    Returns: 
    None
    """
    # Iterate over id_dict added head and tail sequences pairs to data_queue
    for pair_id, (key_1, key_2) in index_dict.iteritems():
        # Feed queue
        data_queue.put({'id':pair_id, 
                        'seq_1':seq_dict_1[key_1], 
                        'seq_2':seq_dict_2[key_2]})
        
    # Add sentinel object for each processQueue process
    for __ in range(nproc):
        data_queue.put(None)

    return None


def collectRecQueue(result_queue, result_count, collect_dict, task_label, seq_file, out_args):
    """
    Assembles results from a queue of individual sequence results and manages log/file I/O

    Arguments: 
    result_queue = a multiprocessing.Queue holding processQueue results
    result_count = the total number of results expected
    collect_dict = a multiprocessing.Manager.dict to store return values
    task_label = the task label used to tag the output files
    seq_file = the input sequence file name
    out_args = common output argument dictionary from parseCommonArgs
    
    Returns: 
    None
    (adds 'log' and 'out_files' to collect_dict)
    """
    # Defined valid alignment output handle
    pass_handle = getOutputHandle(seq_file, 
                                  '%s-pass' % task_label, 
                                  out_dir=out_args['out_dir'], 
                                  out_name=out_args['out_name'], 
                                  out_type=out_args['out_type'])
    # Defined failed alignment output handle
    if out_args['clean']:   
        fail_handle = None
    else:  
        fail_handle = getOutputHandle(seq_file, 
                                      '%s-fail' % task_label, 
                                      out_dir=out_args['out_dir'], 
                                      out_name=out_args['out_name'], 
                                      out_type=out_args['out_type'])
    # Define log handle
    if out_args['log_file'] is None:  
        log_handle = None
    else:  
        log_handle = open(out_args['log_file'], 'w')
        
    # Iterator over results queue until sentinel object reached
    start_time = time()
    seq_count = pass_count = fail_count = 0
    for result in iter(result_queue.get, None): 
        # Print progress for previous iteration
        printProgress(seq_count, result_count, 0.05, start_time) 
        
        # Update counts
        seq_count += 1

        # Write log
        printLog(result['log'], handle=log_handle)

        # Write sequences
        if result['out_seq'] is not None and result['pass']:
            pass_count += 1
            SeqIO.write(result['out_seq'], pass_handle, out_args['out_type'])
        else:
            fail_count += 1
            if fail_handle is not None:
                SeqIO.write(result['in_seq'], fail_handle, out_args['out_type'])
 
    # Print total counts
    printProgress(seq_count, result_count, 0.05, start_time) 

    # Update return list
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['SEQUENCES'] = seq_count
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    collect_dict['log'] = log
    collect_dict['out_files'] = [pass_handle.name]

    # Close file handles
    pass_handle.close()
    if fail_handle is not None:  fail_handle.close()
    if log_handle is not None:  log_handle.close()
        
    return None


def collectSetQueue(result_queue, result_count, collect_dict, task_label, seq_file, out_args):
    """
    Pulls from results queue, assembles results and manages log and file IO

    Arguments: 
    result_queue = a multiprocessing.Queue holding processQueue results
    result_count = the total number of results expected
    collect_dict = a multiprocessing.Manager.dict to store return values
    task_label = the task label used to tag the output files
    seq_file = the sample sequence file name
    out_args = common output argument dictionary from parseCommonArgs
    
    Returns: 
    None
    (adds 'log' and 'out_files' to collect_dict)
    """
    # Defined valid alignment output handle
    pass_handle = getOutputHandle(seq_file, 
                                  '%s-pass' % task_label, 
                                  out_dir=out_args['out_dir'], 
                                  out_name=out_args['out_name'], 
                                  out_type=out_args['out_type'])
    # Defined failed alignment output handle
    if out_args['clean']:  
        fail_handle = None
    else:  
        fail_handle = getOutputHandle(seq_file, 
                                      '%s-fail'  % task_label, 
                                      out_dir=out_args['out_dir'], 
                                      out_name=out_args['out_name'], 
                                      out_type=out_args['out_type'])
    # Define log handle
    if out_args['log_file'] is None:  
        log_handle = None
    else:  
        log_handle = open(out_args['log_file'], 'w')
    
    # Iterator over results queue until sentinel object reached
    start_time = time()
    set_count = seq_count = single_count = pass_count = fail_count = 0
    for result in iter(result_queue.get, None): 
        # Print progress for previous iteration
        printProgress(set_count, result_count, 0.05, start_time) 
        
        # Update counts for current iteration
        set_count += 1
        seq_count += result['seq_count']
        if result['seq_count'] == 1:  single_count += 1
        
        # Write log
        printLog(result['log'], handle=log_handle)

        # Write alignments
        if result['out_list'] is not None and result['pass']:
            pass_count += 1
            SeqIO.write(result['out_list'], pass_handle, out_args['out_type'])
        else:
            fail_count += 1
            if fail_handle is not None:
                SeqIO.write(result['in_list'], fail_handle, out_args['out_type'])
 
    # Print total counts
    printProgress(set_count, result_count, 0.05, start_time)

    # Update return values    
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['SEQUENCES'] = seq_count
    log['SINGLETONS'] = single_count
    log['SETS'] = set_count
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    collect_dict['log'] = log
    collect_dict['out_files'] = [pass_handle.name]

    # Close file handles
    pass_handle.close()
    if fail_handle is not None:  fail_handle.close()
    if log_handle is not None:  log_handle.close()
    
    return None


def collectPairQueue(result_queue, result_count, collect_dict, task_label, 
                     seq_file_1, seq_file_2, out_args):
    """
    Pulls from results queue, assembles results and manages log and file IO

    Arguments: 
    result_queue = a multiprocessing.Queue holding processQueue results
    result_count = the total number of results expected
    collect_dict = a multiprocessing.Manager.dict to store return values
    task_label = the task label used to tag the output files
    seq_file_1 = the first sequence file name
    seq_file_2 = the second sequence file name
    out_args = common output argument dictionary from parseCommonArgs
    
    Returns: 
    None
    (adds 'log' and 'out_files' to collect_dict)
    """
    # Defined valid assembly output handle
    pass_handle = getOutputHandle(seq_file_1, 
                                  '%s-pass' % task_label, 
                                  out_dir=out_args['out_dir'], 
                                  out_name=out_args['out_name'], 
                                  out_type=out_args['out_type'])
    # Defined failed assembly output handles
    if out_args['clean']:
        fail_handle_1 = None
        fail_handle_2 = None
    else:
        # Define output name
        if out_args['out_name'] is None:
            out_name_1 = out_name_2 = None
        else: 
            out_name_1 = '%s-1' % out_args['out_name']
            out_name_2 = '%s-2' % out_args['out_name']
        fail_handle_1 = getOutputHandle(seq_file_1, 
                                        '%s-fail' % task_label, 
                                        out_dir=out_args['out_dir'], 
                                        out_name=out_name_1, 
                                        out_type=out_args['out_type'])
        fail_handle_2 = getOutputHandle(seq_file_2, 
                                        '%s-fail' % task_label, 
                                        out_dir=out_args['out_dir'], 
                                        out_name=out_name_2, 
                                        out_type=out_args['out_type'])

    # Define log handle
    if out_args['log_file'] is None:
        log_handle = None
    else:
        log_handle = open(out_args['log_file'], 'w')
        
    # Iterator over results queue until sentinel object reached
    start_time = time()
    pair_count = pass_count = fail_count = 0
    for result in iter(result_queue.get, None): 
        # Print progress for previous iteration
        printProgress(pair_count, result_count, 0.05, start_time)

        # Update counts for iteration
        pair_count += 1

        # Write log
        printLog(result['log'], handle=log_handle)

        # Write assembled sequences
        if result['out_seq'] is not None and result['pass']:
            pass_count += 1
            SeqIO.write(result['out_seq'], pass_handle, out_args['out_type'])
        else:
            fail_count += 1
            if fail_handle_1 is not None and fail_handle_2 is not None:
                SeqIO.write(result['in_seq_1'], fail_handle_1, out_args['out_type'])
                SeqIO.write(result['in_seq_2'], fail_handle_2, out_args['out_type'])
 
    # Print total counts
    printProgress(pair_count, result_count, 0.05, start_time)

    # Update return values
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['PAIRS'] = pair_count
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    collect_dict['log'] = log
    collect_dict['out_files'] = [pass_handle.name]

    # Close file handles
    pass_handle.close()
    if fail_handle_1 is not None:  fail_handle_1.close()
    if fail_handle_2 is not None:  fail_handle_2.close()
    if log_handle is not None:  log_handle.close()

    
    return None


def printLog(record, handle=sys.stdout, inset=None):
    """
    Formats a dictionary into an IgPipeline log string

    Arguments: 
    record = a dict or OrderedDict of {field names: values}
    handle = the file handle to write the log to;
             if None do not write to file
    inset = minimum field name inset;
            if None automatically space field names
                    
    Returns:
    a formatted multi-line string in IgPipeline log format
    """
    # Return empty string if empty dictionary passed
    if not record:
        return ''
    
    # Determine inset
    max_len = max(map(len, record))
    inset = max(max_len, inset)
    
    # Assemble log string
    record_str = ''
    if isinstance(record, OrderedDict):
        key_list = record.keys()
    else:
        key_list = sorted(record)
    for key in key_list:
        record_str += '%s%s> %s\n' % (' ' * (inset - len(key)), key, record[key])
    
    # Write log record
    if handle is not None:
        try:
            handle.write('%s\n' % record_str)
        except IOError as e:
            sys.stderr.write('I/O error writing to log file: %s\n' % e.strerror)

    return record_str


def printProgress(current, total=None, step=None, start_time=None, end=False):
    """
    Prints a progress bar to standard out
    
    Arguments:
    current = the count of completed tasks 
    total = the total task count;
            if None do not print percentage
    step = a float defining the fractional progress increment to print if total is defined;
           an int defining the progress increment to print at if total is not defined;
           if None always output the progress
    start_time = task start time returned by time.time();
                 if None do not add run time to progress
    end = if True print final log (add newline)
    
    Returns:
    None
    """
    try:
        # Check update condition
        if total is None:
            update = (current % step == 0)
        else:
            update = (current % math.ceil(step*total) == 0)
    except:
        # Return on modulo by zero error
        return None
    else:
        # Check end condition and return if no update needed
        if current == total:
            end = True
        if not end and not update:
            return None
        
    # Define progress bar
    if total is not None and total != 0:
        p = float(current) / total
        c = format(current, "%i,d" % len(format(total, ",d")))
        bar = 'PROGRESS> %s [%-20s] %3.0f%% (%s)' \
              % (strftime('%H:%M:%S'), '#' * int(p*20), p*100, c)
    else:
        bar = 'PROGRESS> %s (%s)' % (strftime('%H:%M:%S'), current)
        
    # Add run time to bar if start_time is specified
    if start_time is not None:
        bar = '%s %.1f min' % (bar, (time() - start_time)/60)
    
    # Print progress bar    
    if current == 0:
        print '%s' % bar,
        sys.stdout.flush()
    elif end:
        print '\r%s\n' % bar
        sys.stdout.flush()
    else:
        print '\r%s' % bar,
        sys.stdout.flush()
    
    return None
        
        
def getCommonParser(seqin=True, seqout=True, paired=False, annotation=True, log=True, multiproc=False):
    """
    Defines an ArgumentParser object with common pRESTO arguments

    Arguments: 
    seqin = if True include sequence input arguments
    seqout = if True include sequence output arguments
    paired = if True defined paired-end sequence input and output arguments
    annotation = if True include annotation arguments
    log = if True include log arguments
    multiproc = if True include multiprocessing arguments
    
    Returns:
    an ArgumentParser object
    """
    parser = ArgumentParser(add_help=False, formatter_class=ArgumentDefaultsHelpFormatter)

    # Sequence processing arguments
    if seqin and not paired:
        parser.add_argument('-s', nargs='+', action='store', dest='seq_files', required=True,
                            help='List of FASTA/FASTQ files containing sequences')
    elif seqin and paired:
        parser.add_argument('-1', nargs='+', action='store', dest='seq_files_1', required=True,
                            help='Ordered list of FASTA/FASTQ files containing head/primary sequences')
        parser.add_argument('-2', nargs='+', action='store', dest='seq_files_2', required=True,
                            help='Ordered list of FASTA/FASTQ files containing tail/secondary sequences')
    if seqout:
        parser.add_argument('--fasta', action='store_const', dest='out_type', const='fasta',
                            help='Specify to force output as FASTA rather than FASTQ')
        parser.add_argument('--clean', action='store_true', dest='clean', 
                            help='If specified do not create files of sequences that fail processing')
            
    # Annotation arguments
    if annotation:
        parser.add_argument('--delim', nargs=3, action='store', dest='delimiter', 
                            type=str, default=default_delimiter, 
                            help='The three delimiters that separate annotation blocks, \
                                  field names and values, and values within a list, respectively')

    # Log arguments
    if log:
        parser.add_argument('--log', action='store', dest='log_file', default=None,
                            help='Specify to write verbose logging to a file')    
    
    # Multiprocessing arguments
    if multiproc:
        parser.add_argument('--nproc', action='store', dest='nproc', type=int, default=mp.cpu_count(),
                            help='The number of computational simultaneous computational processes to execute')
    
    # Universal arguments
    parser.add_argument('--outdir', action='store', dest='out_dir', default=None,
                        help='Changes the output directory to the path specified \
                              Uses input file directory if not specified')
    parser.add_argument('--outname', action='store', dest='out_name', default=None,
                        help='Changes the prefix of the successfully processed output file to the string specified')
    
    
    
    return parser


def parseCommonArgs(args, file_args=None):
    """
    Checks common arguments from getCommonParser and transforms output options to a dictionary

    Arguments: 
    args = argument Namespace defined by ArgumentParser.parse_args
    file_args = a list of non-common file arguments to verify; 
                ['seq_files', '1', '2', 'primer_file', 'out_dir'] 
                are checked by default
                    
    Returns:
    a dictionary copy of args with output arguments embedded in the dictionary out_args
    """ 
    file_types = ['.fasta', '.fastq']
    primer_types = ['.fasta', '.regex']
    args_dict = args.__dict__.copy()
    
    # Verify sequence files
    seq_files = []
    if 'seq_files' in args_dict:
        seq_files.extend(args_dict['seq_files'])
    elif 'seq_files_1' and 'seq_files_2' in args_dict:
        if len(args_dict['seq_files_1']) != len(args_dict['seq_files_2']):
            sys.exit('ERROR:  The -1 and -2 arguments must contain the same number of files')
        for f1, f2 in zip(args_dict['seq_files_1'], args_dict['seq_files_2']):
            if os.path.splitext(f1)[-1].lower() != os.path.splitext(f2)[-1].lower():
                sys.exit('ERROR:  Each pair of files in the -1 and -2 arguments must be the same file type')
        seq_files.extend(args_dict['seq_files_1'])
        seq_files.extend(args_dict['seq_files_2'])
    for f in seq_files:
        if not os.path.isfile(f):
            sys.exit('ERROR:  sequence file %s does not exist' % f)
        if os.path.splitext(f)[-1].lower() not in file_types:
            sys.exit('ERROR:  sequence file %s is not a supported type. Must be one: %s' \
                     % (','.join(file_types), f))

    # Verify primer files
    if 'primer_file' in args_dict:
        primer_files = args_dict['primer_file'] if isinstance(args_dict['primer_file'], list) \
                       else [args_dict['primer_file']]
        for f in primer_files:
            if not os.path.isfile(f):
                sys.exit('ERROR:  primer file %s does not exist' % f)
            if os.path.splitext(f)[-1].lower() not in file_types:
                sys.exit('ERROR:  primer file %s is not a supported type. Must be one: %s' \
                         % (','.join(primer_types), f))
            
    # Verify output directory
    if 'out_dir' in args_dict and args_dict['out_dir'] is not None:
        if os.path.exists(args_dict['out_dir']) and not os.path.isdir(args_dict['out_dir']):
            sys.exit('ERROR:  directory %s exists and is not a directory' % args_dict['out_dir'])
    
    # Verify optional files
    if file_args is not None:
        if not isinstance(file_args, list):  file_args = [file_args]
        for arg in file_args:
            files = args_dict[arg] if isinstance(args_dict[arg], list) \
                    else [args_dict[arg]]
            for f in files:
                if not os.path.isfile(f):
                    sys.exit('ERROR:  file %s does not exist' % f)
    
    # Redfine common output options as out_args dictionary
    out_args = ['log_file', 'delimiter', 'out_dir', 'out_name', 'out_type', 'clean']
    args_dict['out_args'] = {k:args_dict.setdefault(k, None) for k in out_args}
    for k in out_args: del args_dict[k]
    
    return args_dict


if __name__ == '__main__':
    """
    Print module information
    """
    print 'Version: %s %s %s' % (os.path.basename(__file__), __version__, __date__)
    print 'Location: %s' % os.path.dirname(os.path.realpath(__file__))
    #print 'Parent Dir: %s' % path.join(path.dirname(path.realpath(__file__)), path.pardir)
    