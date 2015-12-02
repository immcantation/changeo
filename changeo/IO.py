"""
File I/O and logging functions
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from changeo import __version__, __date__

# Imports
import csv
import os
import sys
from Bio import SeqIO

# Presto and changeo imports
from changeo.Receptor import IgRecord, parseAllele, allele_regex
from presto.IO import getFileType


def getRepo(repo):
    """
    Parses germline repositories

    Arguments:
      repo : String list of directories and/or files
             from which to read germline records

    Returns:
      dictionary : Dictionary of {allele: sequence} germlines
    """
    repo_files = []
    # Iterate over items passed to commandline
    for r in repo:
        # If directory, get fasta files from within
        if os.path.isdir(r):
            repo_files.extend([os.path.join(r, f) for f in os.listdir(r) \
                          if getFileType(f) == 'fasta'])
        # If file, make sure file is fasta
        if os.path.isfile(r) and getFileType(r) == 'fasta':
            repo_files.extend([r])

    # Catch instances where no valid fasta files were passed in
    if len(repo_files) < 1:
        sys.exit('ERROR: No valid germline fasta files were found in %s', repo)

    repo_dict = {}
    for file_name in repo_files:
        with open(file_name, "rU") as file_handle:
            germlines = SeqIO.parse(file_handle, "fasta")
            for g in germlines:
                germ_key = parseAllele(g.description, allele_regex, 'list')
                repo_dict[germ_key] = str(g.seq).upper() # @UndefinedVariable
    return repo_dict


# TODO:  change to require output fields rather than in_file? probably better that way.
def getDbWriter(out_handle, in_file=None, add_fields=None, exclude_fields=None):
    """
    Opens a writer object for an output database file

    Arguments:
      out_handle : File handle to write to
      in_file : the input filename to determine output fields from;
                if None do not define output fields from input file
      add_fields : a list of fields added to the writer not present in the in_file;
                   if None do not add fields
      exclude_fields : a list of fields in the in_file excluded from the writer;
                     if None do not exclude fields

    Returns:
      DictWriter : csv.Dictwriter
    """
    # Get output field names from input file
    if in_file is not None:
        fields = (readDbFile(in_file, ig=False)).fieldnames
    else:
        fields = []
    # Add extra fields
    if add_fields is not None:
        if not isinstance(add_fields, list):  add_fields = [add_fields]
        fields.extend([f for f in add_fields if f not in fields])
    # Remove unwanted fields
    if exclude_fields is not None:
        if not isinstance(exclude_fields, list):  exclude_fields = [exclude_fields]
        fields = [f for f in fields if f not in exclude_fields]

    # Create writer
    try:
        fields = [n.strip().upper() for n in fields]
        # >>> THIS NEEDS TO BE FIXED, extrasaction='ignore' IS A WORKAROUND FOR ADDITIONS TO IgRecord
        db_writer = csv.DictWriter(out_handle, fieldnames=fields, dialect='excel-tab', extrasaction='ignore')
        db_writer.writeheader()
    except:
        sys.exit('ERROR:  File %s cannot be written' % out_handle.name)

    return db_writer


# TODO:  Need to close db_handle?
def readDbFile(db_file, ig=True):
    """
    Reads database files

    Arguments:
      db_file : Tab delimited database file
      ig : If True convert fields to an IgRecord

    Returns:
      iterable : database record iterator
    """
    # Read and check file
    try:
        db_handle = open(db_file, 'rt')
        db_reader = csv.DictReader(db_handle, dialect='excel-tab')
        db_reader.fieldnames = [n.strip().upper() for n in db_reader.fieldnames]
        if ig:
            db_iter = (IgRecord(r) for r in db_reader)
        else:
            db_iter = db_reader
    except IOError:
        sys.exit('ERROR:  File %s cannot be read' % db_file)
    except:
        sys.exit('ERROR:  File %s is invalid' % db_file)

    return db_iter


def countDbFile(db_file):
    """
    Counts the records in database files

    Arguments:
      db_file : Tab delimited database file

    Returns:
      int : Count of records in the database file
    """
    # Count records and check file
    try:
        with open(db_file, 'rt') as db_handle:
            db_records = csv.reader(db_handle, dialect='excel-tab')
            for i, __ in enumerate(db_records):  pass
        db_count = i
    except IOError:
        sys.exit('ERROR:  File %s cannot be read' % db_file)
    except:
        sys.exit('ERROR:  File %s is invalid' % db_file)
    else:
        if db_count == 0:  sys.exit('ERROR:  File %s is empty' % db_file)

    return db_count