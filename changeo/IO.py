"""
File I/O and logging functions
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from changeo import __version__, __date__

# Imports
import csv

# Presto and changeo imports
from changeo.Receptor import IgRecord


# TODO:  change to require output fields rather than in_file? probably better that way.
def getDbWriter(out_handle, in_file=None, add_fields=None, exclude_fields=None):
    """
    Opens a writer object for an output database file

    Arguments:
    out_handle = the file handle to write to
    in_file = the input filename to determine output fields from;
              if None do not define output fields from input file
    add_fields = a list of fields added to the writer not present in the in_file;
                 if None do not add fields
    exclude_fields = a list of fields in the in_file excluded from the writer;
                     if None do not exclude fields

    Returns:
    a writer object
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
    db_file = a tab delimited database file
    ig = if True convert fields to an IgRecord

    Returns:
    a database record iterator
    """
    # Read and check file
    try:
        db_handle = open(db_file, 'rb')
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
    db_file = a tab delimited database file

    Returns:
    the count of records in the database file
    """
    # Count records and check file
    try:
        with open(db_file) as db_handle:
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