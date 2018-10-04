"""
File I/O and parsers
"""
# Info
__author__ = 'Namita Gupta, Jason Anthony Vander Heiden'

# Imports
import csv
import os
import re
import tarfile
import yaml
import zipfile
from itertools import chain, groupby, zip_longest
from tempfile import TemporaryDirectory
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

# Presto and changeo imports
from presto.IO import getFileType, printError, printWarning
from changeo.Defaults import default_csv_size
from changeo.Gene import allele_regex, v_allele_regex, d_allele_regex, j_allele_regex, \
                         parseAllele
from changeo.Receptor import AIRRSchema, ChangeoSchema, Receptor, ReceptorData
from changeo.Alignment import decodeBTOP, encodeCIGAR, padAlignment, gapV, inferJunction, getRegions

# System settings
csv.field_size_limit(default_csv_size)


class TSVReader:
    """
    Simple csv.DictReader wrapper to read format agnostic TSV files.

    Attributes:
      reader (iter): reader object.
      fields (list): field names.
    """

    def __init__(self, handle):
        """
        Initializer

        Arguments:
          handle : handle to an open TSV file

        Returns:
          changeo.IO.TSVReader
        """
        # Arguments
        self.handle = handle
        self.receptor = False
        self.reader = csv.DictReader(self.handle, dialect='excel-tab')
        self.fields = self.reader.fieldnames

    def __iter__(self):
        """
        Iterator initializer

        Returns:
          changeo.IO.TSVReader
        """
        return self

    def __next__(self):
        """
        Next method

        Returns:
          dist : row as a dictionary of field:value pairs.
        """
        # Get next row from reader iterator
        try:
            record = next(self.reader)
        except StopIteration:
            self.handle.close()
            raise StopIteration

        return self._parse(record)

    def _parse(self, record):
        """
        Parses a dictionary of fields

        Arguments:
          record : dict with fields and values to parse.

        Returns:
          dict : parsed dict.
        """
        return record


class TSVWriter:
    """
    Simple csv.DictWriter wrapper to write format agnostic TSV files.
    """
    def __init__(self, handle, fields, header=True):
        """
        Initializer

        Arguments:
          handle : handle to an open output file
          fields : list of output field names
          header : if True write the header on initialization.

        Returns:
          changeo.IO.TSVWriter
        """
        # Arguments
        self.handle = handle
        self.fields = fields
        self.writer = csv.DictWriter(self.handle, fieldnames=self.fields,
                                     dialect='excel-tab', extrasaction='ignore',
                                     lineterminator='\n')
        if header:
            self.writeHeader()

    def writeHeader(self):
        """
        Writes the header

        Returns:
          None
        """
        self.writer.writeheader()

    def writeDict(self, records):
        """
        Writes a row from a dictionary

        Arguments:
          records : dictionary of row data or an iterable of such objects.

        Returns:
          None
        """
        if isinstance(records, dict):
            self.writer.writerow(records)
        else:
            self.writer.writerows(records)


class ChangeoReader(TSVReader):
    """
    An iterator to read and parse Change-O formatted data.
    """
    def __init__(self, handle):
        """
        Initializer

        Arguments:
          handle : handle to an open Change-O formatted file

        Returns:
          changeo.IO.ChangeoReader
        """
        # Arguments
        self.handle = handle
        self.reader = csv.DictReader(self.handle, dialect='excel-tab')
        self.reader.fieldnames = [n.strip().upper() for n in self.reader.fieldnames]
        self.fields = self.reader.fieldnames

    def _parse(self, record):
        """
        Parses a dictionary to a Receptor object

        Arguments:
          record : dict with fields and values in the Change-O format

        Returns:
          changeo.Receptor.Receptor : parsed Receptor object.
        """
        # Parse fields
        result = {}
        for k, v in record.items():
            k = ChangeoSchema.toReceptor(k)
            result[k] = v

        return Receptor(result)


class ChangeoWriter(TSVWriter):
    """
    Writes Change-O formatted data.
    """
    def __init__(self, handle, fields=ChangeoSchema.standard_fields, header=True):
        """
        Initializer

        Arguments:
          handle : handle to an open output file
          fields : list of output field names
          header : if True write the header on initialization.

        Returns:
          changeo.IO.ChangeoWriter
        """
        # Arguments
        self.handle = handle
        self.fields = [n.strip().upper() for n in fields]
        self.writer = csv.DictWriter(self.handle, fieldnames=self.fields,
                                     dialect='excel-tab', extrasaction='ignore',
                                     lineterminator='\n')
        if header:
            self.writeHeader()

    def _parseReceptor(self, record):
        """
        Parses a Receptor object to a Change-O dictionary

        Arguments:
          record : dict with fields and values in the Receptor format.

        Returns:
          dict : parsed dict.
        """
        row = record.toDict()
        # Parse known fields
        result = {}
        for k, v in row.items():
            k = ChangeoSchema.fromReceptor(k)
            result[k] = v

        return result

    def writeReceptor(self, records):
        """
        Writes a row from a Receptor object

        Arguments:
          records : a changeo.Receptor.Receptor object to write or an iterable of such objects.

        Returns:
          None
        """
        if isinstance(records, Receptor):
            row = self._parseReceptor(records)
            self.writer.writerow(row)
        else:
            rows = (self._parseReceptor(r) for r in records)
            self.writer.writerows(rows)


class AIRRReader(TSVReader):
    """
    An iterator to read and parse AIRR formatted data.
    """
    def __init__(self, handle):
        """
        Initializer

        Arguments:
          handle : handle to an open AIRR formatted file
          receptor : if True (default) iteration returns a Receptor object, otherwise it returns a dictionary.

        Returns:
          changeo.IO.AIRRReader
        """
        # Arguments
        self.handle = handle

        # Define reader
        try:
            import airr
            self.reader = airr.io.RearrangementReader(self.handle, debug=False)
        except ImportError as e:
            printError('AIRR library cannot be imported: %s.' % e)

        # Set field list
        self.fields = self.reader.fields

    def _parse(self, record):
        """
        Parses a dictionary of AIRR records to a Receptor object

        Arguments:
          record : dict with fields and values in the AIRR format.

        Returns:
          changeo.Receptor.Receptor : parsed Receptor object.
        """
        # Parse fields
        result = {}
        for k, v in record.items():
            # Rename fields
            k = AIRRSchema.toReceptor(k)
            # Convert start positions to 0-based
            # if k in ReceptorData.start_fields and v is not None and v != '':
            #     v = str(int(v) + 1)
            # Assign new field
            result[k] = v

        # Assign length based on start and end
        for end, (start, length) in ReceptorData.end_fields.items():
            if end in result and result[end] is not None:
                result[length] = int(result[end]) - int(result[start]) + 1

        return Receptor(result)


class AIRRWriter(TSVWriter):
    """
    Writes AIRR formatted data.
    """
    def __init__(self, handle, fields=AIRRSchema.standard_fields):
        """
        Initializer

        Arguments:
          handle : handle to an open output file
          fields : list of output field names

        Returns:
          changeo.IO.AIRRWriter
        """
        # Arguments
        self.handle = handle
        self.fields = [f.lower() for f in fields]

        # Define writer
        try:
            import airr
            self.writer = airr.io.RearrangementWriter(self.handle, fields=self.fields, debug=False)
        except ImportError as e:
            printError('AIRR library cannot be imported: %s.' % e)

    def _parseReceptor(self, record):
        """
        Parses a Receptor object to an AIRR dictionary

        Arguments:
          record : dict with fields and values in the Receptor format

        Returns:
          dict : a parsed dict.
        """
        result = {}
        row = record.toDict()
        for k, v in row.items():
            # Convert start positions to 0-based
            # if k in ReceptorData.start_fields and v is not None and v != '':
            #     v = str(int(v) - 1)
            # Convert field names
            k = AIRRSchema.fromReceptor(k)
            result[k] = v

        return result

    def writeReceptor(self, records):
        """
        Writes a row from a Receptor object

        Arguments:
          records : a changeo.Receptor object to write or iterable of such objects.

        Returns:
          None
        """
        if isinstance(records, Receptor):
            row = self._parseReceptor(records)
            self.writer.write(row)
        else:
            rows = (self._parseReceptor(r) for r in records)
            for r in rows:  self.writer.write(r)


class IMGTReader:
    """
    An iterator to read and parse IMGT output files.
    """
    @staticmethod
    def customFields(scores=False, regions=False, junction=False, schema=None):
        """
        Returns non-standard fields defined by the parser

        Arguments:
          scores : if True include alignment scoring fields.
          regions : if True include IMGT-gapped CDR and FWR region fields.
          junction : if True include detailed junction annotation fields.
          schema : schema class to pass field through for conversion.
                   If None, return changeo.Receptor.Receptor attribute names.

        Returns:
          list : list of field names.
        """
        # Alignment scoring fields
        score_fields = ['v_score',
                        'v_identity',
                        'j_score',
                        'j_identity']

        # FWR amd CDR fields
        region_fields = ['fwr1_imgt',
                         'fwr2_imgt',
                         'fwr3_imgt',
                         'fwr4_imgt',
                         'cdr1_imgt',
                         'cdr2_imgt',
                         'cdr3_imgt']

        # Define default detailed junction field ordering
        junction_fields = ['n1_length',
                           'n2_length',
                           'p3v_length',
                           'p5d_length',
                           'p3d_length',
                           'p5j_length',
                           'd_frame']

        fields = []
        if scores:  fields.extend(score_fields)
        if regions:  fields.extend(region_fields)
        if junction:  fields.extend(junction_fields)

        # Convert field names if schema provided
        if schema is not None:
            fields = [schema.fromReceptor(f) for f in fields]

        return fields

    def __init__(self, summary, gapped, ntseq, junction, receptor=True):
        """
        Initializer

        Arguments:
          summary : handle to an open '1_Summary' IMGT/HighV-QUEST output file.
          gapped : handle to an open '2_IMGT-gapped-nt-sequences' IMGT/HighV-QUEST output file.
          ntseq: handle to an open '3_Nt-sequences' IMGT/HighV-QUEST output file.
          junction : handle to an open '6_Junction' IMGT/HighV-QUEST output file.
          receptor : if True (default) iteration returns an Receptor object, otherwise it returns a dictionary.

        Returns:
          change.Parsers.IMGTReader
        """
        # Arguments
        self.summary = summary
        self.gapped = gapped
        self.ntseq = ntseq
        self.junction = junction
        self.receptor = receptor

        # Open readers
        readers = [csv.DictReader(self.summary, delimiter='\t'),
                   csv.DictReader(self.gapped, delimiter='\t'),
                   csv.DictReader(self.ntseq, delimiter='\t'),
                   csv.DictReader(self.junction, delimiter='\t')]
        self.records = zip(*readers)

    @staticmethod
    def _parseFunctionality(summary):
        """
        Parse functionality information

        Arguments:
          summary : dictionary containing one row of the '1_Summary' file.

        Returns:
          dict : database entries for functionality information.
        """
        # Correct for new functionality column names
        if 'Functionality' not in summary:
            summary['Functionality'] = summary['V-DOMAIN Functionality']
            summary['Functionality comment'] = summary['V-DOMAIN Functionality comment']

        # Orientation parser
        def _revcomp():
            x = {'+': 'F', '-': 'T'}
            return x.get(summary['Orientation'], None)

        # Functionality parser
        def _functional():
            x = summary['Functionality']
            if x.startswith('productive'):
                return 'T'
            elif x.startswith('unproductive'):
                return 'F'
            else:
                return None

        # Junction frame parser
        def _inframe():
            x = {'in-frame': 'T', 'out-of-frame': 'F'}
            return x.get(summary['JUNCTION frame'], None)

        # Stop codon parser
        def _stop():
            return 'T' if 'stop codon' in summary['Functionality comment'] else 'F'

        # Mutated invariant parser
        def _invariant():
            x = summary['Functionality comment']
            y = summary['V-REGION potential ins/del']
            return 'T' if ('missing' in x) or ('missing' in y) else 'F'

        # Mutated invariant parser
        def _indels():
            x = summary['V-REGION potential ins/del']
            y = summary['V-REGION insertions']
            z = summary['V-REGION deletions']
            return 'T' if any([x, y, z]) else 'F'

        result = {}
        # Parse functionality information
        if 'No results' not in summary['Functionality']:
            result['rev_comp'] = _revcomp()
            result['functional'] = _functional()
            result['in_frame'] = _inframe()
            result['stop'] = _stop()
            result['mutated_invariant'] = _invariant()
            result['indels'] = _indels()

        return result

    @staticmethod
    def _parseGenes(summary):
        """
        Parse gene calls

        Arguments:
          summary : dictionary containing one row of the '1_Summary' file.

        Returns:
          dict : database entries for gene calls.
        """
        clean_regex = re.compile('(,)|(\(see comment\))')
        delim_regex = re.compile('\sor\s')

        # Gene calls
        result = {}
        v_call = summary['V-GENE and allele']
        d_call = summary['D-GENE and allele']
        j_call = summary['J-GENE and allele']
        result['v_call'] = delim_regex.sub(',', clean_regex.sub('', v_call)) if v_call else None
        result['d_call'] = delim_regex.sub(',', clean_regex.sub('', d_call)) if d_call else None
        result['j_call'] = delim_regex.sub(',', clean_regex.sub('', j_call)) if j_call else None

        return result

    @staticmethod
    def _parseSequences(gapped, ntseq):
        """
        Parses full length V(D)J sequences

        Arguments:
          gapped : dictionary containing one row of the '2_IMGT-gapped-nt-sequences' file.
          ntseq: dictionary containing one row of the '3_Nt-sequences' file.

        Returns:
          dict : database entries for fill length V(D)J sequences.
        """
        result = {}

        # Extract ungapped sequences
        if ntseq['V-D-J-REGION']:
            result['sequence_vdj'] = ntseq['V-D-J-REGION']
        elif ntseq['V-J-REGION']:
            result['sequence_vdj'] = ntseq['V-J-REGION']
        else:
            result['sequence_vdj'] = ntseq['V-REGION']

        # Extract gapped sequences
        if gapped['V-D-J-REGION']:
            result['sequence_imgt'] = gapped['V-D-J-REGION']
        elif gapped['V-J-REGION']:
            result['sequence_imgt'] = gapped['V-J-REGION']
        else:
            result['sequence_imgt'] = gapped['V-REGION']

        return result

    @staticmethod
    def _parseVPos(gapped, ntseq):
        """
        Parses V alignment positions

        Arguments:
          gapped : dictionary containing one row of the '2_IMGT-gapped-nt-sequences' file.
          ntseq: dictionary containing one row of the '3_Nt-sequences' file.

        Returns:
          dict : database entries for V query and germline alignment positions.
        """
        result = {}
        result['v_seq_start'] = ntseq['V-REGION start']
        result['v_seq_length'] = len(ntseq['V-REGION']) if ntseq['V-REGION'] else 0
        result['v_germ_start_imgt'] = 1
        result['v_germ_length_imgt'] = len(gapped['V-REGION']) if gapped['V-REGION'] else 0

        return result

    @staticmethod
    def _parseJuncPos(junction, db):
        """
        Parses junction N/P and D alignment positions

        Arguments:
          junction : dictionary containing one row of the '6_Junction' file.
          db : database containing V alignment information.

        Returns:
          dict : database entries for junction, N/P and D region alignment positions.
        """
        v_start = db['v_seq_start']
        v_length = db['v_seq_length']

        # First N/P length
        def _np1():
            nb = [junction['P3\'V-nt nb'],
                  junction['N-REGION-nt nb'],
                  junction['N1-REGION-nt nb'],
                  junction['P5\'D-nt nb']]
            return sum(int(i) for i in nb if i)

        # D start
        def _dstart():
            nb = [v_start,
                  v_length,
                  junction['P3\'V-nt nb'],
                  junction['N-REGION-nt nb'],
                  junction['N1-REGION-nt nb'],
                  junction['P5\'D-nt nb']]
            return sum(int(i) for i in nb if i)

        # Second N/P length
        def _np2():
            nb = [junction['P3\'D-nt nb'],
                  junction['N2-REGION-nt nb'],
                  junction['P5\'J-nt nb']]
            return sum(int(i) for i in nb if i)

        result = {}
        # Junction sequence
        result['junction'] = junction['JUNCTION']
        result['junction_aa'] = junction['JUNCTION (AA)']
        result['junction_length'] = len(junction['JUNCTION']) if junction['JUNCTION'] else 0

        # N/P and D alignment positions
        result['np1_length'] = _np1()
        result['d_seq_start'] = _dstart()
        result['d_seq_length'] = int(junction['D-REGION-nt nb'] or 0)
        result['d_germ_start'] = int(junction['5\'D-REGION trimmed-nt nb'] or 0) + 1
        result['d_germ_length'] = int(junction['D-REGION-nt nb'] or 0)
        result['np2_length'] = _np2()

        return result

    @staticmethod
    def _parseJPos(gapped, ntseq, junction, db):
        """
        Parses J alignment positions

        Arguments:
          gapped : dictionary containing one row of the '2_IMGT-gapped-nt-sequences' file.
          ntseq: dictionary containing one row of the '3_Nt-sequences' file.
          junction : dictionary containing one row of the '6_Junction' file.
          db : database containing V, N/P and D alignment information.

        Returns:
          dict : database entries for J region alignment positions.
        """

        # J start
        def _jstart():
            nb = [db['v_seq_start'],
                  db['v_seq_length'],
                  db['np1_length'],
                  db['d_seq_length'],
                  db['np2_length']]
            return sum(int(i) for i in nb if i)

        # J region alignment positions
        result = {}
        result['j_seq_start'] = _jstart()
        result['j_seq_length'] = len(ntseq['J-REGION']) if ntseq['J-REGION'] else 0
        result['j_germ_start'] = int(junction['5\'J-REGION trimmed-nt nb'] or 0) + 1
        result['j_germ_length'] = len(gapped['J-REGION']) if gapped['J-REGION'] else 0

        return result

    @staticmethod
    def _parseScores(summary):
        """
        Parse alignment scores

        Arguments:
          summary : dictionary containing one row of the '1_Summary' file.

        Returns:
          dict : database entries for alignment scores.
        """
        result = {}

        # V score
        try:
            result['v_score'] = float(summary['V-REGION score'])
        except (TypeError, ValueError):
            result['v_score'] = None
        # V identity
        try:
            result['v_identity'] = float(summary['V-REGION identity %']) / 100.0
        except (TypeError, ValueError):
            result['v_identity'] = 'None'
        # J score
        try:
            result['j_score'] = float(summary['J-REGION score'])
        except (TypeError, ValueError):
            result['j_score'] = None
        # J identity
        try:
            result['j_identity'] = float(summary['J-REGION identity %']) / 100.0
        except (TypeError, ValueError):
            result['j_identity'] = None

        return result

    @staticmethod
    def _parseJuncDetails(junction):
        """
        Parse detailed junction region information

        Arguments:
          junction : dictionary containing one row of the '6_Junction' file.

        Returns:
          dict : database entries for detailed D, N and P region information.
        """

        # D reading frame
        def _dframe():
            frame = None
            x = junction['D-REGION reading frame']
            if x:
                try:
                    frame = int(x)
                except ValueError:
                    m = re.search(r'reading frame ([0-9])', x).group(1)
                    frame = int(m)
            return frame

        # First N region length
        def _n1():
            nb = [junction['N-REGION-nt nb'], junction['N1-REGION-nt nb']]
            return sum(int(i) for i in nb if i)

        # D Frame and junction fields
        result = {}
        result['d_frame'] = _dframe()
        result['n1_length'] = _n1()
        result['n2_length'] = int(junction['N2-REGION-nt nb'] or 0)
        result['p3v_length'] = int(junction['P3\'V-nt nb'] or 0)
        result['p5d_length'] = int(junction['P5\'D-nt nb'] or 0)
        result['p3d_length'] = int(junction['P3\'D-nt nb'] or 0)
        result['p5j_length'] = int(junction['P5\'J-nt nb'] or 0)

        return result

    def parseRecord(self, summary, gapped, ntseq, junction):
        """
        Parses a single row from each IMTG file.

        Arguments:
          summary : dictionary containing one row of the '1_Summary' file.
          gapped : dictionary containing one row of the '2_IMGT-gapped-nt-sequences' file.
          ntseq : dictionary containing one row of the '3_Nt-sequences' file.
          junction : dictionary containing one row of the '6_Junction' file.

        Returns:
          dict: database entry for the row.
        """
        # Check that rows are syncronized
        id_set = [summary['Sequence ID'],
                  gapped['Sequence ID'],
                  ntseq['Sequence ID'],
                  junction['Sequence ID']]
        if len(set(id_set)) != 1:
            printError('IMGT files are corrupt starting with Summary file record %s.' % id_set[0])

        # Initialize db with query ID and sequence
        db = {'sequence_id': summary['Sequence ID'],
              'sequence_input': summary['Sequence']}

        # Parse required fields
        db.update(IMGTReader._parseFunctionality(summary))
        db.update(IMGTReader._parseGenes(summary))
        db.update(IMGTReader._parseSequences(gapped, ntseq))
        db.update(IMGTReader._parseVPos(gapped, ntseq))
        db.update(IMGTReader._parseJuncPos(junction, db))
        db.update(IMGTReader._parseJPos(gapped, ntseq, junction, db))

        # Parse optional fields
        db.update(IMGTReader._parseScores(summary))
        db.update(getRegions(db['sequence_imgt'], db['junction_length']))
        db.update(IMGTReader._parseJuncDetails(junction))

        return db

    def __iter__(self):
        """
        Iterator initializer.

        Returns:
          changeo.IO.IMGTReader
        """
        return self

    def __next__(self):
        """
        Next method.

        Returns:
          changeo.Receptor.Receptor : parsed IMGT/HighV-QUEST result as an Receptor (receptor=True) or dictionary (receptor=False).
        """
        # Get next set of records from dictionary readers
        try:
            summary, gapped, ntseq, junction = next(self.records)
        except StopIteration:
            raise StopIteration

        db = self.parseRecord(summary, gapped, ntseq, junction)

        if self.receptor:
            return Receptor(db)
        else:
            return db


class IgBLASTReader:
    """
    An iterator to read and parse IgBLAST output files
    """
    # Ordered list of known fields
    @staticmethod
    def customFields(scores=False, regions=False, cdr3=False, schema=None):
        """
        Returns non-standard fields defined by the parser

        Arguments:
          scores : if True include alignment scoring fields.
          regions : if True include IMGT-gapped CDR and FWR region fields.
          cdr3 : if True include IgBLAST CDR3 assignment fields.
          schema : schema class to pass field through for conversion.
                   If None, return changeo.Receptor.Receptor attribute names.

        Returns:
          list : list of field names.
        """
        # IgBLAST scoring fields
        score_fields = ['v_score',
                        'v_identity',
                        'v_evalue',
                        'v_cigar',
                        'd_score',
                        'd_identity',
                        'd_evalue',
                        'd_cigar',
                        'j_score',
                        'j_identity',
                        'j_evalue',
                        'j_cigar']

        # FWR amd CDR fields
        region_fields = ['fwr1_imgt',
                         'fwr2_imgt',
                         'fwr3_imgt',
                         'fwr4_imgt',
                         'cdr1_imgt',
                         'cdr2_imgt',
                         'cdr3_imgt']

        # IgBLAST CDR3 fields
        cdr3_fields = ['cdr3_igblast',
                       'cdr3_igblast_aa']

        fields = []
        if scores:  fields.extend(score_fields)
        if regions:  fields.extend(region_fields)
        if cdr3:  fields.extend(cdr3_fields)

        # Convert field names if schema provided
        if schema is not None:
            fields = [schema.fromReceptor(f) for f in fields]

        return fields

    def __init__(self, igblast, sequences, references, asis_calls=False, receptor=True):
        """
        Initializer.

        Arguments:
          igblast (file): handle to an open IgBLAST output file written with '-outfmt 7 std qseq sseq btop'.
          sequences (dict): dictionary of query sequences;
                            sequence descriptions as keys with original query sequences as SeqRecord values.
          references (dict): dictionary of IMGT gapped germline sequences.
          asis_calls (bool): if True do not parse gene calls for allele names.
          receptor (bool): if True (default) iteration returns an Receptor object, otherwise it returns a dictionary.

        Returns:
          changeo.IO.IgBLASTReader
        """
        # Arguments
        self.igblast = igblast
        self.sequences = sequences
        self.references = references
        self.asis_calls = asis_calls
        self.receptor = receptor

        # Define parsing blocks
        self.groups = groupby(self.igblast, lambda x: not re.match('# IGBLASTN', x))

    @staticmethod
    def _parseQueryChunk(chunk):
        """
        Parse query section

        Arguments:
          chunk : list of strings

        Returns:
          str : query identifier
        """
        # Extract query id from comments
        query = next((x for x in chunk if x.startswith('# Query:')))

        return query.replace('# Query: ', '', 1)

    @staticmethod
    def _parseSummaryChunk(chunk):
        """
        Parse summary section

        Args:
            chunk: list of strings

        Returns:
            dict : summary section.
        """
        # Mapping for field names in the summary section
        summary_map = {'Top V gene match': 'v_match',
                       'Top D gene match': 'd_match',
                       'Top J gene match': 'j_match',
                       'Chain type': 'chain_type',
                       'stop codon': 'stop_codon',
                       'V-J frame': 'vj_frame',
                       'Productive': 'productive',
                       'Strand': 'strand'}

        # Extract column names from comments
        f = next((x for x in chunk if x.startswith('# V-(D)-J rearrangement summary')))
        f = re.search('summary for query sequence \((.+)\)\.', f).group(1)
        columns = [summary_map[x.strip()] for x in f.split(',')]

        # Extract first row as a list
        row = next((x.split('\t') for x in chunk if not x.startswith('#')))

        # Populate template dictionary with parsed fields
        summary = {v: None for v in summary_map.values()}
        summary.update(dict(zip(columns, row)))

        return summary

    @staticmethod
    def _parseSubregionChunk(chunk):
        """
        Parse CDR3 sequences generated by IgBLAST

        Args:
          chunk: list of strings

        Returns:
          dict : nucleotide and amino acid CDR3 sequences
        """
        # Example:
        #   # Sub-region sequence details (nucleotide sequence, translation, start, end)
        #   CDR3  CAACAGTGGAGTAGTTACCCACGGACG QQWSSYPRT	248	287

        # Define column names
        cdr3_map = {'nucleotide sequence': 'cdr3_igblast',
                    'translation': 'cdr3_igblast_aa',
                    'start': 'cdr3_igblast_start',
                    'end': 'cdr3_igblast_end'}

        # Extract column names from comments
        f = next((x for x in chunk if x.startswith('# Sub-region sequence details')))
        f = re.search('sequence details \((.+)\)', f).group(1)
        columns = [cdr3_map[x.strip()] for x in f.split(',')]

        # Extract first CDR3 as a list and remove the CDR3 label
        rows = next((x.split('\t') for x in chunk if x.startswith('CDR3')))[1:]

        # Populate dictionary with parsed fields
        cdr = {v: None for v in columns}
        cdr.update(dict(zip(columns, rows)))

        # Add length
        if cdr.get('cdr3_igblast', None) is not None:
            cdr['cdr3_igblast_length'] = len(cdr['cdr3_igblast'])

        return cdr

    @staticmethod
    def _parseHitsChunk(chunk):
        """
        Parse hits section

        Args:
          chunk: list of strings

        Returns:
          list: hit table as a list of dictionaries
        """
        # Extract column names from comments
        f = next((x for x in chunk if x.startswith('# Fields:')))
        columns = chain(['segment'], f.replace('# Fields:', '', 1).split(','))
        columns = [x.strip() for x in columns]
        # Split non-comment rows into a list of lists
        rows = [x.split('\t') for x in chunk if not x.startswith('#')]
        # Create list of dictionaries containing hits
        hits = [{k: x[i] for i, k in enumerate(columns)} for x in rows]

        return hits

    # Parse summary results
    @staticmethod
    def _parseSummarySection(summary, db, asis_calls=False):
        """
        Parse summary section

        Arguments:
          summary :  summary section dictionary return by parseBlock
          db : initial database dictionary.
          asis_calls : if True do not parse gene calls for allele names.

        Returns:
          dict : db of results.
        """
        result = {}
        # Parse V, D, and J calls
        if not asis_calls:
            v_call = parseAllele(summary['v_match'], v_allele_regex, action='list')
            d_call = parseAllele(summary['d_match'], d_allele_regex, action='list')
            j_call = parseAllele(summary['j_match'], j_allele_regex, action='list')
            result['v_call'] = ','.join(v_call) if v_call else None
            result['d_call'] = ','.join(d_call) if d_call else None
            result['j_call'] = ','.join(j_call) if j_call else None
        else:
            result['v_call'] = None if summary['v_match'] == 'N/A' else summary['v_match']
            result['d_call'] = None if summary['d_match'] == 'N/A' else summary['d_match']
            result['j_call'] = None if summary['j_match'] == 'N/A' else summary['j_match']

        # Parse locus
        result['locus'] = None if summary['chain_type'] == 'N/A' else summary['chain_type']

        # Parse quality information
        result['stop'] = 'T' if summary['stop_codon'] == 'Yes' else 'F'
        result['in_frame'] = 'T' if summary['vj_frame'] == 'In-frame' else 'F'
        result['functional'] = 'T' if summary['productive'] == 'Yes' else 'F'

        # Reverse complement input sequence if required
        if summary['strand'] == '-':
            seq_rc = Seq(db['sequence_input'], IUPAC.ambiguous_dna).reverse_complement()
            result['sequence_input'] = str(seq_rc)
            result['rev_comp'] = 'T'
        else:
            result['rev_comp'] = 'F'

        return result

    @staticmethod
    def _parseSubregionSection(section, sequence):
        """
        Parse subregion section

        Arguments:
          section :  subregion section dictionary return by parseBlock
          sequence : input sequence

        Returns:
          dict : db of results.
        """
        # Extract junction
        junc_start = int(section['cdr3_igblast_start']) - 3
        junc_end = int(section['cdr3_igblast_end']) + 3
        junc_seq = sequence[(junc_start - 1):junc_end]
        junc_len = len(junc_seq)

        # Translation
        junc_tmp = junc_seq.replace('-', 'N').replace('.', 'N')
        if junc_len % 3 > 0:  junc_tmp = junc_tmp[:junc_len - junc_len % 3]
        junc_aa = str(Seq(junc_tmp).translate())

        # Build return values
        return {'junction': junc_seq,
                'junction_aa': junc_aa,
                'junction_length': junc_len,
                'junction_start': junc_start,
                'junction_end': junc_end}

    @staticmethod
    def _parseVHitPos(v_hit):
        """
        Parse V alignment positions

        Arguments:
          v_hit :  V alignment row from the hit table

        Returns:
          dict: db of D starts and lengths
        """
        result = {}
        # Germline positions
        result['v_germ_start_vdj'] = int(v_hit['s. start'])
        result['v_germ_length_vdj'] = int(v_hit['s. end']) - result['v_germ_start_vdj'] + 1
        # Query sequence positions
        result['v_seq_start'] = int(v_hit['q. start'])
        result['v_seq_length'] = int(v_hit['q. end']) - result['v_seq_start'] + 1
        result['indels'] = 'F' if int(v_hit['gap opens']) == 0 else 'T'

        return result

    @staticmethod
    def _parseDHitPos(d_hit, overlap):
        """
        Parse D alignment positions

        Arguments:
          d_hit :  D alignment row from the hit table
          overlap : V-D overlap length

        Returns:
          dict: db of D starts and lengths
        """
        result = {}
        # Query sequence positions
        result['d_seq_start'] = int(d_hit['q. start']) + overlap
        result['d_seq_length'] = max(int(d_hit['q. end']) - result['d_seq_start'] + 1, 0)
        # Germline positions
        result['d_germ_start'] = int(d_hit['s. start']) + overlap
        result['d_germ_length'] = max(int(d_hit['s. end']) - result['d_germ_start'] + 1, 0)

        return result

    @staticmethod
    def _parseJHitPos(j_hit, overlap):
        """
        Parse J alignment positions

        Arguments:
          j_hit :  J alignment row from the hit table
          overlap : D-J or V-J overlap length

        Returns:
          dict: db of J starts and lengths
        """
        result = {}
        result['j_seq_start'] = int(j_hit['q. start']) + overlap
        result['j_seq_length'] = max(int(j_hit['q. end']) - result['j_seq_start'] + 1, 0)
        result['j_germ_start'] = int(j_hit['s. start']) + overlap
        result['j_germ_length'] = max(int(j_hit['s. end']) - result['j_germ_start'] + 1, 0)

        return result

    @staticmethod
    def _appendSeq(seq, hits, start, trim=True):
        """
        Append aligned query sequence segment

        Arguments:
          seq :  sequence to modify.
          hits : hit table row for the sequence.
          start : start position of the query sequence.
          trim : if True then remove insertions from the hit sequence before appending.

        Returns:
          str: modified sequence.
        """
        if 'subject seq' not in hits or 'query seq' not in hits:
            return None

        # Remove insertions
        if trim:
            for m in re.finditer(r'-', hits['subject seq']):
                ins = m.start()
                seq += hits['query seq'][start:ins]
                start = ins + 1

        # Append
        seq += hits['query seq'][start:]

        return seq

    @staticmethod
    def _parseVHits(hits, db):
        """
        Parse V hit sub-table

        Arguments:
          hits :  hit table as a list of dictionaries.
          db : database dictionary containing summary results.

        Returns:
          dict : db of results.
        """
        result = {}
        seq_vdj = db['sequence_vdj']
        seq_trim = db['sequence_trim']
        v_hit = next(x for x in hits if x['segment'] == 'V')

        # Alignment positions
        result.update(IgBLASTReader._parseVHitPos(v_hit))
        # Update VDJ sequence with and without removing insertions
        result['sequence_vdj'] = IgBLASTReader._appendSeq(seq_vdj, v_hit, 0, trim=False)
        result['sequence_trim'] = IgBLASTReader._appendSeq(seq_trim, v_hit, 0, trim=True)

        return result

    @staticmethod
    def _parseDHits(hits, db):
        """
        Parse D hit sub-table

        Arguments:
          hits :  hit table as a list of dictionaries.
          db : database dictionary containing summary and V results.

        Returns:
          dict : db of results.
        """
        result = {}
        seq_vdj = db['sequence_vdj']
        seq_trim = db['sequence_trim']
        d_hit = next(x for x in hits if x['segment'] == 'D')

        # Determine N-region length and amount of J overlap with V or D alignment
        overlap = 0
        if db['v_call']:
            np1_len = int(d_hit['q. start']) - (db['v_seq_start'] + db['v_seq_length'])
            if np1_len < 0:
                result['np1_length'] = 0
                overlap = abs(np1_len)
            else:
                result['np1_length'] = np1_len
                np1_start = db['v_seq_start'] + db['v_seq_length'] - 1
                np1_end = int(d_hit['q. start']) - 1
                if seq_vdj is not None:
                    seq_vdj += db['sequence_input'][np1_start:np1_end]
                if seq_trim is not None:
                    seq_trim += db['sequence_input'][np1_start:np1_end]

        # D alignment positions
        result.update(IgBLASTReader._parseDHitPos(d_hit, overlap))
        # Update VDJ sequence with and without removing insertions
        result['sequence_vdj'] = IgBLASTReader._appendSeq(seq_vdj, d_hit, overlap, trim=False)
        result['sequence_trim'] = IgBLASTReader._appendSeq(seq_trim, d_hit, overlap, trim=True)

        return result

    @staticmethod
    def _parseJHits(hits, db):
        """
        Parse J hit sub-table

        Arguments:
          hits :  hit table as a list of dictionaries.
          db : database dictionary containing summary, V and D results.

        Returns:
          dict : db of results.
        """
        result = {}
        seq_vdj = db['sequence_vdj']
        seq_trim = db['sequence_trim']
        j_hit = next(x for x in hits if x['segment'] == 'J')

        # Determine N-region length and amount of J overlap with V or D alignment
        overlap = 0
        if db['d_call']:
            np2_len = int(j_hit['q. start']) - (db['d_seq_start'] + db['d_seq_length'])
            if np2_len < 0:
                result['np2_length'] = 0
                overlap = abs(np2_len)
            else:
                result['np2_length'] = np2_len
                n2_start = db['d_seq_start'] + db['d_seq_length'] - 1
                n2_end = int(j_hit['q. start']) - 1
                if seq_vdj is not None:
                    seq_vdj += db['sequence_input'][n2_start:n2_end]
                if seq_trim is not None:
                    seq_trim += db['sequence_input'][n2_start:n2_end]
        elif db['v_call']:
            np1_len = int(j_hit['q. start']) - (db['v_seq_start'] + db['v_seq_length'])
            if np1_len < 0:
                result['np1_length'] = 0
                overlap = abs(np1_len)
            else:
                result['np1_length'] = np1_len
                np1_start = db['v_seq_start'] + db['v_seq_length'] - 1
                np1_end = int(j_hit['q. start']) - 1
                if seq_vdj is not None:
                    seq_vdj += db['sequence_input'][np1_start: np1_end]
                if seq_trim is not None:
                    seq_trim += db['sequence_input'][np1_start: np1_end]
        else:
            result['np1_length'] = 0

        # J alignment positions
        result.update(IgBLASTReader._parseJHitPos(j_hit, overlap))
        # Update VDJ sequence with and without removing insertions
        result['sequence_vdj'] = IgBLASTReader._appendSeq(seq_vdj, j_hit, overlap, trim=False)
        result['sequence_trim'] = IgBLASTReader._appendSeq(seq_trim, j_hit, overlap, trim=True)

        return result

    @staticmethod
    def _parseHitScores(hits, segment):
        """
        Parse alignment scores

        Arguments:
          hits :  hit table as a list of dictionaries.
          segment : segment name; one of 'v', 'd' or 'j'.

        Returns:
          dict : scores
        """
        result = {}
        s_hit = next(x for x in hits if x['segment'] == segment.upper())

        # Score
        try:
            result['%s_score' % segment] = float(s_hit['bit score'])
        except (TypeError, ValueError):
            result['%s_score' % segment] = None
        # Identity
        try:
            result['%s_identity' % segment] = float(s_hit['% identity']) / 100.0
        except (TypeError, ValueError):
            result['%s_identity' % segment] = None
        # E-value
        try:
            result['%s_evalue' % segment] = float(s_hit['evalue'])
        except (TypeError, ValueError):
            result['%s_evalue' % segment] = None
        # BTOP
        try:
            result['%s_btop' % segment] = s_hit['BTOP']
        except (KeyError, TypeError, ValueError):
            result['%s_btop' % segment] = None
        # CIGAR
        try:
            align = decodeBTOP(s_hit['BTOP'])
            align = padAlignment(align, int(s_hit['q. start']) - 1, int(s_hit['s. start']) - 1)
            result['%s_cigar' % segment] = encodeCIGAR(align)
        except (KeyError, TypeError, ValueError):
            result['%s_cigar' % segment] = None

        return result

    def parseBlock(self, block):
        """
        Parses an IgBLAST result into separate sections

        Arguments:
          block (iter): an iterator from itertools.groupby containing a single IgBLAST result.

        Returns:
          dict: a parsed results block;
                with the keys 'query' (sequence identifier as a string),
                'summary' (dictionary of the alignment summary),
                'subregion' (dictionary of IgBLAST CDR3 sequences), and
                'hits' (VDJ hit table as a list of dictionaries).
                Returns None if the block has no data that can be parsed.
        """
        # Parsing info
        #
        #   Columns for non-hit-table sections
        #     'V-(D)-J rearrangement summary': (Top V gene match, Top D gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand)
        #     'V-(D)-J junction details': (V end, V-D junction, D region, D-J junction, J start)
        #     'Alignment summary': (from, to, length, matches, mismatches, gaps, percent identity)
        #     'subregion': (nucleotide sequence, translation)
        #
        #   Ignored sections
        #     'junction': '# V-(D)-J junction details'
        #     'v_alignment': '# Alignment summary'
        #
        #   Hit table fields for -outfmt "7 std qseq sseq btop"
        #     0:  segment
        #     1:  query id
        #     2:  subject id
        #     3:  % identity
        #     4:  alignment length
        #     5:  mismatches
        #     6:  gap opens
        #     7:  gaps
        #     8:  q. start
        #     9:  q. end
        #    10:  s. start
        #    11:  s. end
        #    12:  evalue
        #    13:  bit score
        #    14:  query seq
        #    15:  subject seq
        #    16:  btop
        # Map of valid block parsing keys and functions
        chunk_map = {'query': ('# Query:', self._parseQueryChunk),
                     'summary': ('# V-(D)-J rearrangement summary', self._parseSummaryChunk),
                     'subregion': ('# Sub-region sequence details', self._parseSubregionChunk),
                     'hits': ('# Hit table', self._parseHitsChunk)}

        # Parsing chunks
        results = {}
        for match, chunk in groupby(block, lambda x: x != '\n'):
            if match:
                # Strip whitespace and convert to list
                chunk = [x.strip() for x in chunk]

                # Parse non-query sections
                chunk_dict = {k: f(chunk) for k, (v, f) in chunk_map.items() if chunk[0].startswith(v)}
                results.update(chunk_dict)

        return results if results else None

    def parseSections(self, sections):
        """
        Parses an IgBLAST sections into a db dictionary

        Arguments:
            sections : dictionary of parsed sections from parseBlock.

        Returns:
          dict : db entries.
        """
        # Initialize dictionary with input sequence and id
        db = {}
        if 'query' in sections:
            query = sections['query']
            db['sequence_id'] = query
            db['sequence_input'] = str(self.sequences[query].seq)

        # Parse summary section
        if 'summary' in sections:
            db.update(self._parseSummarySection(sections['summary'], db, asis_calls=self.asis_calls))

        # Parse hit table
        if 'hits' in sections:
            db['sequence_vdj'] = ''
            db['sequence_trim'] = ''
            if db['v_call']:
                db.update(self._parseVHits(sections['hits'], db))
                db.update(self._parseHitScores(sections['hits'], 'v'))
            if db['d_call']:
                db.update(self._parseDHits(sections['hits'], db))
                db.update(self._parseHitScores(sections['hits'], 'd'))
            if db['j_call']:
                db.update(self._parseJHits(sections['hits'], db))
                db.update(self._parseHitScores(sections['hits'], 'j'))

        # Create IMGT-gapped sequence
        if ('v_call' in db and db['v_call']) and ('sequence_trim' in db and db['sequence_trim']):
            try:
                imgt_dict = gapV(db['sequence_trim'],
                                 v_germ_start=db['v_germ_start_vdj'],
                                 v_germ_length=db['v_germ_length_vdj'],
                                 v_call=db['v_call'],
                                 references=self.references,
                                 asis_calls=self.asis_calls)
            except KeyError as e:
                imgt_dict = {'sequence_imgt': None,
                             'v_germ_start_imgt': None,
                             'v_germ_length_imgt': None}
                printWarning(e)
            db.update(imgt_dict)
            del db['sequence_trim']

        # Add junction
        if 'subregion' in sections and 'cdr3_igblast_start' in sections['subregion']:
            junc_dict = self._parseSubregionSection(sections['subregion'], db['sequence_input'])
            db.update(junc_dict)
        elif ('j_call' in db and db['j_call']) and ('sequence_imgt' in db and db['sequence_imgt']):
            junc_dict = inferJunction(db['sequence_imgt'],
                                      j_germ_start=db['j_germ_start'],
                                      j_germ_length=db['j_germ_length'],
                                      j_call=db['j_call'],
                                      references=self.references,
                                      asis_calls=self.asis_calls)
            db.update(junc_dict)

        # Add IgBLAST CDR3 sequences
        if 'subregion' in sections:
            # Sequences already parsed into dict by parseBlock
            db.update(sections['subregion'])
        else:
            # Section does not exist (ie, older version of IgBLAST or CDR3 not found)
            db.update({'cdr3_igblast': None, 'cdr3_igblast_aa': None})

        # Add FWR and CDR regions
        db.update(getRegions(db.get('sequence_imgt', None), db.get('junction_length', None)))

        return db

    def __iter__(self):
        """
        Iterator initializer.

        Returns:
          changeo.IO.IgBLASTReader
        """
        return self

    def __next__(self):
        """
        Next method.

        Returns:
          changeo.Receptor.Receptor : parsed IMGT/HighV-QUEST result as an Receptor (receptor=True) or dictionary (receptor=False).
        """
        # Get next block from groups iterator
        try:
            match = False
            block = None
            while not match:
                match, block = next(self.groups)
        except StopIteration:
            raise StopIteration

        # Parse block
        sections = self.parseBlock(block)
        db = self.parseSections(sections)

        if self.receptor:
            return Receptor(db)
        else:
            return db


class IHMMuneReader:
    """
    An iterator to read and parse iHMMune-Align output files.
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
    ihmmune_fields = ['SEQUENCE_ID',
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
                      'V_SEQ_START',
                      'V_SEQ_LENGTH',
                      'A_SCORE']

    # Ordered list of known fields
    @staticmethod
    def customFields(scores=False, regions=False, schema=None):
        """
        Returns non-standard Receptor attributes defined by the parser

        Arguments:
          scores : if True include alignment scoring fields.
          regions : if True include IMGT-gapped CDR and FWR region fields.
          schema : schema class to pass field through for conversion.
                   If None, return changeo.Receptor.Receptor attribute names.

        Returns:
          list : list of field names.
        """
        # Alignment scoring fields
        score_fields = ['vdj_score']

        # FWR amd CDR fields
        region_fields = ['fwr1_imgt',
                         'fwr2_imgt',
                         'fwr3_imgt',
                         'fwr4_imgt',
                         'cdr1_imgt',
                         'cdr2_imgt',
                         'cdr3_imgt']

        fields = []
        if scores:  fields.extend(score_fields)
        if regions:  fields.extend(region_fields)

        # Convert field names if schema provided
        if schema is not None:
            fields = [schema.fromReceptor(f) for f in fields]

        return fields

    def __init__(self, ihmmune, sequences, references, receptor=True):
        """
        Initializer

        Arguments:
          ihmmune (file): handle to an open iHMMune-Align output file.
          sequences (dict): dictionary with sequence descriptions as keys mapping to the SeqRecord containing
                            the original query sequences.
          references (dict): dictionary of IMGT gapped germline sequences.
          receptor (bool): if True (default) iteration returns an Receptor object, otherwise it returns a dictionary

        Returns:
          changeo.IO.IHMMuneReader
        """
        # Arguments
        self.ihmmune = ihmmune
        self.sequences = sequences
        self.references = references
        self.receptor = receptor

        # Open reader
        self.records = csv.DictReader(self.ihmmune, fieldnames=IHMMuneReader.ihmmune_fields,
                                      delimiter=';', quotechar='"')

    @staticmethod
    def _parseFunctionality(record):
        """
        Parse functionality information

        Arguments:
          record : dictionary containing a single row from the iHMMune-Align ouptut.

        Returns:
          dict : database entries containing functionality information.
        """
        # Orientation
        def _revcomp():
            return 'F' if int(record['RC']) == 0 else 'T'

        # Functional
        def _functional():
            if not record['V_CALL'] or \
                    record['V_CALL'].startswith('NA - ') or \
                    record['J_INFRAME'] != 'true' or \
                    not record['J_CALL'] or \
                    record['J_CALL'] == 'NO_JGENE_ALIGNMENT' or \
                    int(record['STOP_COUNT']) > 0:
                return 'F'
            else:
                return 'T'

        # Stop codon
        def _stop():
            return 'T' if int(record['STOP_COUNT']) > 0 else 'F'

        # J in-frame
        def _inframe():
            return 'T' if record['J_INFRAME'] == 'true' else 'F'

        # Indels
        def _indels():
            check = [x is not None and 'indels' in x \
                     for x in [record['V_CALL'], record['D_CALL'], record['J_CALL']]]
            return 'T' if any(check) else 'F'

        # Parse functionality
        result = {'rev_comp': _revcomp(),
                  'functional': _functional(),
                  'in_frame': _inframe(),
                  'stop': _stop(),
                  'indels':  _indels()}

        return result

    @staticmethod
    def _parseGenes(record):
        """
        Parse gene calls

        Arguments:
          record : dictionary containing a single row from the iHMMune-Align ouptut.

        Returns:
          dict : database entries for gene calls.
        """
        # Extract allele calls
        v_call = parseAllele(record['V_CALL'], v_allele_regex, action='list')
        d_call = parseAllele(record['D_CALL'], d_allele_regex, action='list')
        j_call = parseAllele(record['J_CALL'], j_allele_regex, action='list')

        # Build return object
        result = {'v_call': ','.join(v_call) if v_call else None,
                  'd_call': ','.join(d_call) if d_call else None,
                  'j_call': ','.join(j_call) if j_call else None}

        return result

    @staticmethod
    def _parseNPHit(record):
        """
        Parse N/P region alignment information

        Arguments:
          record : dictionary containing a single row from the iHMMune-Align ouptut.

        Returns:
          dict : database entries containing N/P region lengths.
        """
        # N/P lengths
        result = {'np1_length': len(record['NP1_SEQ']),
                  'np2_length': len(record['NP2_SEQ'])}

        return result

    @staticmethod
    def _parseVHit(record, db):
        """
        Parse V alignment information

        Arguments:
          record : dictionary containing a single row from the iHMMune-Align ouptut.
          db : database containing V and D alignment information.

        Returns:
          dict : database entries containing V call and alignment positions.
        """
        # Default return
        result = {'v_seq_start': None,
                  'v_seq_length': None,
                  'v_germ_start_vdj': None,
                  'v_germ_length_vdj': None}

        # Find V positions
        if db['v_call']:
            # Query positions
            result['v_seq_start'] = int(record['V_SEQ_START'])
            result['v_seq_length'] = len(record['V_SEQ'].strip('.'))
            # Germline positions
            result['v_germ_start_vdj'] = 1
            result['v_germ_length_vdj'] = result['v_seq_length']

        return result

    def _parseDHit(record, db):
        """
        Parse D alignment information

        Arguments:
          record : dictionary containing a single row from the iHMMune-Align ouptut.
          db : database containing V alignment information.


        Returns:
          dict : database entries containing D call and alignment positions.
        """

        # D start position
        def _dstart():
            nb = [db['v_seq_start'],
                  db['v_seq_length'],
                  db['np1_length']]
            return sum(int(i) for i in nb if i)

        # Default return
        result = {'d_seq_start': None,
                  'd_seq_length': None,
                  'd_germ_start': None,
                  'd_germ_length': None}

        if db['d_call']:
            # Query positions
            result['d_seq_start'] = _dstart()
            result['d_seq_length'] = len(record['D_SEQ'].strip('.'))
            # Germline positions
            result['d_germ_start'] = len(record['D_SEQ']) - len(record['D_SEQ'].lstrip('.'))
            result['d_germ_length'] = result['d_seq_length']

        return result

    @staticmethod
    def _parseJHit(record, db):
        """
        Parse J alignment information

        Arguments:
          record : dictionary containing a single row from the iHMMune-Align ouptut.
          db : database containing V and D alignment information.

        Returns:
          dict : database entries containing J call and alignment positions.
        """

        # J start position
        def _jstart():
            # J positions
            nb = [db['v_seq_start'],
                  db['v_seq_length'],
                  db['np1_length'],
                  db['d_seq_length'],
                  db['np2_length']]
            return sum(int(i) for i in nb if i)

        # Default return
        result = {'j_seq_start': None,
                  'j_seq_length': None,
                  'j_germ_start': None,
                  'j_germ_length': None}

        # Find J region
        if db['j_call']:
            # Query positions
            result['j_seq_start'] = _jstart()
            result['j_seq_length'] = len(record['J_SEQ'].strip('.'))
            # Germline positions
            result['j_germ_start'] = len(record['J_SEQ']) - len(record['J_SEQ'].lstrip('.'))
            result['j_germ_length'] = result['j_seq_length']

        return result

    @staticmethod
    def _assembleVDJ(record, db):
        """
        Build full length V(D)J sequence

        Arguments:
          record : dictionary containing a single row from the iHMMune-Align ouptut.
          db : database containing V and D alignment information.

        Returns:
          dict : database entries containing the full length V(D)J sequence.
        """
        segments = [record['V_SEQ'].strip('.') if db['v_call'] else '',
                    record['NP1_SEQ'] if db['np1_length'] else '',
                    record['D_SEQ'].strip('.') if db['d_call'] else '',
                    record['NP2_SEQ'] if db['np2_length'] else '',
                    record['J_SEQ'].strip('.') if db['j_call'] else '']

        return {'sequence_vdj': ''.join(segments)}

    @staticmethod
    def _parseScores(record):
        """
        Parse alignment scores

        Arguments:
          record : dictionary containing a single row from the iHMMune-Align ouptut.

        Returns:
          dict : database entries for alignment scores.
        """
        result = {}
        try:
            result['vdj_score'] = float(record['HMM_SCORE'])
        except (TypeError, ValueError):
            result['vdj_score'] = None

        return result

    def parseRecord(self, record):
        """
        Parses a single row from each IMTG file.

        Arguments:
          record : dictionary containing one row of iHMMune-Align file.

        Returns:
          dict : database entry for the row.
        """
        # Extract query ID and sequence
        query = record['SEQUENCE_ID']
        db = {'sequence_id': query,
              'sequence_input': str(self.sequences[query].seq)}

        # Check for valid alignment
        if not record['V_CALL'] or \
                record['V_CALL'].startswith('NA - ') or \
                record['V_CALL'].startswith('State path'):
            db['functional'] = None
            db['v_call'] = None
            db['d_call'] = None
            db['j_call'] = None
            return db

        # Parse record
        db.update(IHMMuneReader._parseFunctionality(record))
        db.update(IHMMuneReader._parseGenes(record))
        db.update(IHMMuneReader._parseNPHit(record))
        db.update(IHMMuneReader._parseVHit(record, db))
        db.update(IHMMuneReader._parseDHit(record, db))
        db.update(IHMMuneReader._parseJHit(record, db))
        db.update(IHMMuneReader._assembleVDJ(record, db))

        # Create IMGT-gapped sequence
        if 'v_call' in db and db['v_call'] and 'sequence_vdj' in db and db['sequence_vdj']:
            try:
                imgt_dict = gapV(db['sequence_vdj'],
                                 v_germ_start=db['v_germ_start_vdj'],
                                 v_germ_length=db['v_germ_length_vdj'],
                                 v_call=db['v_call'],
                                 references=self.references)
            except KeyError as e:
                imgt_dict = {'sequence_imgt': None,
                             'v_germ_start_imgt': None,
                             'v_germ_length_imgt': None}
                printWarning(e)
            db.update(imgt_dict)

        # Infer IMGT junction
        if ('j_call' in db and db['j_call']) and ('sequence_imgt' in db and db['sequence_imgt']):
            junc_dict = inferJunction(db['sequence_imgt'],
                                      j_germ_start=db['j_germ_start'],
                                      j_germ_length=db['j_germ_length'],
                                      j_call=db['j_call'],
                                      references=self.references)
            db.update(junc_dict)

        # Overall alignment score
        db.update(IHMMuneReader._parseScores(record))

        # FWR and CDR regions
        db.update(getRegions(db.get('sequence_imgt', None), db.get('junction_length', None)))
        
        return db

    def __iter__(self):
        """
        Iterator initializer.

        Returns:
          changeo.IO.IHMMuneReader
        """
        return self

    def __next__(self):
        """
        Next method.

        Returns:
          changeo.Receptor.Receptor : parsed IMGT/HighV-QUEST result as an Receptor (receptor=True) or dictionary (receptor=False).
        """
        # Get next set of records from dictionary readers
        try:
            record = None
            while not record:
                record = next(self.records)
        except StopIteration:
            raise StopIteration

        db = self.parseRecord(record)

        if self.receptor:
            return Receptor(db)
        else:
            return db


def readGermlines(references, asis=False):
    """
    Parses germline repositories

    Arguments:
      references (list): list of strings specifying directories and/or files from which to read germline records.
      asis (bool): if True use sequence ID as record name and do not parse headers for allele names.

    Returns:
      dict: Dictionary of germlines in the form {allele: sequence}.
    """
    repo_files = []
    # Iterate over items passed to commandline
    for r in references:
        if os.path.isdir(r):
            # If directory, get fasta files from within
            repo_files.extend([os.path.join(r, f) for f in os.listdir(r) \
                          if getFileType(f) == 'fasta'])
        elif os.path.isfile(r) and getFileType(r) == 'fasta':
            # If file, make sure file is fasta
            repo_files.extend([r])

    # Catch instances where no valid fasta files were passed in
    if len(repo_files) < 1:
        printError('No valid germline fasta files (.fasta, .fna, .fa) were found at %s.' % ','.join(references))

    repo_dict = {}
    for file_name in repo_files:
        with open(file_name, 'rU') as file_handle:
            germlines = SeqIO.parse(file_handle, 'fasta')
            for g in germlines:
                germ_key = parseAllele(g.description, allele_regex, 'first') if not asis else g.id
                repo_dict[germ_key] = str(g.seq).upper()

    return repo_dict


def extractIMGT(imgt_output):
    """
    Extract necessary files from IMGT/HighV-QUEST results.

    Arguments:
      imgt_output : zipped file or unzipped folder output by IMGT/HighV-QUEST.

    Returns:
      tuple : (temporary directory handle, dictionary with names of extracted IMGT files).
    """
    # Map of IMGT file names
    imgt_names = ('1_Summary', '2_IMGT-gapped', '3_Nt-sequences', '6_Junction')
    imgt_keys = ('summary', 'gapped', 'ntseq', 'junction')

    # Open temporary directory and intialize return dictionary
    temp_dir = TemporaryDirectory()

    # Zip input
    if zipfile.is_zipfile(imgt_output):
        imgt_zip = zipfile.ZipFile(imgt_output, 'r')
        # Extract required files
        imgt_files = sorted([n for n in imgt_zip.namelist() \
                             if os.path.basename(n).startswith(imgt_names)])
        imgt_zip.extractall(temp_dir.name, imgt_files)
        # Define file dictionary
        imgt_dict = {k: os.path.join(temp_dir.name, f) for k, f in zip_longest(imgt_keys, imgt_files)}
    # Folder input
    elif os.path.isdir(imgt_output):
        folder_files = []
        for root, dirs, files in os.walk(imgt_output):
            folder_files.extend([os.path.join(os.path.abspath(root), f) for f in files])
        # Define file dictionary
        imgt_files = sorted([n for n in folder_files \
                             if os.path.basename(n).startswith(imgt_names)])
        imgt_dict = {k: f for k, f in zip_longest(imgt_keys, imgt_files)}
    # Tarball input
    elif tarfile.is_tarfile(imgt_output):
        imgt_tar = tarfile.open(imgt_output, 'r')
        # Extract required files
        imgt_files = sorted([n for n in imgt_tar.getnames() \
                             if os.path.basename(n).startswith(imgt_names)])
        imgt_tar.extractall(temp_dir.name, [imgt_tar.getmember(n) for n in imgt_files])
        # Define file dictionary
        imgt_dict = {k: os.path.join(temp_dir.name, f) for k, f in zip_longest(imgt_keys, imgt_files)}
    else:
        printError('Unsupported IGMT output file. Must be either a zipped file (.zip), LZMA compressed tarfile (.txz) or a folder.')

    # Check extraction for errors
    if len(imgt_dict) != len(imgt_names):
        printError('Extra files or missing necessary file IMGT output %s.' % imgt_output)

    return temp_dir, imgt_dict


def countDbFile(file):
    """
    Counts the records in database files

    Arguments:
      file : tab-delimited database file.

    Returns:
      int : count of records in the database file.
    """
    # Count records and check file
    try:
        with open(file, 'rt') as db_handle:
            db_records = csv.reader(db_handle, dialect='excel-tab')
            for i, __ in enumerate(db_records):  pass
        db_count = i
    except IOError:
        printError('File %s cannot be read.' % file)
    except:
        printError('File %s is invalid.' % file)
    else:
        if db_count == 0:  printError('File %s is empty.' % file)

    return db_count


def getDbFields(file, add=None, exclude=None, reader=TSVReader):
    """
    Get field names from a db file

    Arguments:
      file : db file to pull base fields from.
      add : fields to append to the field set.
      exclude : fields to exclude from the field set.
      reader : reader class.

    Returns:
        list : list of field names
    """
    try:
        with open(file, 'rt') as handle:
            fields = reader(handle).fields
    except IOError:
        printError('File %s cannot be read.' % file)
    except:
        printError('File %s is invalid.' % file)

    # Add extra fields
    if add is not None:
        if not isinstance(add, list):  add = [add]
        fields.extend([f for f in add if f not in fields])
    # Remove unwanted fields
    if exclude is not None:
        if not isinstance(exclude, list):  exclude = [exclude]
        fields = [f for f in fields if f not in exclude]

    return fields


def getFormatOperators(format):
    """
    Simple wrapper for fetching the set of operator classes for a data format

    Arguments:
      format (str): name of the data format.

    Returns:
      tuple: a tuple with the reader class, writer class, and schema definition class.
    """
    # Format options
    if format == 'changeo':
        reader = ChangeoReader
        writer = ChangeoWriter
        schema = ChangeoSchema
    elif format == 'airr':
        reader = AIRRReader
        writer = AIRRWriter
        schema = AIRRSchema
    else:
        raise ValueError

    return reader, writer, schema


def splitName(file):
    """
    Extract the extension from a file name

    Arguments:
      file (str): file name.

    Returns:
      tuple : tuple of the file directory, basename and extension.
    """
    directory, filename = os.path.split(file)
    basename, extension = os.path.splitext(filename)
    extension = extension.lower().lstrip('.')

    return directory, basename, extension


def getOutputName(file, out_label=None, out_dir=None, out_name=None, out_type=None):
    """
    Creates and output filename from an existing filename

    Arguments:
      file : filename to base output file name on.
      out_label : text to be inserted before the file extension;
                  if None do not add a label.
      out_type : the file extension of the output file;
                 if None use input file extension.
      out_dir : the output directory;
                if None use directory of input file
      out_name : the short filename to use for the output file;
                 if None use input file short name.

    Returns:
      str: file name.
    """
    # Get filename components
    directory, basename, extension = splitName(file)

    # Define output directory
    if out_dir is None:
        out_dir = directory
    else:
        out_dir = os.path.abspath(out_dir)
        if not os.path.exists(out_dir):  os.mkdir(out_dir)
    # Define output file prefix
    if out_name is None:  out_name = basename
    # Define output file extension
    if out_type is None:  out_type = extension

    # Define output file name
    if out_label is None:
        out_file = os.path.join(out_dir, '%s.%s' % (out_name, out_type))
    else:
        out_file = os.path.join(out_dir, '%s_%s.%s' % (out_name, out_label, out_type))

    # Return file name
    return out_file


def getOutputHandle(file, out_label=None, out_dir=None, out_name=None, out_type=None):
    """
    Opens an output file handle

    Arguments:
      file : filename to base output file name on.
      out_label : text to be inserted before the file extension;
                  if None do not add a label.
      out_type : the file extension of the output file;
                 if None use input file extension.
      out_dir : the output directory;
                if None use directory of input file
      out_name : the short filename to use for the output file;
                 if None use input file short name.

    Returns:
      file : File handle
    """
    out_file = getOutputName(file, out_label=out_label, out_dir=out_dir,
                             out_name=out_name, out_type=out_type)

    # Open and return handle
    try:
        return open(out_file, mode='w')
    except:
        printError('File %s cannot be opened.' % out_file)


def checkFields(attributes, header, schema=ChangeoSchema):
    """
    Checks that a file header contains a required set of Receptor attributes

    Arguments:
        attributes (list): list of Receptor attributes to check for.
        header (list): list of fields names in the file header.
        schema (object): schema object to convert field names to Receptor attributes.

    Returns:
        bool: True if all attributes mapping fields are found.

    Raises:
        LookupError:
    """
    columns = [schema.fromReceptor(f) for f in attributes]
    missing = [x for x in columns if x not in header]

    if len(missing) > 0:
        raise LookupError('Missing required fields: %s' % ', '.join(missing))

    return True


def yamlDict(file):
    """
    Returns a dictionary from a yaml file

    Arguments:
      file (str): simple yaml file with rows in the form 'argument: value'.

    Returns:
      dict: dictionary of key:value pairs in the file.
    """
    try:
        yaml_dict = dict(yaml.load(open(file, 'r')))
    except:
        printError('YAML file is invalid.')

    return yaml_dict