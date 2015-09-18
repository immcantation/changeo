"""
Multiprocessing functions
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from changeo import __version__, __date__

# Imports
import os
import sys
from collections import OrderedDict
from time import time

# Presto and changeo imports
from presto.IO import getFileType, getOutputHandle, printProgress, printLog
from changeo.IO import readDbFile, countDbFile, getDbWriter
from changeo.Receptor import IgRecord


class DbData:
    """
    A class defining IgRecord data objects for worker processes
    """
    # Instantiation
    def __init__(self, key, records):
        self.id = key
        self.data = records
        self.valid = (key is not None and records is not None)

    # Boolean evaluation
    def __bool__(self):
        return self.valid

    # Length evaluation
    def __len__(self):
        if isinstance(self.data, IgRecord):
            return 1
        elif self.data is None:
            return 0
        else:
            return len(self.data)


class DbResult:
    """
    A class defining IgRecord result objects for collector processes
    """
    # Instantiation
    def __init__(self, key, records):
        self.id = key
        self.data = records
        self.results = None
        self.valid = False
        self.log = OrderedDict([('ID', key)])
        #if isinstance(values, list):
        #    for v in values:  setattr(self, v, None)
        #else:
        #    setattr(self, values, None)

    # Boolean evaluation
    def __bool__(self):
        return self.valid

    # Length evaluation
    def __len__(self):
        if isinstance(self.results, IgRecord):
            return 1
        elif self.results is None:
            return 0
        else:
            return len(self.results)

    # Set data_count to number of data records
    @property
    def data_count(self):
        if isinstance(self.data, IgRecord):
            return 1
        elif self.data is None:
            return 0
        else:
            return len(self.data)


def feedDbQueue(alive, data_queue, db_file, group_func=None, group_args={}):
    """
    Feeds the data queue with Ig records

    Arguments:
      alive : multiprocessing.Value boolean controlling whether processing continues
              if False exit process
      data_queue : multiprocessing.Queue to hold data for processing
      db_file : Database file
      group_func : Function to use for grouping records
      group_args : Dictionary of arguments to pass to group_func

    Returns:
      None
    """
    # Open input file and perform grouping
    try:
        # Iterate over Ig records and assign groups
        db_iter = readDbFile(db_file)
        if group_func is not None:
            group_dict = group_func(db_iter, **group_args)
            group_iter = iter(group_dict.items())
        else:
            group_iter = ((r.id, r) for r in db_iter)
    except:
        alive.value = False
        raise

    # Add groups to data queue
    try:
        # Iterate over groups and feed data queue
        while alive.value:
            # Get data from queue
            if data_queue.full():  continue
            else:  data = next(group_iter, None)
            # Exit upon reaching end of iterator
            if data is None:  break

            # Feed queue
            data_queue.put(DbData(*data))
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
    except:
        #sys.stderr.write('Exception in feeder queue feeding step\n')
        alive.value = False
        raise

    return None


def processDbQueue(alive, data_queue, result_queue, process_func, process_args={}):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments:
      alive : multiprocessing.Value boolean controlling whether processing
            continues; when False function returns
      data_queue : multiprocessing.Queue holding data to process
      result_queue : multiprocessing.Queue to hold processed results
      process_func : Function to use for filtering sequences
      process_args : Dictionary of arguments to pass to process_func

    Returns:
      None
    """
    try:
        # Iterator over data queue until sentinel object reached
        while alive.value:
            # Get data from queue
            if data_queue.empty():  continue
            else:  data = data_queue.get()
            # Exit upon reaching sentinel
            if data is None:  break

            # Perform work
            result = process_func(data, **process_args)

            # Feed results to result queue
            result_queue.put(result)
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
    except:
        alive.value = False
        sys.stderr.write('Error processing sequence with ID: %s.\n' % str(data.id))
        raise

    return None


def collectDbQueue(alive, result_queue, collect_queue, db_file, task_label, out_args,
                   add_fields=None):
    """
    Pulls from results queue, assembles results and manages log and file IO

    Arguments:
      alive : multiprocessing.Value boolean controlling whether processing
              continues; when False function returns
      result_queue : multiprocessing.Queue holding worker results
      collect_queue : multiprocessing.Queue to store collector return values
      db_file : Database file name
      task_label : Task label used to tag the output files
      out_args : Common output argument dictionary from parseCommonArgs
      add_fields : List of fields added to the writer not present in the in_file;
                 if None do not add fields

    Returns:
      None : Adds a dictionary with key value pairs to collect_queue containing
            'log' defining a log object,
            'out_files' defining the output file names
    """
    try:
        result_count = countDbFile(db_file)

        # Define output format
        out_type = getFileType(db_file) if out_args['out_type'] is None \
                   else out_args['out_type']

        # Defined valid alignment output handle
        pass_handle = getOutputHandle(db_file,
                                      '%s-pass' % task_label,
                                      out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'],
                                      out_type=out_type)
        pass_writer = getDbWriter(pass_handle, db_file, add_fields=add_fields)
        # Defined failed alignment output handle
        if out_args['failed']:
            fail_handle = getOutputHandle(db_file,
                                          '%s-fail'  % task_label,
                                          out_dir=out_args['out_dir'],
                                          out_name=out_args['out_name'],
                                          out_type=out_type)
            fail_writer = getDbWriter(fail_handle, db_file)
        else:
            fail_handle = None

        # Define log handle
        if out_args['log_file'] is None:
            log_handle = None
        else:
            log_handle = open(out_args['log_file'], 'w')
    except:
        alive.value = False
        raise

    try:
        # Iterator over results queue until sentinel object reached
        start_time = time()
        set_count = rec_count = pass_count = fail_count = 0
        while alive.value:
            # Get result from queue
            if result_queue.empty():  continue
            else:  result = result_queue.get()
            # Exit upon reaching sentinel
            if result is None:  break

            # Print progress for previous iteration
            printProgress(pass_count, result_count, 0.05, start_time)

            # Update counts for current iteration
            set_count += 1
            rec_count += result.data_count

            # Write log
            printLog(result.log, handle=log_handle)

            # Write alignments
            if result:
                pass_count += result.data_count
                for rec in result.results:
                    pass_writer.writerow(rec.toDict())
            else:
                fail_count += result.data_count
                if fail_handle is not None:
                    for rec in result.data:
                        pass_writer.writerow(rec.toDict())
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None

        # Print total counts
        printProgress(pass_count, result_count, 0.05, start_time)

        # Update return values
        log = OrderedDict()
        log['OUTPUT'] = os.path.basename(pass_handle.name)
        log['RECORDS'] = rec_count
        log['GROUPS'] = set_count
        log['PASS'] = pass_count
        log['FAIL'] = fail_count
        collect_dict = {'log':log, 'out_files': [pass_handle.name]}
        collect_queue.put(collect_dict)

        # Close file handles
        pass_handle.close()
        if fail_handle is not None:  fail_handle.close()
        if log_handle is not None:  log_handle.close()
    except:
        alive.value = False
        raise

    return None