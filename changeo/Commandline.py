"""
Commandline interface
"""

# Info
__author__    = 'Jason Anthony Vander Heiden, Namita Gupta'
from changeo import __version__, __date__

# Imports
import os
import sys
import multiprocessing as mp
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, \
                     RawDescriptionHelpFormatter

# Changeo imports
from presto.IO import printWarning, printError
from changeo.Defaults import choices_format, default_format
from changeo.IO import gzipOutputName
from changeo.Receptor import AIRRSchema, ChangeoSchema


class CommonHelpFormatter(RawDescriptionHelpFormatter, ArgumentDefaultsHelpFormatter):
    """
    Custom argparse.HelpFormatter
    """
    # TODO:  add some sort of list formating for arguments with choices
    # TODO:  override argument position.
    # def __init__(self, prog, indent_increment=2, max_help_position=10, width=None):
    #    super(self.__class__, self).__init__(self, prog,
    #                                         indent_increment=indent_increment,
    #                                         max_help_position=max_help_position,
    #                                         width=width)

    # TODO:  remove multiple inheritance and clean up default value printing.
    # From ArgumentDefaultsHelpFormatter
    # def _get_help_string(self, action):
    #     help = action.help
    #     if '%(default)' not in action.help:
    #         if action.default is not SUPPRESS:
    #             defaulting_nargs = [OPTIONAL, ZERO_OR_MORE]
    #             if action.option_strings or action.nargs in defaulting_nargs:
    #                 help += ' (default: %(default)s)'
    #     return help
    pass


def getCommonArgParser(db_in=True, db_out=True, out_file=True, failed=True, log=True,
                       format=True, multiproc=False, add_help=True, gzip_output=True):
    """
    Defines an ArgumentParser object with common pRESTO arguments

    Arguments:
      db_in (bool): if True include tab delimited database input arguments.
      db_out (bool): if True include explicit output file name argument.
      out_file (bool): if True add explicit output file name arguments.
      failed (bool): if True include arguments for output of failed results.
      log (bool): if True include log arguments.
      format (bool): input and output type arguments.
      multiproc (bool): if True include multiprocessing arguments.
      gzip_output (bool): if True include the --gzip-output argument. Should be
                          False for tools whose output is not produced through
                          getOutputHandle/openFile (e.g. wrapped external
                          executables), since the flag would otherwise be
                          silently ignored.

    Returns:
      argparse.ArgumentParser : an argument parser.
    """
    parser = ArgumentParser(formatter_class=CommonHelpFormatter, add_help=False)

    # Add help and version arguments
    if add_help:
        group_help = parser.add_argument_group('help')
        group_help.add_argument('--version', action='version',
                                version='%(prog)s:' + ' %s %s' %(__version__, __date__))
        group_help.add_argument('-h', '--help', action='help', help='show this help message and exit')

    # Set standard group
    group = parser.add_argument_group('standard arguments')

    # Database arguments
    if db_in:
        group.add_argument('-d', nargs='+', action='store', dest='db_files', required=True,
                           help='A list of tab delimited database files.')
    if db_out:
        # Place holder for the future
        pass

    # Output filename
    if out_file:
        group.add_argument('-o', nargs='+', action='store', dest='out_files', default=None,
                           help='''Explicit output file name. Note, this argument cannot be used with 
                                 the --failed, --outdir, or --outname arguments. If unspecified, then
                                 the output filename will be based on the input filename(s).''')

    # Universal arguments
    group.add_argument('--outdir', action='store', dest='out_dir', default=None,
                       help='''Specify to changes the output directory to the location specified.
                            The input file directory is used if this is not specified.''')
    group.add_argument('--outname', action='store', dest='out_name', default=None,
                       help='''Changes the prefix of the successfully processed output file
                            to the string specified. May not be specified with multiple
                            input files.''')
    if gzip_output:
        group.add_argument('--gzip-output', action='store_true', dest='gzip_output', default=None,
                           help='''Specify to force gzip compressed output files. Output is
                                automatically gzip compressed when the input file(s) are
                                themselves gzip compressed (.gz), even without this argument.''')

    # Log arguments
    if log:
        group.add_argument('--log', action='store', dest='log_file', default=None,
                           help='''Specify to write verbose logging to a file. May not be
                                specified with multiple input files.''')

    # Failed result arguments
    if failed:
        group.add_argument('--failed', action='store_true', dest='failed',
                           help='''If specified create files containing records that
                                fail processing.''')
    # Format arguments
    if format:
        group.add_argument('--format', action='store', dest='format',
                           default=default_format, choices=choices_format,
                           help='''Output format. Also specifies the input format for tools accepting 
                                tab delimited AIRR Rearrangement or Change-O files.''')

    # Multiprocessing arguments
    if multiproc:
        group.add_argument('--nproc', action='store', dest='nproc', type=int, default=mp.cpu_count(),
                           help='''The number of simultaneous computational processes to execute
                                (CPU cores to utilize).''')

    return parser


def _samePath(path_a, path_b):
    """
    Determines whether two paths refer to the same file

    Arguments:
      path_a (str): first file path.
      path_b (str): second file path.

    Returns:
      bool: True if the paths resolve to the same file. Uses os.path.samefile
           (which follows symlinks and normalizes case/aliasing per the
           filesystem) when both paths exist; otherwise falls back to
           comparing normalized absolute paths, since one of the paths
           (typically the output) may not exist yet.
    """
    if os.path.exists(path_a) and os.path.exists(path_b):
        try:
            return os.path.samefile(path_a, path_b)
        except OSError:
            pass

    return os.path.normpath(os.path.abspath(path_a)) == os.path.normpath(os.path.abspath(path_b))


def parseCommonArgs(args, in_arg=None, in_types=None, in_list=False):
    """
    Checks common arguments from getCommonArgParser and transforms output options to a dictionary

    Arguments: 
      args : Argument Namespace defined by ArgumentParser.parse_args.
      in_arg : String defining a non-standard input file argument to verify;
               by default 'db_files' and 'seq_files' are supported in that order.
      in_types : List of types (file extensions as strings) to allow for files in file_arg;
                 if None do not check type.
      in_list : if True allow multiple input files with the out_name and log arguments.
                    
    Returns:
      dict : Dictionary copy of args with output arguments embedded in the dictionary out_args
    """ 
    db_types = ['.tab', '.tsv', '.txt']
    seq_types = ['.fasta', 'fna', '.fa', '.fastq', '.fq']
    if in_types is not None:
        in_types = [f.lower for f in in_types]
    args_dict = args.__dict__.copy()
    
    # Count input files
    if 'seq_files' in args_dict:
        input_count = len(args_dict['seq_files']  or [])
        input_files = args_dict['seq_files']
    elif 'db_files' in args_dict:
        input_count = len(args_dict['db_files'] or [])
        input_files = args_dict['db_files']
    elif in_arg is not None and in_arg in args_dict:
        input_count = len(args_dict[in_arg] or [])
        input_files = args_dict[in_arg]
    else:
        printError('Cannot determine input file argument.')

    # Verify sequence files
    if 'seq_files' in args_dict and args_dict['seq_files']:
        for f in args_dict['seq_files']:
            if not os.path.isfile(f):
                printError('Sequence file %s does not exist.' % f)
            if os.path.splitext(f)[-1].lower() not in seq_types:
                printError('Sequence file %s is not a supported type. Must be one: %s.' \
                           % (f, ', '.join(seq_types)))

    # Verify database files
    if 'db_files' in args_dict and args_dict['db_files']:
        for f in args_dict['db_files']:
            if not os.path.isfile(f):
                printError('Database file %s does not exist.' % f)
            # Strip a trailing .gz extension before checking the underlying file type
            f_check = f[:-3] if f.lower().endswith('.gz') else f
            if os.path.splitext(f_check)[-1].lower() not in db_types:
                printError('Database file %s is not a supported type. Must be one: %s (optionally gzip compressed with a .gz extension).' \
                           % (f, ', '.join(db_types)))

    # Verify non-standard input files
    if in_arg is not None and in_arg in args_dict and args_dict[in_arg]:
        files = args_dict[in_arg] if isinstance(args_dict[in_arg], list) else [args_dict[in_arg]]
        for f in files:
            if not os.path.exists(f):
                printError('Input %s does not exist.' % f)
            if in_types is not None and os.path.splitext(f)[-1].lower() not in in_types:
                printError('Input %s is not a supported type. Must be one: %s.' \
                         % (f, ', '.join(in_types)))

    # Verify output file arguments and exit if anything is hinky
    if args_dict.get('out_files', None) is not None \
            or args_dict.get('out_file', None) is not None:
        if args_dict.get('out_dir', None) is not None:
            printError('The -o argument may not be specified with the --outdir argument.')
        if args_dict.get('out_name', None) is not None:
            printError('The -o argument may not be specified with the --outname argument.')
        if args_dict.get('failed', False):
            printError('The -o argument may not be specified with the --failed argument.')
    if args_dict.get('out_files', None) is not None:
        if len(args_dict['out_files']) != input_count:
            printError('The -o argument requires one output file name per input file.')
        # Resolve the .gz suffix that gzipOutputName will add at write time, so the
        # collision and overwrite checks below see the file name that will actually be opened.
        resolved_out_files = [gzipOutputName(f, args_dict.get('gzip_output', None)) \
                              for f in args_dict['out_files']]
        for f in resolved_out_files:
            if any(_samePath(f, i) for i in input_files):
                printError('Output files and input files cannot have the same names.')
        for f in resolved_out_files:
            if os.path.isfile(f):
                printWarning('Output file %s already exists and will be overwritten.' % f)
    if args_dict.get('out_file', None) is not None:
        resolved_out_file = gzipOutputName(args_dict['out_file'], args_dict.get('gzip_output', None))
        if any(_samePath(resolved_out_file, i) for i in input_files):
            printError('Output files and input files cannot have the same names.')
        if os.path.isfile(resolved_out_file):
            printWarning('Output file %s already exists and will be overwritten.' % resolved_out_file)

    # Exit if output names or log files are specified with multiple input files
    if args_dict.get('out_name', None) is not None \
            and input_count > 1 and not in_list:
        printError('The --outname argument may not be specified with multiple input files.')
    if args_dict.get('log_file', None) is not None \
            and input_count > 1 and not in_list:
        printError('The --log argument may not be specified with multiple input files.')
    
    # Verify output directory
    if 'out_dir' in args_dict and args_dict['out_dir']:
        if os.path.exists(args_dict['out_dir']) and not os.path.isdir(args_dict['out_dir']):
            printError('Path %s exists but it is not a directory.' % args_dict['out_dir'])

    # Redefine common output options as out_args dictionary
    out_args = ['log_file', 'out_dir', 'out_name', 'out_type', 'failed', 'gzip_output']
    args_dict['out_args'] = {k:args_dict.setdefault(k, None) for k in out_args}
    for k in out_args: del args_dict[k]
    
    return args_dict


def checkArgs(parser):
    """
    Checks that arguments have been provided and prints help if they have not.

    Arguments:
      parser : An argparse.ArgumentParser defining the commandline arguments.

    Returns:
      boolean : True if arguments are present. Prints help and exits if not.
    """
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    return True

def setDefaultFields(args, defaults, format='airr'):
    """
    Sets default field arguments by format

    Arguments:
      args (dict): parsed argument dictionary.
      defaults (dict): default variables to set with with keys as argument variables and values
                       as AIRR field names.
      format (str): one of 'changeo' or 'airr' which defines the file format.

    Returns:
      dict: modified input args.
    """
    if format == 'changeo':
        defaults = {k: ChangeoSchema.fromReceptor(AIRRSchema.toReceptor(v)) \
                    for k, v in defaults.items()}
    for f in defaults:
        if args[f] is None:  args[f] = defaults[f]

    return(args)
