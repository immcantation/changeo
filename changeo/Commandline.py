"""
Commandline interface functions
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
# TODO:  Can remove if we remove the `annotation` argument from `getCommonArgParser`
default_delimiter = ('|', '=', ',')


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


def getCommonArgParser(seq_in=True, seq_out=True, paired=False, db_in=False, db_out=False,
                       failed=True, log=True, annotation=True, multiproc=False):
    """
    Defines an ArgumentParser object with common pRESTO arguments

    Arguments:
      seq_in : If True include sequence input arguments
      seq_out : If True include sequence output arguments
      paired : If True defined paired-end sequence input and output arguments
      db_in : If True include tab delimited database input arguments
      db_out : If True include tab delimited database output arguments
      failed : If True include arguments for output of failed results
      log : If True include log arguments
      annotation : If True include annotation arguments
      multiproc : If True include multiprocessing arguments

    Returns:
      ArgumentParser : An ArgumentParser object
    """
    parser = ArgumentParser(add_help=False, formatter_class=CommonHelpFormatter)

    # Database arguments
    if db_in:
        parser.add_argument('-d', nargs='+', action='store', dest='db_files', required=True,
                        help='A list of tab delimited database files.')
    if db_out:
        # Place holder for the future
        pass

    # Sequence arguments
    if seq_in and not paired:
        parser.add_argument('-s', nargs='+', action='store', dest='seq_files', required=True,
                            help='A list of FASTA/FASTQ files containing sequences to process.')
    elif seq_in and paired:
        parser.add_argument('-1', nargs='+', action='store', dest='seq_files_1', required=True,
                            help='''An ordered list of FASTA/FASTQ files containing
                                 head/primary sequences.''')
        parser.add_argument('-2', nargs='+', action='store', dest='seq_files_2', required=True,
                            help='''An ordered list of FASTA/FASTQ files containing
                                 tail/secondary sequences.''')
    if seq_out:
        parser.add_argument('--fasta', action='store_const', dest='out_type', const='fasta',
                            help='Specify to force output as FASTA rather than FASTQ.')

    # Failed result arguments
    if failed:
        parser.add_argument('--failed', action='store_true', dest='failed',
                            help='''If specified create files containing records that
                                  fail processing.''')

    # Log arguments
    if log:
        parser.add_argument('--log', action='store', dest='log_file', default=None,
                            help='''Specify to write verbose logging to a file. May not be
                                  specified with multiple input files.''')

    # Annotation arguments
    if annotation:
        parser.add_argument('--delim', nargs=3, action='store', dest='delimiter',
                            type=str, default=default_delimiter,
                            help='''A list of the three delimiters that separate annotation
                                 blocks, field names and values, and values within a field,
                                 respectively.''')

    # Multiprocessing arguments
    if multiproc:
        parser.add_argument('--nproc', action='store', dest='nproc', type=int, default=mp.cpu_count(),
                            help='The number of simultaneous computational processes to execute \
                                  (CPU cores to utilized).')

    # Universal arguments
    parser.add_argument('--outdir', action='store', dest='out_dir', default=None,
                        help='Specify to changes the output directory to the location specified. \
                              The input file directory is used if this is not specified.')
    parser.add_argument('--outname', action='store', dest='out_name', default=None,
                        help='Changes the prefix of the successfully processed output file \
                              to the string specified. May not be specified with multiple \
                              input files.')

    return parser


def parseCommonArgs(args, in_arg=None, in_types=None):
    """
    Checks common arguments from getCommonArgParser and transforms output options to a dictionary

    Arguments: 
      args : Argument Namespace defined by ArgumentParser.parse_args
      in_arg : String defining a non-standard input file argument to verify;
               by default ['db_files', 'seq_files', 'seq_files_1', 'seq_files_2', 'primer_file']
               are supported in that order
      in_types : List of types (file extensions as strings) to allow for files in file_arg
                 if None do not check type
                    
    Returns:
      dict : Dictionary copy of args with output arguments embedded in the dictionary out_args
    """ 
    db_types = ['.tab', '.tsv']
    seq_types = ['.fasta', 'fna', '.fa', '.fastq', '.fq']
    primer_types = ['.fasta', '.regex']
    if in_types is not None:  in_types = [f.lower for f in in_types]
    args_dict = args.__dict__.copy()
    
    # Count input files
    if 'db_files' in args_dict:      
        input_count = len(args_dict['db_files'] or [])
    elif 'seq_files' in args_dict:       
        input_count = len(args_dict['seq_files']  or [])
    elif all([k in args_dict for k in ('seq_files_1', 'seq_files_2')]):
        input_count = len(args_dict['seq_files_1']  or [])
    elif 'primer_file' in args_dict:   
        input_count = 1
    elif in_arg is not None and in_arg in args_dict: 
        input_count = len(args_dict[in_arg] or [])
    else:
        sys.exit('ERROR:  Cannot determine input file argument')
    
    # Exit if output names or log files are specified with multiple input files    
    if args_dict.get('out_name', None) is not None and input_count > 1:
        sys.exit('ERROR:  The --outname argument may not be specified with multiple input files')
    if args_dict.get('log_file', None) is not None and input_count > 1:
        sys.exit('ERROR:  The --log argument may not be specified with multiple input files')
        
    # Verify database files
    if 'db_files' in args_dict and args_dict['db_files']:
        for f in args_dict['db_files']:
            if not os.path.isfile(f):
                sys.exit('ERROR:  Database file %s does not exist' % f)
            if os.path.splitext(f)[-1].lower() not in db_types:
                sys.exit('ERROR:  Database file %s is not a supported type. Must be one: %s' \
                         % (f, ', '.join(db_types)))
    
    # Verify single-end sequence files
    if 'seq_files' in args_dict and args_dict['seq_files']:
        for f in args_dict['seq_files']:
            if not os.path.isfile(f):
                sys.exit('ERROR:  Sequence file %s does not exist' % f)
            if os.path.splitext(f)[-1].lower() not in seq_types:
                sys.exit('ERROR:  Sequence file %s is not a supported type. Must be one: %s' \
                         % (f, ', '.join(seq_types)))
    
    # Verify paired-end sequence files
    if all([k in args_dict and args_dict[k] for k in ('seq_files_1', 'seq_files_2')]):
        if len(args_dict['seq_files_1']) != len(args_dict['seq_files_2']):
            sys.exit('ERROR:  The -1 and -2 arguments must contain the same number of files')
        for f1, f2 in zip(args_dict['seq_files_1'], args_dict['seq_files_2']):
            if os.path.splitext(f1)[-1].lower() != os.path.splitext(f2)[-1].lower():
                sys.exit('ERROR:  Each pair of files in the -1 and -2 arguments must be the same file type')
        for f in (args_dict['seq_files_1'] + args_dict['seq_files_2']):
            if not os.path.isfile(f):
                sys.exit('ERROR:  Sequence file %s does not exist' % f)
            if os.path.splitext(f)[-1].lower() not in seq_types:
                sys.exit('ERROR:  Sequence file %s is not a supported type. Must be one: %s' \
                         % (f, ', '.join(seq_types)))
                    
    # Verify primer file
    if 'primer_file' in args_dict and args_dict['primer_file']:
        primer_file = args_dict['primer_file']
        if not os.path.isfile(primer_file):
            sys.exit('ERROR:  Primer file %s does not exist' % primer_file)
        if os.path.splitext(primer_file)[-1].lower() not in primer_types:
            sys.exit('ERROR:  Primer file %s is not a supported type. Must be one: %s' \
                     % (primer_file, ', '.join(primer_types)))
    
    # Verify non-standard input files
    if in_arg is not None and in_arg in args_dict and args_dict[in_arg]:
        files = args_dict[in_arg] if isinstance(args_dict[in_arg], list) \
                else [args_dict[in_arg]]
        for f in files:
            if not os.path.exists(f):
                sys.exit('ERROR:  Input %s does not exist' % f)
            if in_types is not None and os.path.splitext(f)[-1].lower() not in in_types:
                sys.exit('ERROR:  Input %s is not a supported type. Must be one: %s' \
                         % (f, ', '.join(in_types)))
    
    # Verify output directory
    if 'out_dir' in args_dict and args_dict['out_dir']:
        if os.path.exists(args_dict['out_dir']) and not os.path.isdir(args_dict['out_dir']):
            sys.exit('ERROR:  Directory %s exists but it is not a directory' % args_dict['out_dir'])

    # Redefine common output options as out_args dictionary
    out_args = ['log_file', 'delimiter', 'separator', 
                'out_dir', 'out_name', 'out_type', 'failed']
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