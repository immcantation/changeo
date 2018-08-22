#!/usr/bin/env python3
"""
Assign V(D)J gene annotations
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from changeo import __version__, __date__

# Imports
import os
import shutil
import sys
from argparse import ArgumentParser
from collections import OrderedDict
from pkg_resources import parse_version
from textwrap import dedent
from time import time

# Presto imports
from presto.IO import printLog, printMessage
from changeo.Applications import runIgBLAST, getIgBLASTVersion
from changeo.Commandline import CommonHelpFormatter, checkArgs, getCommonArgParser, parseCommonArgs
from changeo.Defaults import default_igblast_exec, default_out_args

# Defaults
choices_format = ('blast', 'airr')
choices_organism = ('human', 'mouse')
choices_loci = ('ig', 'tr')

default_igdata = '~/share/igblast'
default_format = 'blast'
default_organism = 'human'
default_loci = 'ig'


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
    # Get in_file components
    dir_name, file_name = os.path.split(file)
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

    # Return file name
    return out_file


def assignIgBLAST(seq_file, igdata=default_igdata, loci='ig', organism='human', format=default_format,
                  igblast_exec=default_igblast_exec, out_file=None, out_args=default_out_args,
                  nproc=None):
    """
    Performs clustering on sets of sequences

    Arguments:
      seq_file (str): the sample sequence file name.
      igdata (str): path to the IgBLAST database directory (IGDATA environment).
      loci (str): receptor type; one of 'ig' or 'tr'.
      organism (str): species name.
      format (str): output format. One of 'blast' or 'airr'.
      exec (str): the path to the igblastn executable.
      out_file (str): output file name. Automatically generated from the input file if None.
      out_args (dict): common output argument dictionary from parseCommonArgs.
      nproc (int): the number of processQueue processes;
              if None defaults to the number of CPUs.

    Returns:
      str: the output file name
    """
    # Check format argument
    try:
        out_type = {'blast': 'fmt7', 'airr': 'tsv'}[format]
    except KeyError:
        sys.exit('Error: Invalid output format %s' % format)

    # Get IgBLAST version
    version = getIgBLASTVersion(exec=igblast_exec)
    if parse_version(version) < parse_version('1.6'):
        sys.exit('Error: IgBLAST version is %s and 1.6 or higher is required.' % version)
    if format == 'airr' and parse_version(version) < parse_version('1.9'):
        sys.exit('Error: IgBLAST version is %s and 1.9 or higher is required for AIRR format support.' % version)

    # Print parameter info
    log = OrderedDict()
    log['START'] = 'AssignGenes'
    log['COMMAND'] = 'igblast'
    log['VERSION'] = version
    log['FILE'] = os.path.basename(seq_file)
    log['ORGANISM'] = organism
    log['LOCI'] = loci
    log['NPROC'] = nproc
    printLog(log)

    # Open output writer
    if out_file is None:
        out_file = getOutputName(seq_file, out_label='igblast', out_dir=out_args['out_dir'],
                                 out_name=out_args['out_name'], out_type=out_type)

    # Run IgBLAST clustering
    start_time = time()
    printMessage('Running IgBLAST', start_time=start_time, width=25)
    console_out = runIgBLAST(seq_file, igdata, loci=loci, organism=organism, output=out_file,
                             format=format, threads=nproc, exec=igblast_exec)
    printMessage('Done', start_time=start_time, end=True, width=25)

    # Print log
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(out_file)
    log['END'] = 'AssignGenes'
    printLog(log)

    return out_file


def getArgParser():
    """
    Defines the ArgumentParser

    Arguments:
    None
                      
    Returns: 
    an ArgumentParser object
    """
    # Define output file names and header fields
    fields = dedent(
             '''
             output files:
                 igblast
                    Reference alignment results from IgBLAST.
             ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            formatter_class=CommonHelpFormatter, add_help=False)
    group_help = parser.add_argument_group('help')
    group_help.add_argument('--version', action='version',
                            version='%(prog)s:' + ' %s-%s' %(__version__, __date__))
    group_help.add_argument('-h', '--help', action='help', help='show this help message and exit')
    subparsers = parser.add_subparsers(title='subcommands', dest='command', metavar='',
                                       help='Assignment operation')
    # TODO:  This is a temporary fix for Python issue 9253
    subparsers.required = True

    # Parent parser
    parent_parser = getCommonArgParser(db_in=False, format=False, multiproc=True)

    # Subparser to run IgBLAT
    parser_igblast = subparsers.add_parser('igblast', parents=[parent_parser],
                                           formatter_class=CommonHelpFormatter, add_help=False,
                                           help='Executes IgBLAST.',
                                           description='Executes IgBLAST.')
    group_igblast = parser_igblast.add_argument_group('alignment arguments')
    group_igblast.add_argument('-s', nargs='+', action='store', dest='seq_files', required=True,
                               help='A list of FASTA files containing sequences to process.')
    group_igblast.add_argument('-b', action='store', dest='igdata', required=True,
                               help='IgBLAST database directory (IGDATA).')
    group_igblast.add_argument('--organism', action='store', dest='organism', default=default_organism,
                               choices=choices_organism, help='Organism name.')
    group_igblast.add_argument('--loci', action='store', dest='loci', default=default_loci,
                               choices=choices_loci, help='The receptor type.')
    group_igblast.add_argument('--format', action='store', dest='format', default=default_format,
                               choices=choices_format,
                               help='''Specify the output format. The "blast" will result in
                                    the IgBLAST "-outfmt 7 std qseq sseq btop" output format.
                                    Specifying "airr" will output the AIRR TSV format provided by
                                    the IgBLAST argument "-outfmt 19".''')
    group_igblast.add_argument('--exec', action='store', dest='igblast_exec',
                              default=default_igblast_exec,
                              help='Path to the igblastn executable.')
    parser_igblast.set_defaults(func=assignIgBLAST)

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

    # Check if a valid clustering executable was specified
    if not shutil.which(args_dict['igblast_exec']):
        parser.error('%s executable not found' % args_dict['igblast_exec'])

    # Clean arguments dictionary
    del args_dict['seq_files']
    if 'out_files' in args_dict: del args_dict['out_files']
    del args_dict['func']
    del args_dict['command']

    # Call main function for each input file
    for i, f in enumerate(args.__dict__['seq_files']):
        args_dict['seq_file'] = f
        args_dict['out_file'] = args.__dict__['out_files'][i] \
            if args.__dict__['out_files'] else None
        args.func(**args_dict)
