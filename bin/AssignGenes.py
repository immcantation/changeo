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
from textwrap import dedent
from time import time
from Bio import SeqIO

# Presto imports
from presto.IO import printLog, printMessage, printProgress
from changeo.Applications import runIgBLAST
from changeo.Commandline import CommonHelpFormatter, checkArgs, getCommonArgParser, parseCommonArgs
from changeo.Defaults import default_igblast_exec, default_out_args

# Defaults
default_igdata = '~/share/igblast'

def assignIgBLAST(seq_file, igdata=default_igdata, locus='ig', organism='human',
                  exec=default_igblast_exec,
                  out_file=None, out_args=default_out_args, nproc=None):
    """
    Performs clustering on sets of sequences

    Arguments:
      seq_file : the sample sequence file name.
      igdata (str): path to the IgBLAST database directory (IGDATA environment).
      locus (str): receptor type; one of 'ig' or 'tr'.
      organism (str): species name.
      exec : the path to the igblastn executable.
      out_file : output file name. Automatically generated from the input file if None.
      out_args : common output argument dictionary from parseCommonArgs.
      nproc : the number of processQueue processes;
              if None defaults to the number of CPUs.

    Returns:
      str: the output file name
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'AssignGenes'
    log['COMMAND'] = 'igblast'
    log['FILE'] = os.path.basename(seq_file)
    log['NPROC'] = nproc
    printLog(log)

    out_file = ''
    igdata = ''

    # Run IgBLAST clustering
    start_time = time()
    printMessage('Running IgBLAST', start_time=start_time, width=25)
    console_out = runIgBLAST(seq_file, igdata, locus=locus, organism=organism, output=out_file,
                             threads=nproc, exec=exec)
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
                 fmt7
                    IgBLAST output.
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
    parent_parser = getCommonArgParser(multiproc=True)

    # Subparser to run IgBLAT
    parser_igblast = subparsers.add_parser('igblast', parents=[parent_parser],
                                           formatter_class=CommonHelpFormatter, add_help=False,
                                           help='Executes IgBLAST.',
                                           description='Executes IgBLAST.')
    group_igblast = parser_igblast.add_argument_group('alignment arguments')
    group_igblast.add_argument('-r', nargs='+', action='store', dest='repo', required=False,
                               help='''List of folders and/or fasta files containing
                                    IMGT-gapped germline sequences corresponding to the
                                    set of germlines used for the alignment. Specifying 
                                    this argument is not required, but doing so will add IMGT-gaps 
                                    to SEQUENCE_IMGT (sequence_alignment) and rebuild 
                                    GERMLINE_IMGT (germline_alignment). Requires the 
                                    V_GERM_START_IMGT (v_germline_start) and 
                                    V_GERM_LENGTH_IMGT (v_germline_end) fields.''')
    group_igblast.add_argument('--format', action='store', dest='format', default=default_format,
                               choices=choices_format,
                               help='''Specify the input format. Output will always be AIRR TSV, but 
                                    specifying AIRR input will allow IMGT-gapping of the alignment 
                                    columns.''')
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

    # Call cluster main function for each input file
    del args_dict['seq_files']
    del args_dict['func']
    del args_dict['command']
    for f in args.__dict__['seq_files']:
        args_dict['seq_file'] = f
        args.func(**args_dict)
