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
from argparse import ArgumentParser
from collections import OrderedDict
from pkg_resources import parse_version
from textwrap import dedent
from time import time
import re
import Bio
from Bio import SeqIO

# Presto imports
from presto.IO import printLog, printMessage, printError, printWarning
from changeo.Defaults import default_igblastn_exec, default_igblastp_exec, default_out_args
from changeo.Applications import runIgBLASTN, runIgBLASTP, getIgBLASTVersion
from changeo.Commandline import CommonHelpFormatter, checkArgs, getCommonArgParser, parseCommonArgs
from changeo.IO import getOutputName

# Defaults
choices_format = ('blast', 'airr')
choices_loci = ('ig', 'tr')
choices_organism = ('human', 'mouse', 'rabbit', 'rat', 'rhesus_monkey')
default_format = 'blast'
default_loci = 'ig'
default_organism = 'human'
default_igdata = '~/share/igblast'


def assignIgBLAST(seq_file, amino_acid=False, igdata=default_igdata, loci='ig', organism='human',
                  vdb=None, ddb=None, jdb=None, cdb=None, format=default_format,
                  igblast_exec=default_igblastn_exec, out_file=None,
                  out_args=default_out_args, nproc=None):
    """
    Performs clustering on sets of sequences

    Arguments:
      seq_file (str): the sample sequence file name.
      amino_acid : if True then run igblastp. igblastn is assumed if False.
      igdata (str): path to the IgBLAST database directory (IGDATA environment).
      loci (str): receptor type; one of 'ig' or 'tr'.
      organism (str): species name.
      vdb (str): name of a custom V reference in the database folder to use.
      ddb (str): name of a custom D reference in the database folder to use.
      jdb (str): name of a custom J reference in the database folder to use.
      cdb (str): name of a custom C reference in the database folder to use.
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
        printError('Invalid output format %s.' % format)

    # Get IgBLAST version
    version = getIgBLASTVersion(exec=igblast_exec)
    if parse_version(version) < parse_version('1.6'):
        printError('IgBLAST version is %s and 1.6 or higher is required.' % version)
    if format == 'airr' and parse_version(version) < parse_version('1.9'):
        printError('IgBLAST version is %s and 1.9 or higher is required for AIRR format support.' % version)

    # Print parameter info
    log = OrderedDict()
    log['START'] = 'AssignGenes'
    log['COMMAND'] = 'igblast-aa' if amino_acid else 'igblast'
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
    # convert to FASTA if needed
    infile = open(seq_file, 'r')
    test = infile.read()[0]
    if test == "@":
        printMessage("Running conversion from FASTQ to FASTA")
        fasta_out_dir, filename = os.path.split(out_file)
        out_fasta_file = os.path.split(seq_file)[1]
        out_fasta_file = os.path.join(fasta_out_dir,'%s.fasta' % os.path.splitext(out_fasta_file)[0])
        with open(out_fasta_file, "w") as out_handle:
            records = SeqIO.parse(seq_file, 'fastq')
            if parse_version(Bio.__version__) >= parse_version('1.71'):
                # Biopython >= v1.71
                SeqIO.write(records, out_handle, format='fasta-2line')
            else:
                # Biopython < v1.71
                writer = SeqIO.FastaIO.FastaWriter(out_handle, wrap=None)
                writer.write_file(records)
        seq_file = out_fasta_file

    # Run IgBLAST clustering
    start_time = time()
    printMessage('Running IgBLAST', start_time=start_time, width=25)
    if not amino_acid:
        console_out = runIgBLASTN(seq_file, igdata, loci=loci, organism=organism,
                                  vdb=vdb, ddb=ddb, jdb=jdb, cdb=cdb, output=out_file,
                                  format=format, threads=nproc, exec=igblast_exec)
    else:
        console_out = runIgBLASTP(seq_file, igdata, loci=loci, organism=organism,
                                  vdb=vdb, output=out_file,
                                  threads=nproc, exec=igblast_exec)
    printMessage('Done', start_time=start_time, end=True, width=25)

    # Get number of processed sequences
    if (format == 'blast'):
        with open(out_file, 'rb') as f:
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b'\n':
                f.seek(-2, os.SEEK_CUR)
            pass_info = f.readline().decode()
        num_seqs_match = re.search('(# BLAST processed )(\d+)( .*)', pass_info)
        num_sequences = num_seqs_match.group(2)
    else:
        f = open(out_file, 'rb')
        lines = 0
        buf_size = 1024 * 1024
        read_f = f.raw.read
        buf = read_f(buf_size)
        while buf:
            lines += buf.count(b'\n')
            buf = read_f(buf_size)
        num_sequences = lines - 1

    # Print log
    log = OrderedDict()
    log['PASS'] = num_sequences
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
                            version='%(prog)s:' + ' %s %s' %(__version__, __date__))
    group_help.add_argument('-h', '--help', action='help', help='show this help message and exit')
    subparsers = parser.add_subparsers(title='subcommands', dest='command', metavar='',
                                       help='Assignment operation')
    # TODO:  This is a temporary fix for Python issue 9253
    subparsers.required = True

    # Parent parser
    parent_parser = getCommonArgParser(db_in=False, log=False, failed=False, format=False, multiproc=True)

    # Subparser to run igblastn
    parser_igblast = subparsers.add_parser('igblast', parents=[parent_parser],
                                           formatter_class=CommonHelpFormatter, add_help=False,
                                           help='Executes igblastn.',
                                           description='Executes igblastn.')
    group_igblast = parser_igblast.add_argument_group('alignment arguments')
    group_igblast.add_argument('-s', nargs='+', action='store', dest='seq_files', required=True,
                               help='A list of FASTA files containing sequences to process.')
    group_igblast.add_argument('-b', action='store', dest='igdata', required=True,
                               help='IgBLAST database directory (IGDATA).')
    group_igblast.add_argument('--organism', action='store', dest='organism', default=default_organism,
                               choices=choices_organism, help='Organism name.')
    group_igblast.add_argument('--loci', action='store', dest='loci', default=default_loci,
                               choices=choices_loci, help='The receptor type.')
    group_igblast.add_argument('--vdb', action='store', dest='vdb', default=None,
                               help='''Name of the custom V reference in the IgBLAST database folder.
                                    If not specified, then a default database name with the form
                                    imgt_<organism>_<loci>_v will be used.''')
    group_igblast.add_argument('--ddb', action='store', dest='ddb', default=None,
                               help='''Name of the custom D reference in the IgBLAST database folder.
                                    If not specified, then a default database name with the form
                                    imgt_<organism>_<loci>_d will be used.''')
    group_igblast.add_argument('--jdb', action='store', dest='jdb', default=None,
                               help='''Name of the custom J reference in the IgBLAST database folder.
                                    If not specified, then a default database name with the form
                                    imgt_<organism>_<loci>_j will be used.''')
    group_igblast.add_argument('--cdb', action='store', dest='cdb', default=None,
                               help='''Name of the custom C reference in the IgBLAST database folder.
                                    If not specified, then a default database name with the form 
                                    imgt_<organism>_<loci>_c will be used. Note, this argument will be 
                                    ignored for IgBLAST versions below 1.18.0.''')
    group_igblast.add_argument('--format', action='store', dest='format', default=default_format,
                               choices=choices_format,
                               help='''Specify the output format. The "blast" will result in
                                    the IgBLAST "-outfmt 7 std qseq sseq btop" output format.
                                    Specifying "airr" will output the AIRR TSV format provided by
                                    the IgBLAST argument "-outfmt 19".''')
    group_igblast.add_argument('--exec', action='store', dest='igblast_exec',
                               default=default_igblastn_exec,
                               help='Path to the igblastn executable.')
    parser_igblast.set_defaults(func=assignIgBLAST, amino_acid=False)

    # Subparser to run igblastp
    parser_igblast_aa = subparsers.add_parser('igblast-aa', parents=[parent_parser],
                                           formatter_class=CommonHelpFormatter, add_help=False,
                                           help='Executes igblastp.',
                                           description='Executes igblastp.')
    group_igblast_aa = parser_igblast_aa.add_argument_group('alignment arguments')
    group_igblast_aa.add_argument('-s', nargs='+', action='store', dest='seq_files', required=True,
                                  help='A list of FASTA files containing sequences to process.')
    group_igblast_aa.add_argument('-b', action='store', dest='igdata', required=True,
                                  help='IgBLAST database directory (IGDATA).')
    group_igblast_aa.add_argument('--organism', action='store', dest='organism', default=default_organism,
                                  choices=choices_organism, help='Organism name.')
    group_igblast_aa.add_argument('--loci', action='store', dest='loci', default=default_loci,
                                  choices=choices_loci, help='The receptor type.')
    group_igblast_aa.add_argument('--vdb', action='store', dest='vdb', default=None,
                                  help='''Name of the custom V reference in the IgBLAST database folder.
                                       If not specified, then a default database name with the form
                                       imgt_aa_<organism>_<loci>_v will be used.''')
    group_igblast_aa.add_argument('--exec', action='store', dest='igblast_exec',
                                  default=default_igblastp_exec,
                                  help='Path to the igblastp executable.')
    parser_igblast_aa.set_defaults(func=assignIgBLAST, amino_acid=True, ddb=None, jdb=None, format='blast')


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
