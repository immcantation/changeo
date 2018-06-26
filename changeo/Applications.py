"""
Application wrappers
"""

# Info
__author__ = 'Jason Anthony Vander Heiden'

# Imports
import os
import sys
import subprocess
from subprocess import check_output, STDOUT, CalledProcessError

# Defaults
default_tbl2asn_exec='tbl2asn'
default_igphyml_exec = 'igphyml'


def runASN(fasta, template=None, tbl2asn_exec=default_tbl2asn_exec):
    """
    Executes tbl2asn to generate Sequin files

    Arguments:
      fasta : fsa file.
      template : sbt file.
      tbl2asn_exec : the name or path to the tbl2asn executable.

    Returns:
      str : tbl2asn console output.
    """
    # Basic command that requires .fsa and .tbl files in the same directory
    # tbl2asn -i records.fsa -a s -V vb -t template.sbt

    # Define tb2asn command
    cmd = [tbl2asn_exec,
           '-i', os.path.abspath(fasta),
           '-a', 's',
           '-V', 'vb']
    if template is not None:
        cmd.extend(['-t', os.path.abspath(template)])

    # Execute tbl2asn
    try:
        stdout_str = check_output(cmd, stderr=STDOUT, shell=False,
                                  universal_newlines=True)
    except CalledProcessError as e:
        sys.stderr.write('\nError running command: %s\n' % ' '.join(cmd))
        sys.exit(e.output)

    if 'Unable to read any FASTA records' in stdout_str:
        sys.stderr.write('\n%s failed: %s\n' % (' '.join(cmd), stdout_str))

    return stdout_str


def runIgPhyML(rep_file, rep_dir, model='HLP17', motifs='FCH',
               threads=1, exec=default_igphyml_exec):
    """
    Run IgPhyML

    Arguments:
      rep_file (str): repertoire tsv file.
      rep_dir (str): directory containing input fasta files.
      model (str): model to use.
      motif (str): motifs argument.
      threads : number of threads.
      exec : the path to the IgPhyMl executable.

    Returns:
      str: name of the output tree file.
    """
    # cd rep_dir
    # igphyml --repfile rep_file -m HLP17 --motifs FCH --omegaOpt e,e --run_id test -o tlr --threads 4 --minSeq 2

    # Define igphyml command
    cmd = [exec,
           '--repfile', rep_file,
           '-m', model,
           '--motifs', motifs,
           '--omegaOpt',  'e,e',
           '-o', 'tlr',
           '--minSeq', '2',
           '--threads', str(threads)]

    # Run IgPhyMl
    try:
        stdout_str = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=False,
                                             universal_newlines=True, cwd=rep_dir)
    except subprocess.CalledProcessError as e:
        sys.stderr.write('\nError running command: %s\n' % ' '.join(cmd))
        sys.exit(e.output)

    return None