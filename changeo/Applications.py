"""
Application wrappers
"""

# Info
__author__ = 'Jason Anthony Vander Heiden'

# Imports
import os
import re
from subprocess import check_output, STDOUT, CalledProcessError

# Presto and changeo imports
from presto.IO import printError, printWarning
from changeo.Defaults import default_igblast_exec, default_tbl2asn_exec, default_igphyml_exec

# Defaults
default_igblast_output = 'legacy'


def runASN(fasta, template=None, exec=default_tbl2asn_exec):
    """
    Executes tbl2asn to generate Sequin files

    Arguments:
      fasta (str): fsa file name.
      template (str): sbt file name.
      exec (str): the name or path to the tbl2asn executable.

    Returns:
      str: tbl2asn console output.
    """
    # Basic command that requires .fsa and .tbl files in the same directory
    # tbl2asn -i records.fsa -a s -V vb -t template.sbt

    # Define tbl2asn command
    cmd = [exec,
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
        printError('Running command: %s\n%s' % (' '.join(cmd), e.output))

    if 'Unable to read any FASTA records' in stdout_str:
        printError('%s failed: %s' % (' '.join(cmd), stdout_str))

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
        stdout_str = check_output(cmd, stderr=STDOUT, shell=False,
                                  universal_newlines=True, cwd=rep_dir)
    except CalledProcessError as e:
        printError('Running command: %s\n%s' % (' '.join(cmd), e.output))

    return None


def runIgBLAST(fasta, igdata, loci='ig', organism='human', vdb=None, ddb=None, jdb=None, output=None,
               format=default_igblast_output, threads=1, exec=default_igblast_exec):
    """
    Runs IgBLAST on a sequence file

    Arguments:
      fasta (str): fasta file containing sequences.
      igdata (str): path to the IgBLAST database directory (IGDATA environment).
      loci (str): receptor type; one of 'ig' or 'tr'.
      organism (str): species name.
      vdb (str): name of a custom V reference in the database folder to use.
      ddb (str): name of a custom D reference in the database folder to use.
      jdb (str): name of a custom J reference in the database folder to use.
      output (str): output file name. If None, automatically generate from the fasta file name.
      format (str): output format. One of 'blast' or 'airr'.
      threads (int): number of threads for igblastn.
      exec (str): the name or path to the igblastn executable.

    Returns:
      str: IgBLAST console output.

    """
    # export IGDATA
    # declare -A SEQTYPE
    # SEQTYPE[ig] = "Ig"
    # SEQTYPE[tr] = "TCR"
    # GERMLINE_V = "imgt_${SPECIES}_${RECEPTOR}_v"
    # GERMLINE_D = "imgt_${SPECIES}_${RECEPTOR}_d"
    # GERMLINE_J = "imgt_${SPECIES}_${RECEPTOR}_j"
    # AUXILIARY = "${SPECIES}_gl.aux"
    # IGBLAST_DB = "${IGDATA}/database"
    # IGBLAST_CMD = "igblastn \
    #     -germline_db_V ${IGBLAST_DB}/${GERMLINE_V} \
    #     -germline_db_D ${IGBLAST_DB}/${GERMLINE_D} \
    #     -germline_db_J ${IGBLAST_DB}/${GERMLINE_J} \
    #     -auxiliary_data ${IGDATA}/optional_file/${AUXILIARY} \
    #     -ig_seqtype ${SEQTYPE[${RECEPTOR}]} -organism ${SPECIES} \
    #     -domain_system imgt -outfmt '7 std qseq sseq btop'"
    #
    # # Set run commmand
    # OUTFILE =$(basename ${READFILE})
    # OUTFILE = "${OUTDIR}/${OUTFILE%.fasta}.fmt7"
    # IGBLAST_VER =$(${IGBLAST_CMD} -version | grep 'Package' | sed s / 'Package: ' //)
    # IGBLAST_RUN = "${IGBLAST_CMD} -query ${READFILE} -out ${OUTFILE} -num_threads ${NPROC}"

    try:
        outfmt = {'blast': '7 std qseq sseq btop', 'airr': '19'}[format]
    except KeyError:
        printError('Invalid output format %s.' % format)

    try:
        seqtype = {'ig': 'Ig', 'tr': 'TCR'}[loci]
    except KeyError:
        printError('Invalid receptor type %s.' % loci)

    # Set auxilary data
    auxilary = os.path.join(igdata, 'optional_file', '%s_gl.aux' % organism)
    # Set V database
    if vdb is not None:  v_germ = os.path.join(igdata, 'database', vdb)
    else:  v_germ = os.path.join(igdata, 'database', 'imgt_%s_%s_v' % (organism, loci))
    # Set D database
    if ddb is not None:  d_germ = os.path.join(igdata, 'database', ddb)
    else:  d_germ = os.path.join(igdata, 'database', 'imgt_%s_%s_d' % (organism, loci))
    # Set J database
    if jdb is not None:  j_germ = os.path.join(igdata, 'database', jdb)
    else:  j_germ = os.path.join(igdata, 'database', 'imgt_%s_%s_j' % (organism, loci))

    # Define IgBLAST command
    cmd = [exec,
           '-query', os.path.abspath(fasta),
           '-out', os.path.abspath(output),
           '-num_threads', str(threads),
           '-ig_seqtype', seqtype,
           '-organism', organism,
           '-auxiliary_data', str(auxilary),
           '-germline_db_V', str(v_germ),
           '-germline_db_D', str(d_germ),
           '-germline_db_J', str(j_germ),
           '-outfmt', outfmt,
           '-domain_system', 'imgt']

    # Execute IgBLAST
    env = os.environ.copy()
    env['IGDATA'] = igdata
    try:
        stdout_str = check_output(cmd, stderr=STDOUT, shell=False, env=env,
                                  universal_newlines=True)
    except CalledProcessError as e:
        printError('Running command: %s\n%s' % (' '.join(cmd), e.output))

    #if 'Unable to read any FASTA records' in stdout_str:
    #    sys.stderr.write('\n%s failed: %s\n' % (' '.join(cmd), stdout_str))

    return stdout_str


def getIgBLASTVersion(exec=default_igblast_exec):
    """
    Gets the version of the IgBLAST executable

    Arguments:
      exec (str): the name or path to the igblastn executable.

    Returns:
      str: version number.
    """
    # Build commandline
    cmd = [exec, '-version']

    # Run
    try:
        stdout_str = check_output(cmd, stderr=STDOUT, shell=False, universal_newlines=True)
    except CalledProcessError as e:
        printError('Running command: %s\n%s' % (' '.join(cmd), e.output))

    # Extract version number
    match = re.search('(?<=Package: igblast )(\d+\.\d+\.\d+)', stdout_str)
    version = match.group(0)

    return version