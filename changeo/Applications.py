"""
Application wrappers
"""

# Info
__author__ = 'Jason Anthony Vander Heiden'

# Imports
import glob
import os
import re
from packaging.version import parse as parse_version
from subprocess import check_output, STDOUT, CalledProcessError

# Presto and changeo imports
from presto.IO import printError, printWarning
from changeo.Defaults import default_igblastn_exec, default_igblastp_exec, default_tbl2asn_exec, \
                             default_igphyml_exec

# Defaults
default_igblast_output = 'legacy'


def _resolve_db(igdata, segment, organism, loci, explicit=None, aa=False):
    """
    Resolve an IgBLAST database path by probing the filesystem.

    Tries the 'airrc-imgt' prefix first, then 'imgt'. Returns None if neither is
    found and no explicit name was given.

    Arguments:
      igdata (str): path to the IgBLAST database directory.
      segment (str): gene segment; one of 'v', 'd', 'j', or 'c'.
      organism (str): species name.
      loci (str): receptor type; one of 'ig' or 'tr'.
      explicit (str): explicit database name, bypassing prefix resolution.
      aa (bool): if True, use the amino acid infix ('airrc-imgt_aa_' / 'imgt_aa_').

    Returns:
      str: full path to the database, or None if not found.
    """
    if explicit is not None:
        return os.path.join(igdata, 'database', explicit)

    infix = '_aa_' if aa else '_'
    db_dir = os.path.join(igdata, 'database')
    for prefix in ('airrc-imgt', 'imgt'):
        name = '%s%s%s_%s_%s' % (prefix, infix, organism, loci, segment)
        path = os.path.join(db_dir, name)
        if glob.glob(path + '.*') or os.path.isfile(path):
            return path

    return None


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


def runIgBLASTN(fasta, igdata, loci='ig', organism='human', vdb=None, ddb=None, jdb=None, cdb=None,
                output=None, format=default_igblast_output, threads=1, exec=default_igblastn_exec):
    """
    Runs igblastn on a sequence file

    Arguments:
      fasta (str): fasta file containing sequences.
      igdata (str): path to the IgBLAST database directory (IGDATA environment).
      loci (str): receptor type; one of 'ig' or 'tr'.
      organism (str): species name.
      vdb (str): name of a custom V reference in the database folder to use.
      ddb (str): name of a custom D reference in the database folder to use.
      jdb (str): name of a custom J reference in the database folder to use.
      cdb (str): name of a custom C reference in the database folder to use.
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
    # GERMLINE_C = "imgt_${SPECIES}_${RECEPTOR}_c"
    # AUXILIARY = "${SPECIES}_gl.aux"
    # IGBLAST_DB = "${IGDATA}/database"
    # IGBLAST_CMD = "igblastn \
    #     -germline_db_V ${IGBLAST_DB}/${GERMLINE_V} \
    #     -germline_db_D ${IGBLAST_DB}/${GERMLINE_D} \
    #     -germline_db_J ${IGBLAST_DB}/${GERMLINE_J} \
    #     -c_region_db ${IGBLAST_DB}/${GERMLINE_C} \
    #     -auxiliary_data ${IGDATA}/optional_file/${AUXILIARY} \
    #     -ig_seqtype ${SEQTYPE[${RECEPTOR}]} -organism ${SPECIES} \
    #     -domain_system imgt -outfmt '7 std qseq sseq btop'"
    #
    # # Set run commmand
    # OUTFILE =$(basename ${READFILE})
    # OUTFILE = "${OUTDIR}/${OUTFILE%.fasta}.fmt7"
    # IGBLAST_VER =$(${IGBLAST_CMD} -version | grep 'Package' | sed s / 'Package: ' //)
    # IGBLAST_RUN = "${IGBLAST_CMD} -query ${READFILE} -out ${OUTFILE} -num_threads ${NPROC}"

    # Define arguments
    try:
        outfmt = {'blast': '7 std qseq sseq btop', 'airr': '19'}[format]
    except KeyError:
        printError('Invalid output format %s.' % format)

    try:
        seqtype = {'ig': 'Ig', 'tr': 'TCR'}[loci]
    except KeyError:
        printError('Invalid receptor type %s.' % loci)

    # Set auxiliary data
    auxiliary = os.path.join(igdata, 'optional_file', '%s_gl.aux' % organism)
    # Set V database
    v_germ = _resolve_db(igdata, 'v', organism, loci, explicit=vdb)
    if v_germ is None:
        printError('No V database found for organism "%s" and loci "%s" with "airrc-imgt" or "imgt" prefix. '
                   'Specify --vdb, --ddb and --jdb explicitly.' % (organism, loci))
    # Set D database
    d_germ = _resolve_db(igdata, 'd', organism, loci, explicit=ddb)
    if d_germ is None:
        printError('No D database found for organism "%s" and loci "%s" with "airrc-imgt" or "imgt" prefix. '
                   'Specify --vdb, --ddb and --jdb explicitly.' % (organism, loci))
    # Set J database
    j_germ = _resolve_db(igdata, 'j', organism, loci, explicit=jdb)
    if j_germ is None:
        printError('No J database found for organism "%s" and loci "%s" with "airrc-imgt" or "imgt" prefix. '
                   'Specify --vdb, --ddb and --jdb explicitly.' % (organism, loci))
    # Set C database (optional; warn if not found rather than error)
    c_germ = _resolve_db(igdata, 'c', organism, loci, explicit=cdb)
    if c_germ is None:
        printWarning('No C database found for organism "%s" and loci "%s" with "airrc-imgt" or "imgt" prefix. '
                     'Specify --cdb explicitly if C-region annotation is required.' % (organism, loci))
        c_germ = os.path.join(igdata, 'database', 'imgt_%s_%s_c' % (organism, loci))

    # Define IgBLAST command
    cmd = [exec,
           '-query', os.path.abspath(fasta),
           '-out', os.path.abspath(output),
           '-num_threads', str(threads),
           '-ig_seqtype', seqtype,
           '-organism', organism,
           '-auxiliary_data', str(auxiliary),
           '-germline_db_V', str(v_germ),
           '-germline_db_D', str(d_germ),
           '-germline_db_J', str(j_germ),
           '-outfmt', outfmt,
           '-domain_system', 'imgt']

    # Add C-region arguments for igblastn v1.18.0
    version = getIgBLASTVersion(exec=exec)
    if parse_version(version) >= parse_version('1.18.0'):
        cmd = cmd + ['-c_region_db', str(c_germ)]
    #print(cmd)

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


def runIgBLASTP(fasta, igdata, loci='ig', organism='human', vdb=None, output=None,
                threads=1, exec=default_igblastp_exec):
    """
    Runs igblastp on a sequence file

    Arguments:
      fasta (str): fasta file containing sequences.
      igdata (str): path to the IgBLAST database directory (IGDATA environment).
      loci (str): receptor type; one of 'ig' or 'tr'.
      organism (str): species name.
      vdb (str): name of a custom V reference in the database folder to use.
      output (str): output file name. If None, automatically generate from the fasta file name.
      threads (int): number of threads for igblastp.
      exec (str): the name or path to the igblastp executable.

    Returns:
      str: IgBLAST console output.

    """
    # IGBLAST_CMD="igblastp \
    #     -germline_db_V ${IGBLAST_DB}/imgt_aa_${SPECIES}_${RECEPTOR}_v \
    #     -ig_seqtype ${SEQTYPE} -organism ${SPECIES} \
    #     -domain_system imgt -outfmt '7 std qseq sseq btop'"

    try:
        seqtype = {'ig': 'Ig', 'tr': 'TCR'}[loci]
    except KeyError:
        printError('Invalid receptor type %s.' % loci)

    # Set V database
    v_germ = _resolve_db(igdata, 'v', organism, loci, explicit=vdb, aa=True)
    if v_germ is None:
        printError('No V database found for organism "%s" and loci "%s" with "airrc-imgt_aa" or "imgt_aa" prefix. '
                   'Specify --vdb explicitly.' % (organism, loci))

    # Define IgBLAST command
    cmd = [exec,
           '-query', os.path.abspath(fasta),
           '-out', os.path.abspath(output),
           '-num_threads', str(threads),
           '-ig_seqtype', seqtype,
           '-organism', organism,
           '-germline_db_V', str(v_germ),
           '-outfmt', '7 std qseq sseq btop',
           '-domain_system', 'imgt']

    # Execute IgBLAST
    env = os.environ.copy()
    env['IGDATA'] = igdata
    try:
        stdout_str = check_output(cmd, stderr=STDOUT, shell=False, env=env,
                                  universal_newlines=True)
    except CalledProcessError as e:
        printError('Running command: %s\n%s' % (' '.join(cmd), e.output))

    return stdout_str


def getIgBLASTVersion(exec=default_igblastn_exec):
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
    match = re.search(r'(?<=Package: igblast )(\d+\.\d+\.\d+)', stdout_str)
    version = match.group(0)

    return version