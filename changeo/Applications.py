"""
Application wrappers
"""

# Info
__author__ = 'Jason Anthony Vander Heiden'

# Imports
import os
import sys
from subprocess import check_output, STDOUT, CalledProcessError

# Presto and changeo imports
from changeo.Defaults import default_igblast_exec, default_tbl2asn_exec, default_igphyml_exec


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

    # Define tb2asn command
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
        stdout_str = check_output(cmd, stderr=STDOUT, shell=False,
                                  universal_newlines=True, cwd=rep_dir)
    except CalledProcessError as e:
        sys.stderr.write('\nError running command: %s\n' % ' '.join(cmd))
        sys.exit(e.output)

    return None


def runIgBLAST(fasta, igdata, locus='ig', organism='human', output=None,
               threads=1, exec=default_igblast_exec):
    """
    Runs IgBLAST on a sequence file

    Arguments:
      fasta (str): fasta file containing sequences.
      igdata (str): path to the IgBLAST database directory (IGDATA environment).
      locus (str): receptor type; one of 'ig' or 'tr'.
      organism (str): species name.
      output (str): output file name. If None, automatically generate from the fasta file name.
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
        seqtype = {'ig': 'Ig', 'tr': 'TCR'}[locus]
    except KeyError:
        sys.exit('Error: Invalid receptor type %s' % locus)

    v_germ = os.path.join([igdata, 'database', 'imgt_%s_%s_v' % (organism, locus)])
    d_germ = os.path.join([igdata, 'database', 'imgt_%s_%s_v' % (organism, locus)])
    j_germ = os.path.join([igdata, 'database', 'imgt_%s_%s_v' % (organism, locus)])
    auxilary = os.path.join([igdata, 'optional_file', '%s_gl.aux' % organism])

    # Define IgBLAST command
    cmd = [exec,
           '-query', os.path.abspath(fasta),
           '-out', os.path.abspath(output),
           '-num_threads', str(threads),
           '-germline_db_V', str(v_germ),
           '-germline_db_D', str(d_germ),
           '-germline_db_J', str(j_germ),
           '-auxiliary_data', str(auxilary),
           '-ig_seqtype', seqtype,
           '-organism', organism,
           '-domain_system', 'imgt',
           '-outfmt', '"7 std qseq sseq btop"']

    # Execute IgBLAST
    try:
        stdout_str = check_output(cmd, stderr=STDOUT, shell=False,
                                  universal_newlines=True)
    except CalledProcessError as e:
        sys.stderr.write('\nError running command: %s\n' % ' '.join(cmd))
        sys.exit(e.output)

    #if 'Unable to read any FASTA records' in stdout_str:
    #    sys.stderr.write('\n%s failed: %s\n' % (' '.join(cmd), stdout_str))

    return stdout_str