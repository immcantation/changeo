Parsing IgBLAST output
================================================================================

Building the IgBLAST database
--------------------------------------------------------------------------------

.. todo::

.. code-block:: bash
    :linenos:

    getIMGT.sh
    imgt2igblast.sh

Running IgBLAST
--------------------------------------------------------------------------------

.. todo::

.. code-block:: bash
    :linenos:

    # Define IgBLAST command
    IGDATA="${HOME}/apps/igblast-1.4.0"
    IGBLAST_DB="${IGDATA}/database"
    IGBLAST_CMD="${IGDATA}/igblastn \
    -germline_db_V ${IGBLAST_DB}/imgt_human_IG_V \
    -germline_db_D ${IGBLAST_DB}/imgt_human_IG_D \
    -germline_db_J ${IGBLAST_DB}/imgt_human_IG_J \
    -auxiliary_data ${IGDATA}/optional_file/human_gl.aux \
    -domain_system imgt -ig_seqtype Ig -organism human \
    -outfmt '7 std qseq sseq btop'"

    # Align V(D)J segments using IgBLAST
    IGBLAST_RUN="${IGBLAST_CMD} -query ${OUTNAME}_collapse-unique.fasta \
        -out ${OUTNAME}_collapse-unique.fmt7
    eval $IGBLAST_RUN


Processing the output of IgBLAST
--------------------------------------------------------------------------------

.. code-block:: bash
    :linenos:

	#!/bin/bash
	# Parse IgBLAST results
	#
	# Author:  Jason Anthony Vander Heiden, Namita Gupta
	# Date:    2015.12.04
	#
	# Required Arguments:
	#   $1 = IgBLAST output file
	#   $2 = FASTA file that was submitted to IgBLAST
	#   $3 = the folder(s) or file(s) containing the IMGT-gapped reference sequences


	# Capture command line parameters
	IGBLAST_FILE=$(readlink -f $1)
	SEQ_FILE=$(readlink -f $2)
	REF_GAPPED=$(readlink -f $3)

	# Parse IMGT output
	MakeDb.py igblast -i  IGBLAST_FILE \
    -s  SEQ_FILE -r $REF_GAPPED \
    --scores --regions
