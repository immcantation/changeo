Parsing IMGT output
================================================================================

Reducing file size for submission to IMGT/HighV-QUEST
--------------------------------------------------------------------------------

`Via pRESTO <http://presto.readthedocs.org/en/latest/tasks.html#reducing-file-size-for-submission-to-imgt-highv-quest>`__

.. todo::

.. code-block:: bash
    :linenos:

    CollapseSeq
    SplitSeq count

Processing the output of IMGT/HighV-QUEST
--------------------------------------------------------------------------------

.. todo::

.. code-block:: bash
    :linenos:

	#!/bin/bash
	# Parse IMGT results
	# 
	# Author:  Jason Anthony Vander Heiden, Namita Gupta
	# Date:    2015.12.04
	# 
	# Required Arguments:
	#   $1 = IMGT output file 
	#   $2 = FASTA file that was submitted to IMGT
	#   $3 = the folder(s) or file(s) containing the IMGT reference database sequences


	# Capture command line parameters
	IMGT_FILE=$(readlink -f $1)
	SEQ_FILE=$(readlink -f $2)
	GERM_DIR=$(readlink -f $3)

	# Parse IMGT output
	MakeDb.py imgt -i $IMGT_FILE -s $SEQ_FILE


