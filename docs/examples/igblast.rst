.. _IgBLAST:

Parsing IgBLAST output
================================================================================

Example Data
--------------------------------------------------------------------------------

We have hosted a small example dataset resulting from the pre-processing
done in a `pRESTO example workflow <http://presto.readthedocs.org/en/latest/workflows/Stern2014_Workflow.html>`__.
In addition to the sample data fasta file, we include IMGT reference germline
files for running IgBLAST on the data. The files can be downloded from here:

`Change-O example files <http://clip.med.yale.edu/changeo/rtd/Changeo_Example.tar.gz>`__

Building the IgBLAST database
--------------------------------------------------------------------------------

In order to ensure that the receptor sequences can be manipulated to follow
the IMGT numbering scheme, enabling identification of the junction region as
well as framework and complementarity-determining regions, it is essential to
build the IgBLAST reference database from an originally IMGT numbered reference.
The first step is to download the requisite reference database from any source
with IMGT numbered germlines in fasta format
(`IMGT reference directory <http://imgt.org/vquest/refseqh.html>`__).
These fasta files must then be converted into IgBLAST databases. The
following example demonstrates how to make the IgBLAST database for the
included reference fasta files (Human B cell receptor heavy chains)
that should be placed in the IgBLAST directory
(e.g., ``/home/user/apps/igblast-1.4.0/bin``).

.. literalinclude:: scripts/IgBLAST_Commands.sh
   :language: none
   :linenos:
   :lineno-match:
   :lines: 1-8

Once these databases are built for each segment of the receptor, they can be
referenced when running IgBLAST.

Running IgBLAST
--------------------------------------------------------------------------------

To run IgBLAST with the example code below, place the example input fasta
file into the IgBLAST directory (e.g., ``/home/user/apps/igblast-1.4.0/bin``).
The example data is human B cell receptor heavy chains, but the parameters
can be changed according to IgBLAST help to analyze other types of receptors.
To ensure the requisite information is generated a special output format
is specified using the ``-outfmt`` parameter:

.. literalinclude:: scripts/IgBLAST_Commands.sh
   :language: none
   :linenos:
   :lineno-match:
   :lines: 11-19


Processing the output of IgBLAST
--------------------------------------------------------------------------------

IgBLAST output can be parsed via the :program:`igblast` subcommand of
:ref:`MakeDb` to form the standard tab-delimited database file on which any
subsequent Change-O module operates. The original gapped reference fasta files
must be provided to recreate the IMGT numbered sequences and junction region
following the IMGT definition for further analyses:

.. literalinclude:: scripts/IgBLAST_Commands.sh
   :language: none
   :linenos:
   :lineno-match:
   :lines: 21
