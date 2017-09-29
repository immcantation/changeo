.. _IgBLAST:

Parsing IgBLAST output
================================================================================

Example Data
--------------------------------------------------------------------------------

We have hosted a small example data set resulting from the
`Roche 454 example workflow <http://presto.readthedocs.io/en/latest/workflows/Jiang2013_Workflow.html>`__
described in the `pRESTO <http://presto.readthedocs.io>`__ documentation. In addition to the
example FASTA files, we have included the standalone `IgBLAST <http://www.ncbi.nlm.nih.gov/projects/igblast/faq.html#standalone>`__
results. The files can be downloded from here:

`Change-O Example Files <http://clip.med.yale.edu/immcantation/examples/Changeo_Example.tar.gz>`__

Building the IgBLAST database
--------------------------------------------------------------------------------

To ensure that the receptor sequences can be manipulated to follow
the IMGT numbering scheme, enabling identification of the junction region as
well as framework and complementarity-determining regions, it is essential to
build the standalone `IgBLAST <http://www.ncbi.nlm.nih.gov/projects/igblast/faq.html#standalone>`__
reference database from an originally IMGT numbered reference.
The first step is to download the requisite reference database from any source
with IMGT numbered germlines in FASTA format
(`IMGT reference directory <http://imgt.org/vquest/refseqh.html>`__).
These FASTA files must then be converted into IgBLAST databases. The
following example demonstrates how to make the IgBLAST database for the
included reference fasta files (Human B cell receptor heavy chains)
that should be placed in the IgBLAST directory
(eg, ``/home/user/apps/igblast-1.7.0/bin``).

.. literalinclude:: scripts/IgBLAST_Commands.sh
   :language: none
   :linenos:
   :lines: 2-10

Once these databases are built for each segment of the receptor, they can be
referenced when running IgBLAST.

Running IgBLAST Standalone
--------------------------------------------------------------------------------

To run standalone `IgBLAST <http://www.ncbi.nlm.nih.gov/projects/igblast/faq.html#standalone>`__
with the example code below, place the example input FASTA
file into the IgBLAST directory (eg, ``/home/user/apps/igblast-1.7.0/bin``).
The example data is human B cell receptor heavy chains, but the parameters
can be changed according to analyze other types of receptors.
To ensure the requisite information is generated you must specify the special
output format shown below using the ``-outfmt '7 std qseq sseq btop'`` argument:

.. literalinclude:: scripts/IgBLAST_Commands.sh
   :language: none
   :linenos:
   :lines: 12-20

.. seealso::

    For this example, we have suggested moving the FASTA file into the IgBLAST
    directory for the sake of simplicity. See the
    `IgBLAST documentation <http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/scripts/projects/igblast/README>`__
    regarding use of the ``IGDATA`` environment variable to control where IgBLAST
    looks for its internal database files. This will allow you run IgBLAST from any location.

Configuring and Running IgBLAST Using the Wrapper Scripts
--------------------------------------------------------------------------------

A collection of scripts for setting up and running IgBLAST are available in the
`Immcantation repository <https://bitbucket.org/kleinstein/immcantation/src/tip/scripts>`__.
To use the IgBLAST wrapper, ``run_igblast.sh``, copy all the tools in the ``/scripts`` folder
to a location in your ``PATH``. Download and configure the IgBLAST and IMGT reference databases
as follows:

.. literalinclude:: scripts/IgBLAST_Wrapper.sh
   :language: none
   :linenos:
   :lines: 2-6

Then run IgBLAST as follows:

.. literalinclude:: scripts/IgBLAST_Wrapper.sh
   :language: none
   :linenos:
   :lines: 7-8

Processing the output of IgBLAST
--------------------------------------------------------------------------------

Standalone `IgBLAST <http://www.ncbi.nlm.nih.gov/projects/igblast/faq.html#standalone>`__
output is parsed by the :program:`igblast` subcommand of
:ref:`MakeDb` to generate the standardized tab-delimited database file on which all
subsequent Change-O modules operate. In addition to the IgBLAST output
(:option:`-i S43_atleast-2.fmt7  <MakeDb igblast -i>`), both the FASTA files input to
IgBLAST (:option:`-s S43_atleast-2.fasta  <MakeDb igblast -s>`) and the IMGT-gapped reference
sequences (:option:`-r Human_IGH[VDJ].fasta <MakeDb igblast -r>`) must be provided to :ref:`MakeDb`:

.. literalinclude:: scripts/IgBLAST_Commands.sh
   :language: none
   :lines: 22-23

The optional (:option:`--regions <MakeDb igblast --regions>`) and
(:option:`--scores <MakeDb igblast --scores>`) arguments add extra columns to the output
database containing IMGT-gapped CDR/FWR regions and alignment metrics, respectively.

.. warning::

    The references sequences you provide to :ref:`MakeDb` must contain IMGT-gapped
    V-segment references, and these reference must be the same sequences used to
    build the IgBLAST reference database. If your IgBLAST germlines are not IMGT-gapped
    and/or they are not identical to those provided to :ref:`MakeDb`, sequences assigned
    missing germlines will fail the parsing operation and the CDR3/junction sequences
    will not be correct.
