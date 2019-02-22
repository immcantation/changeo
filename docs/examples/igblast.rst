.. _IgBLAST:

Using IgBLAST
================================================================================

Example data
--------------------------------------------------------------------------------

We have hosted a small example data set resulting from the
`Roche 454 example workflow <http://presto.readthedocs.io/en/latest/workflows/Jiang2013_Workflow.html>`__
described in the `pRESTO <http://presto.readthedocs.io>`__ documentation. In addition to the
example FASTA files, we have included the standalone `IgBLAST <https://ncbi.github.io/igblast>`__
results. The files can be downloded from here:

`Change-O Example Files <http://clip.med.yale.edu/immcantation/examples/Changeo_Example.tar.gz>`__

Configuring IgBLAST
--------------------------------------------------------------------------------

A collection of scripts for setting up the standalone IgBLAST database from the
IMGT reference sequences are available on the
`Immcantation repository <https://bitbucket.org/kleinstein/immcantation/src/tip/scripts>`__.
To use these scripts, copy all the tools in the ``/scripts`` folder to a location
in your ``PATH``. At a minimum, you'll need the following scripts:

1. ``fetch_igblastdb.sh``
2. ``fetch_imgtdb.sh``
3. ``clean_imgtdb.py``
4. ``imgt2igblast.sh``

Download and configure the IgBLAST and IMGT reference databases
as follows:

.. literalinclude:: scripts/IgBLAST_Wrapper.sh
   :language: none
   :linenos:
   :lines: 2-6

.. note::

    Several Immcantation tools require the observed V(D)J sequence
    (``SEQUENCE_IMGT``) and associated germline fields (``GERMLINE_IMGT``
    or ``GERMLINE_IMGT_D_MASK``) to have gaps inserted to conform to the
    IMGT numbering scheme. Thus, when a tool such as :ref:`MakeDb` or
    :ref:`CreateGermlines` requires a reference sequence set as input,
    it will required the IMGT-gapped reference set. Meaning,
    the reference sequences that were downloaded using the ``fetch_imgtdb.sh``
    script, or downloaded manually from the
    `IMGT reference directory <http://imgt.org/vquest/refseqh.html>`__,
    rather than the final upgapped reference set required by IgBLAST.

.. seealso::

    The provided scripts download only the mouse and human IMGT reference databases.
    See the `IgBLAST documentation <https://ncbi.github.io/igblast>`__ for instructions
    on how to build the database in a more general case. Shown below is an example of how
    to performed the same steps as the Immcantation scripts using a separately
    downloaded IMGT reference set and the scripts provided by IgBLAST. You must have all of
    the associated commands in your ``PATH`` and the appropriate directories created:

    .. literalinclude:: scripts/IgBLAST_Commands.sh
       :language: none
       :linenos:
       :lines: 2-13

    Once these databases are built for each segment they can be referenced when
    running IgBLAST.

Running IgBLAST
--------------------------------------------------------------------------------

Change-O provides a simple wrapper script to run IgBLAST with the required
options as the :program:`igblast` subcommand of :ref:`AssignGenes`. This wrapper
can be run as follows using the database built using the Immcantation scripts::

    AssignGenes.py igblast -s S43_atleast-2.fasta -b ~/share/igblast \
        --organism human --loci ig --format blast

The optional :option:`--format blast <AssignGenes igblast --format>` argument
defines the output format of IgBLAST. The default, ``blast``, is the
blocked tabular output provided by specifying the ``-outfmt '7 std qseq sseq btop'``
argument to IgBLAST. Specifying :option:`--format airr <AssignGenes igblast --format>`
will output a tab-delimited file compliant with the
`AIRR Rearrangement schema <https://airr-standards.readthedocs.io/en/latest/datarep/rearrangements.html>`__
defined by the `AIRR Community <https://www.antibodysociety.org/about-the-airr-community>`__.
AIRR format support requires IgBLAST v1.9.0 or higher.

The :option:`-b ~/share/igblast <AssignGenes igblast -b>` argument specifies the
path containing the ``database``, ``internal_data``, and ``optional_file``
directories required by IgBLAST. This option sets the ``IGDATA`` environment variable
that controls where IgBLAST looks for internal database files. See the
`IgBLAST documentation <https://ncbi.github.io/igblast>`__ for more details
regarding the ``IGDATA`` environment variable.

.. seealso::

    The :ref:`AssignGenes` IgBLAST wrapper provides limited functionality.
    For more control, IgBLAST should be run directly. The only strict
    requirement for compatibility with Changeo-O is that the output must
    either be an AIRR tab-delimited file (``--outfmt 19``) or a blast-style
    tabular output with the optional query sequence, subject sequence and BTOP fields
    (``-outfmt '7 std qseq sseq btop'``). An example of how to run IgBLAST
    directly is shown below:

    .. literalinclude:: scripts/IgBLAST_Commands.sh
       :language: none
       :linenos:
       :lines: 15-24

Processing the output of IgBLAST
--------------------------------------------------------------------------------

Standalone IgBLAST blast-style tabular output is parsed by the :program:`igblast` subcommand of :ref:`MakeDb`
to generate the standardized tab-delimited database file on which all subsequent Change-O modules operate.
In addition to the IgBLAST output (:option:`-i S43_atleast-2.fmt7  <MakeDb igblast -i>`),
both the FASTA files input to IgBLAST (:option:`-s S43_atleast-2.fasta  <MakeDb igblast -s>`)
and the IMGT-gapped reference sequences
(:option:`-r IMGT_Human_IGHV.fasta IMGT_Human_IGHD.fasta IMGT_Human_IGHJ.fasta <MakeDb igblast -r>`)
must be provided to :ref:`MakeDb`::

    MakeDb.py igblast -i S43_atleast-2.fmt7 -s S43_atleast-2.fasta \
        -r IMGT_Human_IGHV.fasta IMGT_Human_IGHD.fasta IMGT_Human_IGHJ.fasta \
        --regions --scores

The optional :option:`--regions <MakeDb igblast --regions>` and
:option:`--scores <MakeDb igblast --scores>` arguments add extra columns to the output
database containing IMGT-gapped CDR/FWR regions and alignment metrics, respectively.

.. warning::

    The references sequences you provide to :ref:`MakeDb` must contain IMGT-gapped
    V segment references, and these reference must be the same sequences used to
    build the IgBLAST reference database. If your IgBLAST germlines are not IMGT-gapped
    and/or they are not identical to those provided to :ref:`MakeDb`, then sequences
    which were assigned missing germlines will fail the parsing operation and the
    junction (CDR3) sequences will not be correct.
