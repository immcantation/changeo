.. _Filtering:

Filtering records
================================================================================

The :ref:`ParseDb` tool provides a basic set of operations for manipulating
Change-O database files from the commandline, including removing or updating
rows and columns.

.. _Filtering-Functional:

Removing non-productive sequences
--------------------------------------------------------------------------------

After building a Change-O database from either :ref:`IMGT/HighV-QUEST <IMGT>` or
:ref:`IgBLAST <IgBLAST>` output, you may wish to subset your data to only
productive sequences. This can be done in one of two roughly equivalent ways
using the :ref:`ParseDb` tool:

.. code-block:: none
    :linenos:

    ParseDb.py select -d HD13M_db-pass.tsv -f productive -u T
    ParseDb.py split -d HD13M_db-pass.tsv -f productive

The first line above uses the :program:`select` subcommand to output a single file
labeled ``parse-select`` containing only records with the value of ``T``
(:option:`-u T <ParseDb select -u>`) in the ``productive`` column
(:option:`-f productive <ParseDb select -f>`).

Alternatively, the second line above uses the :program:`split` subcommand to output
multiple files with each file containing records with one of the values found in the
``productive`` column (:option:`-f productive <ParseDb split -f>`). This will
generate two files labeled ``productive-T`` and ``productive-F``.

Removing disagreements between the C-region primers and the reference alignment
--------------------------------------------------------------------------------

If you have data that includes both heavy and light chains in the same library,
the V-segment and J-segment alignments from IMGT/HighV-QUEST or IgBLAST may not
always agree with the isotype assignments from the C-region primers. In these cases,
you can filter out such reads with the :program:`select` subcommand of :ref:`ParseDb`.
An example function call using an imaginary file ``db.tsv`` is provided below:

.. code-block:: none
    :linenos:

    ParseDb.py select -d db.tsv -f v_call j_call c_call -u "IGH" \
        --logic all --regex --outname heavy
    ParseDb.py select -d db.tsv -f v_call j_call c_call -u "IG[LK]" \
        --logic all --regex --outname light

These commands will require that all of the ``v_call``, ``j_call`` and ``c_call``
fields (:option:`-f v_call j_call c_call <ParseDb select -f>` and
:option:`--logic all <ParseDb select --logic>`) contain the string ``IGH`` (lines 1-2)
or one of ``IGK`` or ``IGL`` (lines 3-4). The :option:`--regex <ParseDb select --regex>`
argument allows for partial matching and interpretation of regular expressions. The
output from these two commands are two files, one containing only heavy chains
(``heavy_parse-select.tsv``) and one containg only light chains (``light_parse-select.tsv``).


Exporting records to FASTA files
--------------------------------------------------------------------------------

You may want to use external tools, or tools from `pRESTO <presto.readthedocs.io>`__,
on your Change-O result files. The :ref:`ConvertDb` tool provides two options for
exporting data from tab-delimited files to FASTA format.

Standard FASTA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :program:`fasta` subcommand allows you to export sequences and annotations to
FASTA formatted files in the
`pRESTO annototation scheme <http://presto.readthedocs.io/en/stable/overview.html#annotation-scheme>`__::

    ConvertDb.py fasta -d HD13M_db-pass.tsv --if sequence_id \
        --sf sequence_alignment --mf v_call duplicate_count

Where the column containing the sequence identifier is specified by
:option:`--if sequence_id <ConvertDb fasta --if>`, the nucleotide sequence column is
specified by :option:`--sf sequence_id <ConvertDb fasta --sf>`, and additional annotations
to be added to the sequence header are specified by
:option:`--mf v_call duplicate_count <ConvertDb fasta --mf>`.

BASELINe FASTA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :program:`baseline` subcommand generates a FASTA derivative format required by the
`BASELINe <http://selection.med.yale.edu/baseline>`__ web tool. Generating these
files is similar to building standard FASTA files, but requires a few more options.
An example function call using an imaginary file ``db.tsv`` is provided below::

    ConvertDb.py baseline -d db.tsv --if sequence_id \
        --sf sequence_alignment --mf v_call duplicate_count \
        --cf clone_id --gf germline_alignment_d_mask

The additional arguments required by the :program:`baseline` subcommand include the
clonal grouping (:option:`--cf clone_id <ConvertDb baseline --cf>`) and germline sequence
(:option:`--gf germline_alignment_d_mask <ConvertDb baseline --gf>`) columns added by
the :ref:`DefineClones <Cloning>` and :ref:`CreateGermlines <Germlines>` tasks,
respectively.

.. note::

    The :program:`baseline` subcommand requires the ``CLONE`` column to be sorted.
    :ref:`DefineClones` generates a sorted ``CLONE`` column by default. However,
    you needed to alter the order of the ``CLONE`` column at some point,
    then you can re-sort the clonal assignments using the :program:`sort`
    subcommand of :ref:`ParseDb`. An example function call using an imaginary 
    file ``db.tsv`` is provided below::

        ParseDb.py sort -d db.tsv -f clone_id

    Which will sort records by the value in the ``clone_id`` column
    (:option:`-f clone_id <ParseDb sort -f>`).
