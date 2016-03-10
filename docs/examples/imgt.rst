Parsing IMGT output
================================================================================

Example Data
--------------------------------------------------------------------------------

We have hosted a small example dataset resulting from the pre-processing
done in a `pRESTO example workflow <http://presto.readthedocs.org/en/latest/workflows/Stern2014_Workflow.html>`__.
In addition to the sample data fasta file, we include `IMGT/HighV-QUEST <http://imgt.org/HighV-QUEST>`__
results. The files can be downloded from here:

`MakeDb example files <http://clip.med.yale.edu/changeo/rtd/MakeDb_Example.tar.gz>`__

Reducing file size for submission to IMGT/HighV-QUEST
--------------------------------------------------------------------------------

`IMGT/HighV-QUEST <http://imgt.org/HighV-QUEST>`__ currently limits the size of
uploaded files to 500,000 sequences. To accomodate this limit, you can use
the :program:`count` subcommand of `SplitSeq <http://presto.readthedocs.org/en/latest/tools/SplitSeq.html#splitseq>`__ to divide your files into
small pieces.

.. code-block:: none

    SplitSeq.py count -s file.fastq -n 500,000 --fasta

The ``-n 500,000`` argument sets the maximum number of
sequences in each file and the ``--fasta``
tells the tool to output a FASTA, rather than FASTQ, formatted file.

.. note::

    You can usually avoid the necessity of reducing file sizes by removing
    duplicate sequences first using the `CollapseSeq <http://presto.readthedocs.org/en/latest/tools/CollapseSeq.html#collapseseq>`__ tool.

.. seealso::

    `pRESTO <http://presto.readthedocs.org/en/latest/tasks.html#reducing-file-size-for-submission-to-imgt-highv-quest>`__


Processing the output of IMGT/HighV-QUEST
--------------------------------------------------------------------------------

Standard output from `IMGT/HighV-QUEST <http://imgt.org/HighV-QUEST>`__ can be
parsed via the :program:`imgt` subcommand of :ref:`MakeDb` to form the standard
tab-delimited database file on which any subsequent Change-O module operates.
This parsing requires the fasta file (:option:`-s <MakeDb imgt -s>`)
used as input to HighV-QUEST and either the compressed output file or the
uncompressed folder with the ``1_Summary``, ``2_IMGT-gapped``, ``3_Nt-sequences``, and
``6_Junction`` files (:option:`-i <MakeDb imgt -i>`). In this example, the standard
information is parsed in addition to the IMGT-defined framework and
complementarity-determining regions (:option:`--regions <MakeDb imgt --regions>`).

.. code-block:: none

   MakeDb.py imgt -s MS12_atleast-2.fasta -i MS12_atleast-2.txz --regions


