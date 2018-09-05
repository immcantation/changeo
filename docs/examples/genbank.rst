.. _GenBank:

Generating MiAIRR compliant GenBank/TLS submissions
================================================================================

MiAIRR
--------------------------------------------------------------------------------

The MiAIRR standard (minimal information about adaptive immune receptor repertoires)
is a minimal reporting standard for experiments using sequencing-based technologies
to study adaptive immune receptors (T and B cell receptors). The current version (1.0)
of the standard was published in `Rubelt et al, 2017 <https://www.nature.com/articles/ni.3873>`__
and accepted by the general assembly at the annual AIRR Community meeting in December 2017.

MiAIRR recommends submission of raw read data to the `Sequence Read Archive (SRA) <https://www.ncbi.nlm.nih.gov/sra>`__
and submission of processed and annotated data to the
`Targeted Locus Study (TLS) <https://www.ncbi.nlm.nih.gov/genbank/tls>`__ section of
`GenBank <https://www.ncbi.nlm.nih.gov/genbank>`__.

This example will cover generation of files for submission to TLS starting from Change-O
formatted data.  For complete details of the required and optional elements of the TLS submission
see the `AIRR Standards <http://docs.airr-community.org/en/latest/miairr/miairr_ncbi_overview.html>`__
documentation site.

Example Data
--------------------------------------------------------------------------------
We have hosted a small example data set resulting from the
`Roche 454 example workflow <http://presto.readthedocs.io/en/latest/workflows/Jiang2013_Workflow.html>`__
described in the `pRESTO <http://presto.readthedocs.io>`__ documentation. The files can be
downloded from here:

`Change-O Example Files <http://clip.med.yale.edu/immcantation/examples/Changeo_Example.tar.gz>`__

The following examples use the ``S43_db-pass_parse-select.tab`` database file and
``S43_template.sbt`` file provided in the example bundle, which has already undergone
the :ref:`IgBLAST <IgBLAST>` annotation, parsing, and :ref:`filtering <Filtering-Functional>`
operations.

Generating files for submission
--------------------------------------------------------------------------------

Requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    + An annotated data set in either the :ref:`Change-O <Standard>` or the
      `AIRR Data Representation <http://docs.airr-community.org/en/latest/datarep/rearrangements.html>`__
      formats. Records must have valid V, J and junction region annotations to be suitable for submission.
    + `tbl2asn <https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2>`__ installed and in your ``PATH``.
    + A GenBank submission template file (``.sbt``), generating using the
      `NCBI Template Generator <https://submit.ncbi.nlm.nih.gov/genbank/template/submission>`__.

Cleaning up C region annotations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, we need to rename the C primer annotations to valid constant region gene names so they are
recognized by remote databases. We'll use the :program:`update` subcommand of
:ref:`ParseDb` to rename the values in the ``CPRIMER`` column of the example file and save that
file as ``S43_update.tab``::

    ParseDb.py update -d S43_db-pass_parse-select.tab -o S43_update.tab \
        -f CPRIMER \
        -u IgA-PCR IgD-PCR IgE-PCR IgG-PCR IgM-PCR \
        -t IGHA IGHD IGHE IGHG IGHM

Creating ASN files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ASN submission files are generated using the :program:`genbank` subcommand of
:ref:`ConvertDb` as follows::

    ConvertDb.py genbank -d S43_update.tab \
        --product "immunoglobulin heavy chain" \
        --inf "IgBLAST:1.7.0" \
        --organism "Homo sapiens" \
        --tissue "Peripheral blood" \
        --cell-type "B cell" \
        --isolate S43 \
        --cf CPRIMER \
        --nf DUPCOUNT \
        --asis-id \
        --asn \
        --sbt S43_template.sbt \
        --outdir S43_TLS

The resulting output in the ``S43_TLS`` folder will include a number of files.
The Sequin file ``S43_update_genbank.sqn`` is the file that will be used for submission
and the GenBank record file ``S43_update_genbank.gbf`` is similar to what the submission
will look like once it has been accepted by GenBank.

The command above manually specifies several required and optional annotations.
Alternatively, sample information (``organism``, ``sex``, ``isolate``, ``tissue_type``,
``cell_type``) can be specified in a separate yaml file and provided via the
:option:`-y <ConvertDb genbank -y>` argument. Additional harmonized
`BioSample attributes <https://www.ncbi.nlm.nih.gov/biosample/docs/attributes>`__,
which are not convered by the existing commandline arguments, may be provided
in the yaml file. Note, the yaml file adds only sample features, so it cannot be used
to specify source features (``--product``, ``--mol``, ``--inf`` and ``--db`` arguments), parsing
arguments, or run parameters (`--label``, ``--exec``, etc). Features specified in the yaml
file will override equivalent features specified through the corresponding commandline arguments.

.. note::

    The example shown above automatically runs :program:`tbl2asn`, because the
    :option:`--asn <ConvertDb genbank --asn>` argument was specified. :ref:`ConvertDb`
    can be run without running :program:`tbl2asn`, which will generate only the
    feature table (``S43_update_genbank.tbl``) and fasta (``S43_update_genbank.fsa``) files
    required to run :program:`tbl2asn` manually via the command::

        tbl2asn -p . -a s -V vb -t S43_template.sbt

.. warning::

    There is a known issue with the :option:`--asn <ConvertDb genbank --asn>` argument.
    In some environments, for reasons that are presently unknown, :program:`tbl2asn`
    may fail to recongizing the input fasta file and report an error stating
    **Unable to read any FASTA records**. Running :program:`tbl2asn` manually should
    resolve the issue.


Submitting to GenBank/TLS using SequinMacroSend
--------------------------------------------------------------------------------

After generating the ``.sqn`` files, you can submit them as MiAIRR compliant
GenBank/TLS records using GenBank's
`SequinMacroSend <https://www.ncbi.nlm.nih.gov/LargeDirSubs/dir_submit.cgi>`__ service.

When submitting, simply add the keyword ``AIRR`` to the subject line in the
submission system and it will be routed accordingly.

.. warning::

    Currently, the SequinMacroSend system cannot accept files over 512MB in size.
    For submissions over the size limit, you must split them into smaller files
    and note in the submission comments that they are a part of a split submission.
    See the :program:`split` subcommand of  :ref:`ParseDb` for one approach to
    logically dividing large submissions.