Overview
================================================================================

Change-O performs analyses of lymphocyte receptor sequences following alignment
against the germline reference. It includes tools for standardizing the
output of alignment software, clonal assignment, germline reconstruction, and
basic database manipulation. Change-O was designed to be simple to use, but it
does require some familiarity with commandline applications.
To maximize flexibility, Change-O employs a simple tab-delimited database format
with :ref:`standardized column names <Standard>`, allowing easy use of Change-O
output with external environments and interoperability with the
`Alakazam <http://kleinstein.bitbucket.org/alakazam>`_,
`Shazam <http://kleinstein.bitbucket.org/shazam>`_, and
`TIgGER <http://kleinstein.bitbucket.org/tigger>`_ R packages.
A brief description of each tool is shown in the table below.

.. _FeatureTable:

.. csv-table::
   :file: tools/tool_summary.tsv
   :delim: tab
   :header-rows: 1
   :widths: 15, 5, 75