.. _Standard:

Data Standard
================================================================================

All Change-O tools supports both the classical Change-O standard and the new
`Adaptive Immune Receptor Repertoire (AIRR) <https://docs.airr-community.org/en/latest/index.html>`__
standard developed by the AIRR Community (`AIRR-C <https://www.antibodysociety.org/the-airr-community/>`__).

As of v1.0.0 the default file format will be the AIRR-C format, as described by the Rearrangement 
Schema. To learn about this format, the valid field names and their expected values, visit the 
AIRR-C `Rearrangement Schema documentation site <https://docs.airr-community.org/en/latest/datarep/rearrangements.html>`__.
    
The Change-O standard is a tab-delimited file format (``.tab``) with a set of predefine
column names. The standardized column names used by the Change-O format are shown in the table below.
Most tools do not require every column. The columns required by, and added by, each
individual tool are described in the :ref:`commandline usage <Usage>` documentation.
If a columns contains multiple entries, such as ambiguous V
segment assignments, these nested extries are delimited by commas.
The ordering of the columns does not matter.

An API for input and output of the Change-O format is provide in
:class:`changeo.IO.ChangeoReader` and :class:`changeo.IO.ChangeoWriter`, respectively.
See :class:`changeo.Receptor.Receptor` for a description of the data structure used by
the reader and writer classes.

.. csv-table::
   :file: tables/column_descriptions.tsv
   :delim: tab
   :header-rows: 1
   :widths: 15, 15, 10, 75
