.. _Standard:

Data Standards
================================================================================

All Change-O tools supports both the legacy Change-O standard and the new
`Adaptive Immune Receptor Repertoire (AIRR) <https://docs.airr-community.org/en/latest/index.html>`__
standard developed by the AIRR Community (`AIRR-C <https://www.antibodysociety.org/the-airr-community/>`__).

AIRR-C Format
--------------------------------------------------------------------------------

As of v1.0.0, the default file format is the AIRR-C format as described by the
Rearrangement Schema (v1.2). The AIRR-C Rearrangement format is a tab-delimited
file format (``.tsv``) that defines the required and optional annotations for
rearranged adaptive immune receptor sequences.

To learn more about this format, the valid field names and their expected values, visit the AIRR-C
`Rearrangement Schema documentation site <https://docs.airr-community.org/en/stable/datarep/overview.html>`__.

An API for the input and output of the AIRR-C format is provided in the
`AIRR Python package <https://docs.airr-community.org/en/stable/packages/airr-python/overview.html>`__.
Wrappers for this package are provided in the API as :class:`changeo.IO.AIRRReader`
and :class:`changeo.IO.AIRRWriter`.

Change-O Format
--------------------------------------------------------------------------------

The legacy Change-O standard is a tab-delimited file format (``.tab``) with a set
of predefined column names. The standardized column names used by the Change-O format
are shown in the table below. Most tools do not require every column. The columns
required by and added by each individual tool are described in the
:ref:`commandline usage <Usage>` documentation. If a column contains multiple
entries, such as ambiguous V gene assignments, these nested entries are delimited
by commas. The ordering of the columns does not matter.

An API for the input and output of the Change-O format is provided in
:class:`changeo.IO.ChangeoReader` and :class:`changeo.IO.ChangeoWriter` respectively.

.. csv-table::
   :file: tables/column_descriptions.tsv
   :delim: tab
   :header-rows: 1
   :widths: 15, 15, 10, 60
