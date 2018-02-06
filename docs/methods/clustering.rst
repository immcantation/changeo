.. _Clustering:

Clonal clustering methods
================================================================================

The :ref:`DefineClones` tool provides multiple different approaches to assigning
Ig sequences into clonal groups.

Clustering by V-gene, J-gene and junction length
--------------------------------------------------------------------------------

All methods provided by the :program:`bygroup` subcommand of :ref:`DefineClones`
first partition sequences based on common IGHV gene, IGHJ gene, and
junction region length. These groups are then further subdivided into
clonally related groups based on the following distance metrics on the
junction region. The specified distance metric
(:option:`--model <DefineClones bygroup --model>`) is then
used to perform hierarchical clustering under the specified linkage
(:option:`--link <DefineClones bygroup --link>`) clustering. Clonal groups are
defined by trimming the resulting dendrogram at the specified threshold
(:option:`--dist <DefineClones bygroup --dist>`).

Amino acid model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :option:`aa <DefineClones bygroup --model>` distance model is the Hamming distance
between junction amino acid sequences.

Hamming distance model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :option:`ham <DefineClones bygroup --model>` distance model is the Hamming
distance between junction nucleotide sequences.

Human and mouse 1-mer models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :option:`hh_s1f <DefineClones bygroup --model>` and
:option:`mk_rs5nf <DefineClones bygroup --model>` distance models are single
nucleotide distance matrices derived from averaging and symmetrizing the human 5-mer
targeting model in :cite:`Yaari2013` and the mouse 5-mer targeting model in
:cite:`Cui2016`. The are broadly similar to a transition/transversion model.

.. _`HH_S1F substitution matrix`:

Human 1-mer substitution matrix:

.. csv-table::
   :file: ../tables/hh_s1f_substitution.tsv
   :delim: tab
   :header-rows: 1
   :stub-columns: 1
   :widths: 15, 10, 10, 10, 10, 10

.. _`MK_RS1NF substitution matrix`:

Mouse 1-mer substitution matrix:

.. csv-table::
   :file: ../tables/mk_rs1nf_substitution.tsv
   :delim: tab
   :header-rows: 1
   :stub-columns: 1
   :widths: 15, 10, 10, 10, 10, 10

Human and mouse 5-mer models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :option:`hh_s5f <DefineClones bygroup --model>` and
:option:`mk_rs5nf <DefineClones bygroup --model>` distance models are based on
the human 5-mer targeting model in :cite:`Yaari2013` and mouse 5-mer
argeting models in :cite:`Cui2016`, respectively. The targeting
matrix :math:`T` has 5-mers across the columns and the nucleotide to
which the center base of the 5-mer mutates as the rows. The value for a
given nucleotide, 5-mer pair :math:`T[i,j]` is the product of the
likelihood of that 5-mer to be mutated :math:`mut(j)` and the
likelihood of the center base mutating to the given nucleotide
:math:`sub(j\rightarrow i)`. This matrix of probabilities is converted
into a distance matrix :math:`D` via the following steps:

#. :math:`D = -log10(T)`

#. :math:`D` is then divided by the mean of values in :math:`D`

#. All distances in :math:`D` that are infinite (probability of zero),
   distances on the diagonal (no change), and NA distances are set to 0.

Since the distance matrix :math:`D` is not symmetric, the
:option:`--sym <DefineClones bygroup --sym>` argument
can be specified to calculate either the average (avg) or minimum (min)
of :math:`D(j\rightarrow i)` and :math:`D(i\rightarrow j)`.
The distances defined by :math:`D` for each nucleotide difference are
summed for all 5-mers in the junction to yield the distance between the
two junction sequences.


Clustering by the full sequence
---------------------------------------------------------------------------------

The :program:`chen2010` and :program:`ademokun2011` methods provided by
:ref:`DefineClones` cluster sequences based on the full length sequence, by
imposing penalties for V-gene and/or J-gene mismatches.

Chen et al, 2010 method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :program:`chen2010` method of :ref:`DefineClones` is directly from
:cite:`Chen2010`. The distance metric is a normalized edit distance
(:math:`NED_VJ`) calculated as:

.. math:: NED\_VJ = \frac{LD+S_V+S_J}{L}

where :math:`LD` is the unnormalized Levenshtein distance, :math:`S_V`
is the mismatch penalty for the V-gene (0 if same gene,
1 if allele differs, 3 if gene differs, and 5 if family differs),
:math:`S_J` is the mismatch penalty for J-gene (0 if same
gene, 1 if allele differs, 3 if gene differs). :math:`L` is the CDR3
alignment length. Given this distance metric, sequences are clustered
using hierarchical clustering with average linkage. The
resulting dendrogram is trimmed at 0.32.

Ademokun et al, 2011 method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :program:`ademokun2011` method of :ref:`DefineClones` is directly
from :cite:`Ademokun2011`. The distance metric is a minimum edit distance
normalized to the length of the shorter sequence up to a maximum of 1 in 5
(or a total of 10) mismatches or indels. Distance is set to 1 for sequences
with more than the maximum number of mismatches or sequences with different
V-gene families. This metric is then used to do complete
linkage hierarchical clustering. The resulting dendrogram is trimmed at
0.25.

.. bibliography:: ../references.bib
