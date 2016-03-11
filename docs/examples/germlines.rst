Building germline sequences
================================================================================

:ref:`CreateGermlines` is used to infer the germline V(D)J receptor sequence
from which mutations can be inferred. The alignment information is saved
by :ref:`MakeDb`, but the germline gene segment sequences that were used for
the alignment must be passed on via :option:`-r <CreateGermlines -r>`. Because the D call for B cell receptor alignments
is often low confidence, the default germline format :option:`-g dmask <CreateGermlines -g>`
puts Ns in the junction region rather than the inferred germline D region. If
running :ref:`CreateGermlines` after :ref:`DefineClones`, the :option:`--cloned <CreateGermlines --cloned>`
flag will generate a single germline of consensus length for each clone. The
`Change-O example files <http://clip.med.yale.edu/changeo/rtd/Changeo_Example.tar.gz>`__
contain the IMGT germline fasta files as well as a database file on which to run
the command below::

    CreateGermlines.py -d S43_db-pass_parse-select.tab -r IMGT_Human_IGH[VDJ].fasta