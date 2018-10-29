IgPhyML lineage tree analysis
===============================

This is the instruction guide for using the new, repertoire-wide
version of `IgPhyML <https://bitbucket.org/kbhoehn/igphyml>`_. A
description of HLP17, the model underlying IgPhyML, can be found
here:

    Hoehn, K. B., Lunter, G., & Pybus, O. G. (2017). A phylogenetic codon
    substitution model for antibody lineages. Genetics, 206(1), 417-427.
    doi:http://dx.doi.org/10.1534/genetics.116.196303

.. warning::

    The new repertoire-wide version of IgPhyML isn't officially released yet,
    so please let us know if you’re planning to publish anything using it.
    A lot of the features are still in active development, so feedback,
    especially with input/output format and problems encountered, would be
    greatly appreciated.
 
    Also, if you’re analysing single lineages, all of the old features
    still work, so check out the old manual for that
    (see ``docs/IgPhyML_Manual.pdf`` in the IgPhyML source bundle).

Installation
-----------------
 
Download and install the new `IgPhyML <https://bitbucket.org/kbhoehn/igphyml>`_::

    git clone https://bitbucket.org/kbhoehn/igphyml
    cd igphyml
 
Linux
~~~~~~~~~

On Linux operating systems, you can usually just run::

    ./make_phyml_blas_omp

If these commands dont work, you may need to install BLAS and LAPACK,
which provide libraries for better doing matrix exponentiation
operations. In Ubuntu Linux, these are provided in the packages
``libblas-dev`` and ``liblapack-dev``. Other distros probably have
similar package names. If you still can’t get this to work, IgPhyML
can be compiled with OpenMP but without BLAS, though this may negatively
affect performance::
 
    ./make_phyml_omp
 
or without OpenMP as well (though this will slow analysis
considerably)::
 
    ./make_phyml

Once everything is compiled, just add the igphyml/src directory to your
``PATH`` variable, and IgPhyML should work. Importantly, some directory
files are hardcoded in at compilation, so re-compile IgPhyML if you move
the installation directory. Alternatively, you can set the IGPHYML_PATH
environment variable to the location of the igphyml/src/motifs folder for
the same effect.

Mac OS X
~~~~~~~~~~

Installation on Mac OS X is trickier, but possible. The primary issue
is gaining OpenMP support, and installing some GNU command line tools.
The best way is to just install the latest version of ``llvm``
available through ``homebrew``, as well as ``autoconf`` and
``automake``. To do these you’ll need to install
`homebrew <http://brew.sh/index.html>`_. If it’s already installed be
sure it’s at the latest version (``brew update``). You may need to install
Xcode as well. Next, install ``autoconf``, ``automake``, and ``llvm``::

    brew install autoconf
    brew install automake
    brew install llvm

Specify the ``llvm`` version of ``clang`` in ``Makefile.am`` and
``src/Makefile.am`` by adding the line ``CC=<path to llvm clang>``
to the beginning of both files. You will also need to add
``MACOMP=<path to omp.h>`` and ``MACLLVM=<path to llvm lib>`` to
``src/Makefile.am``. For instance, if you’ve install ``llvm 3.9.1``
via homebrew, you will likely need to add the line
``CC=/usr/local/Cellar/llvm/3.9.1/bin/clang``
to ``Makefile.am`` and the lines::

    CC=/usr/local/Cellar/llvm/3.9.1/bin/clang
    MACOMP=/usr/local/Cellar/llvm/3.9.1/lib/clang/3.9.1/include/omp.h
    MACLLVM=/usr/local/Cellar/llvm/3.9.1/lib

to ``src/Makefile.am``.
Your specific path may look different, but you can check locations
of these files and folders by looking around in
``/usr/local/Cellar/llvm/``. The directory structure should be
similar. Run ``./make_blas_phyml_omp``, or other versions, as desired, and add
the ``src`` folder to your ``PATH`` variable.

On some versions of OS X it may be necessary to install XCode command
line tools using::

    xcode-select --install
    cd /Library/Developer/CommandLineTools/Packages/
    open macOS_SDK_headers_for_macOS_<OS X version>.pkg


Quick start
-------------------------------------------------------------------------------

These commands should work as a first pass on many reasonably sized
datasets, but if you really want to understand what’s going on or make
sure what you’re doing makes sense, please check out the rest of the
manual.
 
**Convert Change-O files into IgPhyML inputs**
 
Move to the ``examples`` subfolder and run::

    BuildTrees.py -d example.tab --outname ex --log ex.log --collapse
 
**Build lineage trees using the GY94 model** (fast, doesn’t correct
for hotspots)::
 
    igphyml --repfile ex_lineages.tsv -m GY
 
**Build lineage trees using the HLP17 model** (slower, corrects for
WRC/GYW hotspots)::
 
    igphyml --repfile ex_lineages.tsv -m HLP17
 
Both of these can be parallelized by adding
``--threads <thread count>`` option. Trees files are listed as
``ex/<clone id>.fa_igphyml_tree.txt``, and can be viewed with most
tree viewers (I recommend
`FigTree <http://tree.bio.ed.ac.uk/software/figtree/>`__). Parameter
estimates are in ``ex_lineages.tsv_igphyml_stats.txt``.


Processing Change-O data sets
-------------------------------------------------------------------------------

The process begins with a Change-O formatted :ref:`data file <Standard>`, in
which each sequence has been :ref:`clustered <Cloning>` into a clonal group,
which has subsequently had its unmutated V and J sequence :ref:`predicted <Germlines>`.
 
Use :ref:`BuildTrees` to break this file into separate sequence
alignment files that can be used with IgPhyML. This program will:

1. Filter out nonfunctional sequences.
2. Mask codons split by insertions.
3. Separate clonal groups into separate alignment files (aligned by IMGT site) and information files
4. Create the repertoire files for this dataset.
 
Create IgPhyML input files from ``examples/example.tab``::
 
    cd examples
    BuildTrees.py -d example.tab --outname ex --log ex.log --collapse
 
This will create the directory ``ex`` and the file
``ex_lineages.tsv``. Each ``ex/<clone ID>.fa`` contains the IMGT
mutliple sequence alignemt for a particular clone, and each
``ex/<clone ID>.part.txt`` file contains information about V and J
germline assignments, as well as IMGT unique numbering for each site.
The file ``ex.log`` will contain information about whether or not each
sequence was included in the analysis. The file ``ex_lineages.tsv`` is
the direct input to IgPhyML. Each line represents a clone and shows
the multiple sequence alignment, starting tree topology (N if
ignored), germline sequence ID in alignment file, and partition file
(N if ignored). These repertoire files start with the number of
lineages in the repertoire, and lineages are arranged from most to
least number of sequences. Here, the ``--collapse`` flag is used to
collapse identical sequences. This is highly recommended because
identical sequences slow down calculations without actually affecting
likelihood values in IgPhyML.
 
.. note::

    IgPhyML requires at least three sequences in a lineage, so in
    the case that there is only one observed sequence within a clone, that
    sequence is duplicated. This will not affect the likelihood
    caluclation because these seqeunces will have a branch length of zero,
    but it will affect metrics that take sequence frequency into account.
    You can find further explanation of the different options in the
    :ref:`commandline help <BuildTrees>`,
    including controlling output directories and file names.

IgPhyML Analysis
-------------------------------------------------------------------------------

IgPhyML analysis consists of estimating maximum likelihood (ML) tree
topologies and substitution model parameters for a set of clonal
sequence alignments.

The HLP17 model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The HLP17 model is the heart of IgPhyML and adjusts for features of
affinity maturation that violate the assumptions of most other
phylogenetic models. It uses four sets of parameters to characterize
the types of mutations the occurred over a lineage’s development, and
to help build the tree.
 
:math:`\omega`: Also called dN/dS, or the ratio of nonsynonymous
(amino acid replacement) and synonymous (silent) mutation rates. This
parameter generally relates to clonal selection, with totally neutral
amino acid evolution having an :math:`\omega \approx 1`, negative
selection indicated by :math:`\omega < 1` and diversifying selection
indicated by :math:`\omega > 1`. Generally, when using a partitioned
model (see "Partition models"), we find a lower :math:`\omega` for FWRs than
CDRs, presumably because FWRs are more structurally constrained.
 
:math:`\kappa`: Ratio of transitions (within purines/pyrimidines) to
transversions (between purines/pyrimidines). For normal somatic
hypermutation this ratio is usually :math:`\approx 2`.
 
Motif mutability (e.g. :math:`h^{WRC}`): Mutability parameters for
specified hot- and coldspot motifs. These estimates are equivalent to
the fold-change in mutability for that motif compared to regular
motifs, minus one. So, :math:`h^{WRC} > 0` indicates at hotspot,
:math:`h^{WRC} < 0` indicates a coldspot, and :math:`h^{WRC} = 2`
indicates a 3x increase in *WRC* substitution rate.
 
Codon frequencies (:math:`\pi`): These are calculated using separate
estimates for each nucleotide at each of the three codon positions,
and so are estimated using twelve nucleotide frequency parameters.
These don’t have an immediate interpretation, but are estimated for
each dataset by ML unless fixed to empirical estimates using
``-f empirical``.

Building B cell lineage trees
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before doing any further analysis, I strongly recommend estimating
intitial tree topologies using the GY94 model. This can improve
runtime for HLP17 analysis::
 
    igphyml --repfile ex_lineages.tsv -m GY --outrep ex_lineages.GY.tsv --run_id GY
 
Here, the data files are specifed with ``--repfile``. Topologies are
searched using NNI moves. To do a more thorough topology search, use
``-s SPR``. The flag ``--outrep`` will create a repertoire file that is
identical to the file specified in ``--repfile`` but with the resulting
GY94 topologies specified for each lineage. We can view the ML
parameter estimates for the GY94 fit in
``ex_lineages.tsv_igphyml_stats_GY.txt``, and the tree topologies for
each clone individual lineage in
``ex/<clone id>.fa_igphyml_tree_GY.txt``. I recommend using
`FigTree <http://tree.bio.ed.ac.uk/software/figtree/>`__ to visualize
topologies.
 
To estimate ML tree topologies using the HLP17 model wth a GY94
starting topology, use::
 
    igphyml --repfile ex_lineages.GY.tsv -m HLP17 --run_id HLP --threads 2
 
This will estimate a single :math:`\omega`, :math:`\kappa`, set of
codon frequencies (:math:`\pi`), and WRC/GYW mutability across the
entire repertoire, and search for topologies using NNI moves. You can
see parameter estimates in
``ex_lineages.GY.tsv_igphyml_stats_HLP.txt``, and trees in
``ex/<clone id>.fa_igphyml_tree_HLP.txt``. This command will also
parallelize the calculation across 2 threads using the ``--threads``
flag.

Heirarchical substitution models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Substitution models are specified using the ``-t`` for :math:`\kappa`
(transition/transverion rate), ``--omegaOpt`` for :math:`\omega`
(nonsynonymous/synonymous mutation rate), and ``--motifs`` and
``--hotness`` for specifying the motif mutability models. The default
for all of these is to estimate shared parameter values across all
lineages, which is also specified by ``e``. The default motif model is
symmetric WRC/GYW. So, the following two commands are equivalent::
 
    igphyml --repfile ex_lineages.GY.tsv -m HLP17 -o lr --run_id HLP
 
    igphyml --repfile ex_lineages.GY.tsv -m HLP17 -t e --omegaOpt e,e --motifs WRC_2:0,GYW_0:1 \
        --hotness e,e -o lr --run_id HLP
 
In both cases parameter estimates are recorded in
``ex_lineages.GY.tsv_igphyml_stats_HLP.txt``. Note that here we use
``-o lr``, which will only optimize branch lengths and substitution
parameters. This will keep topologies the same as the GY94, but will
estimate substitution parameters much more quickly. To estimate
mutabilities of all six canonical hotspot motifs, use ``--motifs FCH``,
for ‘Free coldspots and hotspots’, though this will result in extreme
parameter values if there is insufficient information in the
repertoire file.

Partition models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To estimate separate values of :math:`\omega` for CDR/FWR partitions,
specify more than one value in the ``--omegaOpt`` option. For instance::
 
    igphyml --repfile ex_lineages.GY.tsv -m HLP17 --omegaOpt e,e -o lr --run_id HLP
 
Will estimate a separate :math:`\omega` at the repertoire level for
the FWRs (‘Omega 0’) and CDRs (‘Omega 1’) of each lineage. This is the default behavior
if partition files are specified. If partition files are specified and you only
want a single :math:`\omega` use ``--omegaOpt e``.


Optimizing performance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

IgPhyML is a computationally intensive program. There are some ways to
make calculations more practical, however:
 
GY94 starting topologies: Calculations are much faster under the GY94
model (see [top]), so it is usually better to do an initial topolgoy
searching under the GY94 model, and then using those trees as starting
topologies for HLP17 . You can also fix these topologies during HLP17
parameter estimation (``-o lr``) for an even greater speedup, though,
obviously, this will not result in a change in topology from GY94.
 
Enforcing minimum lineage size: Many repertoires often contain huge
numbers of small lineages that can make computations impractical. To
limit the size of lineages being analyzed, specify a cutoff with
``--minSeq``, and note that 1) the germline sequence is added to
sequence files, and 2) single sequence lineages are duplicated (see
"Processing Change-O data sets") and thus have three sequences total. So, to limit analyses to
lineages with at least three observed sequences, use ``--minSeq 4``.
``--minSeq 3`` and ``--minSeq 2`` are identical because single lineages
have duplicated sequences, and ``--minSeq 1`` is useless.
 
Parallelizing computations: It is possible to parallelize likelihood
calulcations using the ``--threads`` option. By default, calculations
are parallelized by tree, so there is no point in using more threads
than you have lineages in your repertoire file. If you are analyzing a
single large lineage, or a repertoire dominated by one lineage and a
couple of much smaller lineages, it may be more efficient to instead
parallelize by site. To do this, add ``--splitByTree 0`` to parallelize
calculations within each tree, and analyze the trees sequentially.
