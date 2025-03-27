Release Notes
===============================================================================

Version 1.3.1:  March 27, 2025
-------------------------------------------------------------------------------

+ Active development has moved from Bitbucket to GitHub (https://github.com/immcantation/changeo).
+ Documentation updates.
+ Various updates to internals to avoid deprecation warnings.
+ Updates in requirements: Python 3.7.0, biopython>=1.81 and packaging>=23.2.

Version 1.3.0:  December 11, 2022
-------------------------------------------------------------------------------

+ Various updates to internals and error messages.

AssignGenes:

+ Added support for ``.fastq`` files. If a ``.fastq`` file is input, then a
  corresponding ``.fasta`` file will be created in output directory.
+ Added support for C region alignment calls provide by IgBLAST v1.18+.

MakeDb:

+ Added support for C region alignment calls provide by IgBLAST v1.18+.


Version 1.2.0:  October 29, 2021
-------------------------------------------------------------------------------

+ Updated dependencies to presto >= v0.7.0.

AssignGenes:

+ Fixed reporting of IgBLAST output counts when specifying ``--format airr``.

BuildTrees:

+ Added support for specifying fixed omega and hotness parameters at the
  commandline.

CreateGermlines:

+ Will now use the first allele in the reference database when duplicate
  allele names are provided. Only appears to affect mouse BCR light chains
  and TCR alleles in the IMGT database when the same allele name differs by
  strain.

MakeDb:

+ Added support for changes in how IMGT/HighV-QUEST v1.8.4 handles special
  characters in sequence identifiers.
+ Fixed the ``imgt`` subcommand incorrectly allowing execution without
  specifying the IMGT/HighV-QUEST output file at the commandline.

ParseDb:

+ Added reporting of output file sizes to the console log of the ``split``
  subcommand.


Version 1.1.0:  June 21, 2021
-------------------------------------------------------------------------------

+ Fixed gene parsing for IMGT temporary designation nomenclature.
+ Updated dependencies to biopython >= v1.77, airr >= v1.3.1, PyYAML>=5.1.

MakeDb:

+ Added the ``--imgt-id-len`` argument to accommodate changes introduced in how
  IMGT/HighV-QUEST truncates sequence identifiers as of v1.8.3 (May 7, 2021).
  The header lines in the fasta files are now truncated to 49 characters. In
  IMGT/HighV-QUEST versions older than v1.8.3, they were truncated to 50 characters.
  ``--imgt-id-len`` default value is 49. Users should specify ``--imgt-id-len 50``
  to analyze IMGT results generated with IMGT/HighV-QUEST versions older than v1.8.3.
+ Added the ``--infer-junction`` argument to ``MakeDb igblast``, to enable the inference
  of the junction sequence when not reported by IgBLAST. Should be used with data from
  IgBLAST v1.6.0 or older; before igblast added the IMGT-CDR3 inference.


Version 1.0.2:  January 18, 2021
-------------------------------------------------------------------------------

AlignRecords:

+ Fixed a bug caused the program to exit when encountering missing sequence
  data. It will now fail the row or group with missing data and continue.

MakeDb:

+ Added support for IgBLAST v1.17.0.

ParseDb:

+ Added a relevant error message when an input field is missing from the data.


Version 1.0.1:  October 13, 2020
-------------------------------------------------------------------------------

+ Updated to support Biopython v1.78.
+ Increased the biopython dependency to v1.71.
+ Increased the presto dependency to 0.6.2.


Version 1.0.0:  May 6, 2020
-------------------------------------------------------------------------------

+ The default output in all tools is now the AIRR Rearrangement standard
  (``--format airr``). Support for the legacy Change-O data standard is still
  provided through the ``--format changeo`` argument to the tools.
+ License changed to AGPL-3.

AssignGenes:

+ Added the ``igblast-aa`` subcommand to run igblastp on amino acid input.

BuildTrees:

+ Adjusted ``RECORDS`` to indicate all sequences in input file.
  ``INITIAL_FILTER`` now shows sequence count after initial
  ``min_seq`` filtering.
+ Added option to skip codon masking: ``--nmask``.
+ Mask ``:``, ``,``, ``)``, and ``(`` in IDs and metadata with ``-``.
+ Can obtain germline from ``GERMLINE_IMGT`` if ``GERMLINE_IMGT_D_MASK``
  not specified.
+ Can reconstruct intermediate sequences with IgPhyML using ``--asr``.

ConvertDb:

+ Fixed a bug in the ``airr`` subcommand that caused the ``junction_length``
  field to be deleted from the output.
+ Fixed a bug in the ``genbank`` subcommand that caused the junction CDS
  to be missing from the ASN output.

CreateGermlines:

+ Added the ``--cf`` argument to allow specification of the clone field.

MakeDb:

+ Added the ``igblast-aa`` subcommand to parse the output of igblastp.
+ Changed the log entry ``FUNCTIONAL`` to ``PRODUCTIVE`` and removed the
  ``IMGT_PASS`` log entry in favor of an informative ``ERROR`` entry
  when sequences fail the junction region validation.
+ Add --regions argument to the ``igblast`` and ``igblast-aa`` subcommands
  to allow specification of the IMGT CDR/FWR region boundaries. Currently,
  the supported specifications are ``default`` (human, mouse) and
  ``rhesus-igl``.


Version 0.4.6:  July 19, 2019
-------------------------------------------------------------------------------

BuildTrees:

+ Added capability of running IgPhyML on outputted data (``--igphyml``) and
  support for passing IgPhyML arguments through BuildTrees.
+ Added the ``--clean`` argument to force deletion of all intermediate files
  after IgPhyML execution.
+ Added the ``--format`` argument to allow specification input and output of
  either the Change-O standard (``changeo``) or AIRR Rearrangement standard
  (``airr``).

CreateGermlines:

+ Fixed a bug causing incorrect reporting of the germline format in the
  console log.

ConvertDb:

+ Removed requirement for the ``NP1_LENGTH`` and ``NP2_LENGTH`` fields from
  the genbank subcommand.

DefineClones:

+ Fixed a biopython warning arising when applying ``--model aa`` to junction
  sequences that are not a multiple of three. The junction will now be
  padded with an appropriate number of Ns (usually resulting in a translation
  to X).

MakeDb:

+ Added the ``--10x`` argument to all subcommands to support merging of
  Cell Ranger annotation data, such as UMI count and C-region assignment,
  with the output of the supported alignment tools.
+ Added inference of the receptor locus from the alignment data to all
  subcommands, which is output in the ``LOCUS`` field.
+ Combined the extended field arguments of all subcommands (``--scores``,
  ``--regions``, ``--cdr3``, and ``--junction``) into a single ``--extended``
  argument.
+ Removed parsing of old IgBLAST v1.5 CDR3 fields
  (``CDR3_IGBLAST``, ``CDR3_IGBLAST_AA``).


Version 0.4.5:  January 9, 2019
-------------------------------------------------------------------------------

+ Slightly changed version number display in commandline help.

BuildTrees:

+ Fixed a bug that caused malformed lineages.tsv output file.

CreateGermlines:

+ Fixed a bug in the CreateGermlines log output causing incorrect missing
  D gene or J gene error messages.

DefineClones:

+ Fixed a bug that caused a missing junction column to cluster sequences
  together.

MakeDb:

+ Fixed a bug that caused failed germline reconstructions to be recorded as
  ``None``, rather than an empty string, in the ``GERMLINE_IMGT`` column.


Version 0.4.4:  October 27, 2018
-------------------------------------------------------------------------------

+ Fixed a bug causing the values of ``_start`` fields to be off by one from
  the v1.2 AIRR Schema requirement when specifying ``--format airr``.


Version 0.4.3:  October 19, 2018
-------------------------------------------------------------------------------

+ Updated airr library requirement to v1.2.1 to fix empty V(D)J start
  coordinate values when specifying ``--format airr`` to tools.
+ Changed pRESTO dependency to v0.5.10.

BuildTrees:

+ New tool.
+ Converts tab-delimited database files into input for
  `IgPhyML <https://bitbucket.org/kbhoehn/igphyml>`_

CreateGermlines:

+ Now verifies that all files/folder passed to the ``-r`` argument exist.


Version 0.4.2:  September 6, 2018
-------------------------------------------------------------------------------

+ Updated support for the AIRR Rearrangement schema to v1.2 and added the
  associated airr library dependency.

AssignGenes:

+ New tool.
+ Provides a simple IgBLAST wrapper as the ``igblast`` subcommand.

ConvertDb:

+ The ``genbank`` subcommand will perform a check for some of the required
  columns in the input file and exit if they are not found.
+ Changed the behavior of the ``-y`` argument in the ``genbank`` subcommand.
  This argument is now featured to sample features only, but allows
  for the inclusion of any BioSample attribute.

CreateGermlines:

+ Will now perform a naive verification that the reference sequences provided
  to the ``-r`` argument are IMGT-gapped. A warning will be issued to standard
  error if the reference sequence fail the check.
+ Will perform a check for some of the required columns in the input file and
  exit if they are not found.

MakeDb:

+ Changed the output of ``SEQUENCE_VDJ`` from the igblast subcommand to retain
  insertions in the query sequence rather than delete them as is done in the
  ``SEQUENCE_IMGT`` field.
+ Will now perform a naive verification that the reference sequences provided
  to the ``-r`` argument are IMGT-gapped. A warning will be issued to standard
  error if the reference sequence fail the check.


Version 0.4.1:  July 16, 2018
-------------------------------------------------------------------------------

+ Fixed installation incompatibility with pip 10.
+ Fixed duplicate newline issue on Windows.
+ All tools will no longer create empty pass or fail files if there are no
  records meeting the appropriate criteria for output.
+ Most tools now allow explicit specification of the output file name via
  the optional ``-o`` argument.
+ Added support for the AIRR standard TSV via the ``--format airr`` argument to
  all relevant tools.
+ Replaced V, D and J ``BTOP`` columns with ``CIGAR`` columns in data standard.
+ Numerous API changes and internal structural changes to commandline tools.

AlignRecords:

+ Fixed a bug arising when space characters are present in the sequence
  identifiers.

ConvertDb:

+ New tool.
+ Includes the airr and changeo subcommand to convert between AIRR and Change-O
  formatted TSV files.
+ The genbank subcommand creates MiAIRR compliant files for submission to
  GenBank/TLS.
+ Contains the baseline and fasta subcommands previously in ParseDb.

CreateGermlines

+ Changed character used to pad clonal consensus sequences from ``.`` to ``N``.
+ Changed tie resolution in clonal consensus from random V/J gene to
  alphabetical by sequence identifier.
+ Added ``--df`` and ``-jf`` arguments for specifying D and J fields,
  respectively.
+ Add initial sorting step with specifying ``--cloned`` so that clonally
  ordered input is no longer required.

DefineClones:

+ Removed the chen2010 and ademokun2011 and made the previous bygroup
  subcommand the default behavior.
+ Renamed the ``--f`` argument to ``--gf`` for consistency with other tools.
+ Added the arguments ``--vf`` and ``-jf`` to allow specification of
  V and J call fields, respectively.

MakeDb:

+ Renamed ``--noparse`` argument to ``--asis-id``.
+ Added ``asis-calls`` argument to igblast subcommand to allow use with
  non-standard gene names.
+ Added the ``GERMLINE_IMGT`` column to the default output.
+ Changed junction inference in igblast subcommand to use IgBLAST's CDR3
  assignment for IgBLAST versions greater than or equal to 1.7.0.
+ Added a verification that the ``SEQUENCE_IMGT`` and ``JUNCTION`` fields
  are in agreement for records to pass.
+ Changed behavior of the igblast subcommand's translation of the junction
  sequence to truncate junction that are not multiples of 3, rather than
  pad to a multiple of 3 (removes trailing X character).
+ The igblast subcommand will now fail records missing the required optional
  fields ``subject seq``, ``query seq`` and ``BTOP``, rather than abort.
+ Fixed bug causing parsing of IgBLAST <= 1.4 output to fail.

ParseDb:

+ Added the merge subcommand which will combine TSV files.
+ All field arguments are now case sensitive to provide support for both
  the Change-O and AIRR data standards.


Version 0.3.12:  February 16, 2018
-------------------------------------------------------------------------------

MakeDb:

+ Fixed a bug wherein specifying multiple simultaneous inputs would cause
  duplication of parsed pRESTO fields to appear in the second and higher
  output files.


Version 0.3.11:  February 6, 2018
-------------------------------------------------------------------------------

MakeDb:

+ Fixed junction inferrence for igblast subcommand when J region is
  truncated.


Version 0.3.10:  February 6, 2018
-------------------------------------------------------------------------------

Fixed incorrect progress bars resulting from files containing empty lines.

DefineClones:

+ Fixed several bugs in the chen2010 and ademokun2011 methods that caused them
  to either fail or incorrectly cluster all sequences into a single clone.
+ Added informative message for out of memory error in chen2010 and
  ademokun2011 methods.


Version 0.3.9:  October 17, 2017
-------------------------------------------------------------------------------

DefineClones:

+ Fixed a bug causing DefineClones to fail when all are sequences removed from
  a group due to missing characters.


Version 0.3.8:  October 5, 2017
-------------------------------------------------------------------------------

AlignRecords:

+ Ressurrected AlignRecords which performs multiple alignment of sequence
  fields.
+ Added new subcommands ``across`` (multiple aligns within columns),
  ``within`` (multiple aligns columns within each row), and ``block``
  (multiple aligns across both columns and rows).

CreateGermlines:

+ Fixed a bug causing CreateGermlines to incorrectly fail records when using
  the argument ``--vf V_CALL_GENOTYPED``.

DefineClones:

+ Added the ``--maxmiss`` argument to the bygroup subcommand of DefineClones
  which set exclusion criteria for junction sequence with ambiguous and
  missing characters. By default, bygroup will now fail all sequences
  with any missing characters in the junction (``--maxmiss 0``).


Version 0.3.7:  June 30, 2017
-------------------------------------------------------------------------------

MakeDb:

+ Fixed an incompatibility with IgBLAST v1.7.0.

CreateGermlines:

+ Fixed an error that occurs when using the ``--cloned`` with an input file
  containing duplicate values in ``SEQUENCE_ID`` that caused some records to
  be discarded.


Version 0.3.6:  June 13, 2017
-------------------------------------------------------------------------------

+ Fixed an overflow error on Windows that caused tools to fatally exit.
+ All tools will now print detailed help if no arguments are provided.


Version 0.3.5:  May 12, 2017
-------------------------------------------------------------------------------

Fixed a bug wherein ``.tsv`` was not being recognized as a valid extension.

MakeDb:

+ Added the ``--cdr3`` argument to the igblast subcommand to extract the
  CDR3 nucleotide and amino acid sequence defined by IgBLAST.
+ Updated the IMGT/HighV-QUEST parser to handle recent column name changes.
+ Fixed a bug in the igblast parser wherein some sequence identifiers were
  not being processed correctly.

DefineClones:

+ Changed the way ``X`` characters are handled in the amino acid Hamming
  distance model to count as a match against any character.


Version 0.3.4:  February 14, 2017
-------------------------------------------------------------------------------

License changed to Creative Commons Attribution-ShareAlike 4.0 International
(CC BY-SA 4.0).

CreateGermlines:

+ Added ``GERMLINE_V_CALL``, ``GERMLINE_D_CALL`` and ``GERMLINE_J_CALL``
  columns to the output when the ``-cloned`` argument is specified. These
  columns contain the consensus annotations when clonal groups contain
  ambiguous gene assignments.
+ Fixed the error message for an invalid repo (``-r``) argument.

DefineClones:

+ Deprecated ``m1n`` and ``hs1f`` distance models, renamed them to
  ``m1n_compat`` and ``hs1f_compat``, and replaced them with ``hh_s1f`` and
  replaced ``mk_rs1nf``, respectively.
+ Renamed the ``hs5f`` distance model to ``hh_s5f``.
+ Added the mouse specific distance model ``mk_rs5nf`` from Cui et al, 2016.

MakeDb:

+ Added compatibility for IgBLAST v1.6.
+ Added the flag ``--partial`` which tells MakeDb to pass incomplete alignment
  results specified.
+ Added missing console log entries for the ihmm subcommand.
+ IMGT/HighV-QUEST, IgBLAST and iHMMune-Align parsers have been cleaned up,
  better documented and moved into the iterable classes
  ``changeo.Parsers.IMGTReader``, ``change.Parsers.IgBLASTReader``, and
  ``change.Parsers.IHMMuneReader``, respectively.
+ Corrected behavior of ``D_FRAME`` annotation from the ``--junction``
  argument to the imgt subcommand such that it now reports no value when no
  value is reported by IMGT, rather than reporting the reading frame as 0 in
  these cases.
+ Fixed parsing of ``IN_FRAME``, ``STOP``, ``D_SEQ_START`` and ``D_SEQ_LENGTH``
  fields from iHMMune-Align output.
+ Removed extraneous score fields from each parser.
+ Fixed the error message for an invalid repo (``-r``) argument.


Version 0.3.3:  August 8, 2016
-------------------------------------------------------------------------------

Increased ``csv.field_size_limit`` in changeo.IO, ParseDb and DefineClones
to be able to handle files with larger number of UMIs in one field.

Renamed the fields ``N1_LENGTH`` to ``NP1_LENGTH`` and ``N2_LENGTH``
to ``NP2_LENGTH``.

CreateGermlines:

+ Added differentiation of the N and P regions the the ``REGION`` log field
  if the N/P region info is present in the input file (eg, from the
  ``--junction`` argument to MakeDb-imgt). If the additional N/P region
  columns are not present, then both N and P regions will be denoted by N,
  as in previous versions.
+ Added the option 'regions' to the ``-g`` argument to create add the
  ``GERMLINE_REGIONS`` field to the output which represents the germline
  positions as V, D, J, N and P characters. This is equivalent to the
  ``REGION`` log entry.

DefineClones:

+ Improved peformance significantly of the ``--act set`` grouping method in
  the bygroup subcommand.

MakeDb:

+ Fixed a bug producing ``D_SEQ_START`` and ``J_SEQ_START`` relative to
  ``SEQUENCE_VDJ`` when they should be relative to ``SEQUENCE_INPUT``.
+ Added the argument ``--junction`` to the imgt subcommand to parse additional
  junction information fields, including N/P region lengths and the D-segment
  reading frame. This provides the following additional output fields:
  ``D_FRAME``, ``N1_LENGTH``, ``N2_LENGTH``, ``P3V_LENGTH``, ``P5D_LENGTH``,
  ``P3D_LENGTH``, ``P5J_LENGTH``.
+ The fields ``N1_LENGTH`` and ``N2_LENGTH`` have been renamed to accommodate
  adding additional output from IMGT under the ``--junction`` flag. The new
  names are ``NP1_LENGTH`` and ``NP2_LENGTH``.
+ Fixed a bug that caused the ``IN_FRAME``, ``MUTATED_INVARIANT`` and
  ``STOP`` field to be be parsed incorrectly from IMGT data.
+ Ouput from iHMMuneAlign can now be parsed via the ``ihmm`` subcommand.
  Note, there is insufficient information returned by iHMMuneAlign to
  reliably reconstruct germline sequences from the output using
  CreateGermlines.


ParseDb:

+ Renamed the clip subcommand to baseline.


Version 0.3.2:  March 8, 2016
-------------------------------------------------------------------------------

Fixed a bug with installation on Windows due to old file paths lingering in
changeo.egg-info/SOURCES.txt.

Updated license from CC BY-NC-SA 3.0 to CC BY-NC-SA 4.0.

CreateGermlines:

+ Fixed a bug producing incorrect values in the ``SEQUENCE`` field on the
  log file.

MakeDb:

+ Updated igblast subcommand to correctly parse records with indels. Now
  igblast must be run with the argument ``outfmt "7 std qseq sseq btop"``.
+ Changed the names of the FWR and CDR output columns added with
  ``--regions`` to ``<region>_IMGT``.
+ Added ``V_BTOP`` and ``J_BTOP`` output when the ``--scores`` flag is
  specified to the igblast subcommand.


Version 0.3.1:  December 18, 2015
-------------------------------------------------------------------------------

MakeDb:

+ Fixed bug wherein the imgt subcommand was not properly recognizing an
  extracted folder as input to the ``-i`` argument.


Version 0.3.0:  December 4, 2015
-------------------------------------------------------------------------------

Conversion to a proper Python package which uses pip and setuptools for
installation.

The package now requires Python 3.4. Python 2.7 is not longer supported.

The required dependency versions have been bumped to numpy 1.9, scipy 0.14,
pandas 0.16 and biopython 1.65.

DbCore:

+ Divided DbCore functionality into the separate modules: Defaults, Distance,
  IO, Multiprocessing and Receptor.

IgCore:

+ Remove IgCore in favor of dependency on pRESTO >= 0.5.0.

AnalyzeAa:

+ This tool was removed. This functionality has been migrated to the alakazam
  R package.

DefineClones:

+ Added ``--sf`` flag to specify sequence field to be used to calculate
  distance between sequences.
+ Fixed bug in wherein sequences with missing data in grouping columns
  were being assigned into a single group and clustered. Sequences with
  missing grouping variables will now be failed.
+ Fixed bug where sequences with "None" junctions were grouped together.

GapRecords:

+ This tool was removed in favor of adding IMGT gapping support to igblast
  subcommand of MakeDb.

MakeDb:

+ Updated IgBLAST parser to create an IMGT gapped sequence and infer the
  junction region as defined by IMGT.
+ Added the ``--regions`` flag which adds extra columns containing FWR and CDR
  regions as defined by IMGT.
+ Added support to imgt subcommand for the new IMGT/HighV-QUEST compression
  scheme (.txz files).


Version 0.2.5:  August 25, 2015
-------------------------------------------------------------------------------

CreateGermlines:

+ Removed default '-r' repository and added informative error messages when
  invalid germline repositories are provided.
+ Updated '-r' flag to take list of folders and/or fasta files with germlines.


Version 0.2.4:  August 19, 2015
-------------------------------------------------------------------------------

MakeDb:

+ Fixed a bug wherein N1 and N2 region indexing was off by one nucleotide
  for the igblast subcommand (leading to incorrect SEQUENCE_VDJ values).

ParseDb:

+ Fixed a bug wherein specifying the ``-f`` argument to the index subcommand
  would cause an error.


Version 0.2.3:  July 22, 2015
-------------------------------------------------------------------------------

DefineClones:

+ Fixed a typo in the default normalization setting of the bygroup subcommand,
  which was being interpreted as 'none' rather than 'len'.
+ Changed the 'hs5f' model of the bygroup subcommand to be centered -log10 of
  the targeting probability.
+ Added the ``--sym`` argument to the bygroup subcommand which determines how
  asymmetric distances are handled.


Version 0.2.2:  July 8, 2015
-------------------------------------------------------------------------------

CreateGermlines:

+ Germline creation now works for IgBLAST output parsed with MakeDb. The
  argument ``--sf SEQUENCE_VDJ`` must be provided to generate germlines from
  IgBLAST output. The same reference database used for the IgBLAST alignment
  must be specified with the ``-r`` flag.
+ Fixed a bug with determination of N1 and N2 region positions.

MakeDb:

+ Combined the ``-z`` and ``-f`` flags of the imgt subcommand into a single flag,
  ``-i``, which autodetects the input type.
+ Added requirement that IgBLAST input be generated using the
  ``-outfmt "7 std qseq"`` argument to igblastn.
+ Modified SEQUENCE_VDJ output from IgBLAST parser to include gaps inserted
  during alignment.
+ Added correction for IgBLAST alignments where V/D, D/J or V/J segments are
  assigned overlapping positions.
+ Corrected N1_LENGTH and N2_LENGTH calculation from IgBLAST output.
+ Added the ``--scores`` flag which adds extra columns containing alignment
  scores from IMGT and IgBLAST output.


Version 0.2.1:  June 18, 2015
-------------------------------------------------------------------------------

DefineClones:

+ Removed mouse 3-mer model, 'm3n'.


Version 0.2.0:  June 17, 2015
-------------------------------------------------------------------------------

Initial public prerelease.

Output files were added to the usage documentation of all scripts.

General code cleanup.

DbCore:

+ Updated loading of database files to convert column names to uppercase.

AnalyzeAa:

+ Fixed a bug where junctions less than one codon long would lead to a
  division by zero error.
+ Added ``--failed`` flag to create database with records that fail analysis.
+ Added ``--sf`` flag to specify sequence field to be analyzed.

CreateGermlines:

+ Fixed a bug where germline sequences could not be created for light chains.

DefineClones:

+ Added a human 1-mer model, 'hs1f', which uses the substitution rates from
  from Yaari et al, 2013.
+ Changed default model to 'hs1f' and default normalization to length for
  bygroup subcommand.
+ Added ``--link`` argument which allows for specification of single, complete,
  or average linkage during clonal clustering (default single).

GapRecords:

+ Fixed a bug wherein non-standard sequence fields could not be aligned.

MakeDb:

+ Fixed bug where the allele 'TRGVA*01' was not recognized as a valid allele.

ParseDb:

+ Added rename subcommand to ParseDb which renames fields.



Version 0.2.0.beta-2015-05-31:  May 31, 2015
-------------------------------------------------------------------------------

Minor changes to a few output file names and log field entries.

ParseDb:

+ Added index subcommand to ParseDb which adds a numeric index field.


Version 0.2.0.beta-2015-05-05:  May 05, 2015
-------------------------------------------------------------------------------

Prerelease for review.
