Release Notes
===============================================================================

Version 0.3.10:  February 6, 2018
-------------------------------------------------------------------------------

Fixed incorrect progress bars resulting from files containing empty lines.

DefineClones:

+ Fixed several bugs in the chen2010 and ademokun2011 methods that caused them
  to eithe fail or incorrectly cluster all sequences into a single clone.
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
