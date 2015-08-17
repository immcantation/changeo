Version 0.3.0.beta:  August 18, 2015
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
  
  
Version 0.2.4:  August 18, 2015
-------------------------------------------------------------------------------

MakeDb:

+ Fixed a bug wherein N1 and N2 region indexing was off by one nucleotide
  for the igblast subcommand (leading to incorrect SEQUENCE_VDJ values).

ParseDb:

+ Fixed a bug wherein specifying the `-f` argument to the index subcommand 
  would cause an error.
  

Version 0.2.3:  July 22, 2015
-------------------------------------------------------------------------------

DefineClones:

+ Fixed a typo in the default normalization setting of the bygroup subcommand, 
  which was being interpreted as 'none' rather than 'len'.
+ Changed the 'hs5f' model of the bygroup subcommand to be centered -log10 of 
  the targeting probability.
+ Added the `--sym` argument to the bygroup subcommand which determines how 
  asymmetric distances are handled.
   

Version 0.2.2:  July 8, 2015
-------------------------------------------------------------------------------

CreateGermlines:

+ Germline creation now works for IgBLAST output parsed with MakeDb. The 
  argument `--sf SEQUENCE_VDJ` must be provided to generate germlines from 
  IgBLAST output. The same reference database used for the IgBLAST alignment
  must be specified with the `-r` flag.
+ Fixed a bug with determination of N1 and N2 region positions.

MakeDb:

+ Combined the `-z` and `-f` flags of the imgt subcommand into a single flag, 
  `-i`, which autodetects the input type.
+ Added requirement that IgBLAST input be generated using the 
  `-outfmt "7 std qseq"` argument to igblastn.
+ Modified SEQUENCE_VDJ output from IgBLAST parser to include gaps inserted 
  during alignment.
+ Added correction for IgBLAST alignments where V/D, D/J or V/J segments are
  assigned overlapping positions.
+ Corrected N1_LENGTH and N2_LENGTH calculation from IgBLAST output.
+ Added the `--scores` flag which adds extra columns containing alignment 
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
+ Added `--failed` flag to create database with records that fail analysis.
+ Added `--sf` flag to specify sequence field to be analyzed.

CreateGermlines:

+ Fixed a bug where germline sequences could not be created for light chains.

DefineClones:

+ Added a human 1-mer model, 'hs1f', which uses the substitution rates from 
  from Yaari et al, 2013.
+ Changed default model to 'hs1f' and default normalization to length for 
  bygroup subcommand.
+ Added `--link` argument which allows for specification of single, complete,
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