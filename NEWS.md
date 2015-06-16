Version 0.2.0:  June 16, 2015
-------------------------------------------------------------------------------

Initial public prerelease.  
Output files were added to the usage documentation of all scripts.  
General code cleanup.  

AnalyzeAa:

+ Fixed a bug where junctions less than one codon long would lead to a 
  division by zero error.

DefineClones:

+ Added a human 1-mer model, `hs1f`, which uses the substitution rates from 
  from Yaari et al, 2013
  
+ Added `--link` argument which allows for specification of single, complete,
  or average linkage during clonal clustering.

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