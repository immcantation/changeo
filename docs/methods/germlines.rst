.. _CreateGermlinesMethod:

Reconstruction of germline sequences from alignment data
================================================================================

The :ref:`CreateGermlines` tool takes the individual segment alignment information for
each sequence and reconstructs a full length germline sequence from the V(D)J 
reference sequences. 

To reconstruct the germline, :ref:`CreateGermlines` trims V(D)J germline segments
and N/P regions by alignment length and concatenates them together. It puts Ns 
in the untemplated N/P regions and optionally masks the D with Ns 
(:option:`-g dmask <CreateGermlines -g>`). :ref:`CreateGermlines` also looks for and 
corrects cases where the alignment tool assigned the same part of the input sequence 
to two different regions (eg, assigning the same nucleotides to N/P and J).

At the end of the germline reconstruction process, each sequence has been assigned
a germline specific to the sequence.

When the (:option:`--cloned <CreateGermlines --cloned>`) flag is specified, the 
process is the same except it is clone specific and results in the 
creation of one germline per clone. :ref:`CreateGermlines` selects first a 
single V, D, and J allele to use as the germline from all the assigned 
annotations in each clone. The selection is made by simple majority rule of all 
the allele calls in the clone. After the germline reconstruction process, all 
sequences belonging to the same clone have been assigned the same germline.
