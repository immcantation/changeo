#!/bin/bash
# Script to run MakeDb and ParseDb on multiple inputs
# 
# Author:  Jason Anthony Vander Heiden
# Date:    2015.03.09
#
# Required Arguments:
#   $1 = A two column tab delimited file mapping IMGT zip files (column 1) 
#        to submitted FASTA files (column 2)
#   $2 = Output directory

while read FILE_MAP
do
    FILE_ARRAY=($FILE_MAP)
    BASENAME=$(basename "${FILE_ARRAY[0]}" ".zip")
    MakeDb.py imgt -z "${FILE_ARRAY[0]}" -s "${FILE_ARRAY[1]}" --outdir $2
    ParseDb.py select -d "$2/${BASENAME}_db-pass.tab" -f FUNCTIONAL -u T
    mv "$2/${BASENAME}_db-pass_parse-select.tab" "$2/${BASENAME}_functional.tab"
done < $1

