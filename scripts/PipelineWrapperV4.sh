#!/bin/bash
# Wrapper script to run the Change-O pipeline script on multiple inputs
# 
# Author:  Jason Anthony Vander Heiden
# Date:    2014.8.24
#
# Required Arguments:
#   $1 = A two column tab delimited file mapping IMGT zip files (column 1) 
#        to submitted FASTA files (column 2)
#   $2 = Output directory

SCRIPT=/home/jason/apps/changeo-0.4/scripts/ClonePipelineV4_IMGT.sh
LOGFILE=ChangeoRunLog.out
NPROC=4

while read FILE_MAP
do
    FILE_ARRAY=($FILE_MAP)
    FOLDER=$(basename ${FILE_ARRAY[0]} | tr -d ".zip")
    
    echo "FOLDER: $FOLDER" | tee -a $LOGFILE 
    echo `date` | tee -a $LOGFILE
    $SCRIPT ${FILE_ARRAY[0]} ${FILE_ARRAY[1]} $2/$FOLDER $NPROC | tee -a $LOGFILE
done < $1