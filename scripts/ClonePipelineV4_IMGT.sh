#!/bin/bash
# Super script to run the Change-O clonal assignment pipeline on IMGT data 
# 
# Author:  Jason Anthony Vander Heiden
# Date:    2014.8.24
# 
# Required Arguments:
#   $1 = IMGT zip file 
#   $2 = FASTA file that was submitted to IMGT
#   $3 = output directory
#   $4 = number of subprocesses for multiprocessing tools

# Define run parameters
IMGT_FILE=$(readlink -f $1)
SEQ_FILE=$(readlink -f $2)
OUTDIR=$(readlink -f $3)
NPROC=$4
LOG_RUNTIMES=false
ZIP_FILES=false
DC_MODEL=h3n
DC_DIST=5
GERMLINE_FOLDER=/mnt/data/germlines/IMGT_Human_2014-08-23
BASENAME=$(basename $IMGT_FILE | tr -d ".zip")

# Define script execution command and log files
mkdir -p $OUTDIR; cd $OUTDIR
RUNLOG="${OUTDIR}/Pipeline.log"
echo '' > $RUNLOG
if $LOG_RUNTIMES; then
    TIMELOG="${OUTDIR}/Time.log"
    echo '' > $TIMELOG 
    RUN="nice -19 /usr/bin/time -o ${TIMELOG} -a -f %C\t%E\t%P\t%Mkb"
else
    RUN="nice -19"
fi
        
# Start
echo "DIRECTORY: ${OUTDIR}"
echo "VERSIONS:"
echo "  $(AAnalysis.py -v 2>&1)"
echo "  $(CreateGermlines.py -v 2>&1)"
echo "  $(DefineClones.py -v 2>&1)"
echo "  $(MakeDb.py -v 2>&1)"
echo "  $(ParseDb.py -v 2>&1)"
echo "  $(SplitDb.py -v 2>&1)"

# Parse IMGT
echo -e "\nSTART"
echo "   1: MakeDb imgt            $(date +'%H:%M %D')"
#$RUN MakeDb.py imgt -z $IMGT_FILE -s $SEQ_FILE --outdir . --log ParseIMGT.log --clean >> $RUNLOG
$RUN MakeDb.py imgt -z $IMGT_FILE -s $SEQ_FILE --outdir . >> $RUNLOG

# Create germlines
echo "   2: CreateGermlines        $(date +'%H:%M %D')"
#$RUN CreateGermlines.py -d "${BASENAME}_db-pass.tab" -r $GERMLINE_FOLDER \
#    --log GermlineLog.log >> $RUNLOG
$RUN CreateGermlines.py -d *_db-pass.tab -r $GERMLINE_FOLDER \
    --log GermlineLog.log >> $RUNLOG


# Assign clones
echo "   3: DefineClones bygroup   $(date +'%H:%M %D')"
#$RUN DefineClones.py bygroup -d "${BASENAME}_db-pass_germ-pass.tab" --model $DC_MODEL --dist $DC_DIST \
#    --log CloneLog.log >> $RUNLOG 
$RUN DefineClones.py bygroup -d *_db-pass_germ-pass.tab --model $DC_MODEL --dist $DC_DIST \
    --log CloneLog.log >> $RUNLOG 

if $ZIP_FILES; then
    tar -cf LogFiles.tar *Log.log
    gzip LogFiles.tar
    rm *Log.log
    
    tar -cf TempFiles.tar *fail.tab *db-pass.tab *germ-pass.tab
    gzip TempFiles.tar
    rm *fail.tab *db-pass.tab *germ-pass.tab
fi

# End
echo -e "DONE\n" 
cd ..