#!/bin/bash
# Super script to run the Change-O clonal assignment pipeline on IMGT data 
# 
# Author:  Jason Anthony Vander Heiden
# Date:    2014.11.22
# 
# Required Arguments:
#   $1 = IMGT zip file 
#   $2 = FASTA file that was submitted to IMGT
#   $3 = the folder containing the IMGT reference database sequences
#   $4 = output directory
#   $5 = output file prefix
#   $6 = number of subprocesses for multiprocessing tools


# Capture command line parameters
IMGT_FILE=$(readlink -f $1)
SEQ_FILE=$(readlink -f $2)
GERM_DIR=$(readlink -f $3)
OUTDIR=$(readlink -f $4)
OUTNAME=$5
NPROC=$6

# Define run parameters
LOG_RUNTIMES=false
ZIP_FILES=false
DC_MODEL=m1n
DC_DIST=5
DC_ACT=set
CG_GERM=dmask
CG_FIELD=V_CALL

# Define script execution command and log files
mkdir -p $OUTDIR; cd $OUTDIR
RUNLOG="Pipeline.log"
echo '' > $RUNLOG
if $LOG_RUNTIMES; then
    TIMELOG="Time.log"
    echo '' > $TIMELOG 
    RUN="nice -19 /usr/bin/time -o ${TIMELOG} -a -f %C\t%E\t%P\t%Mkb"
else
    RUN="nice -19"
fi
        
# Start
echo "DIRECTORY: ${OUTDIR}"
echo "VERSIONS:"
echo "  $(AnalyzeAa.py -v 2>&1)"
echo "  $(CreateGermlines.py -v 2>&1)"
echo "  $(DefineClones.py -v 2>&1)"
echo "  $(MakeDb.py -v 2>&1)"
echo "  $(ParseDb.py -v 2>&1)"
echo "  $(SplitDb.py -v 2>&1)"

# Parse IMGT
echo -e "\nSTART"
STEP=0

printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MakeDb imgt"
#echo $RUN MakeDb.py imgt -z $IMGT_FILE -s $SEQ_FILE --outdir . --clean ">>" $RUNLOG
$RUN MakeDb.py imgt -z $IMGT_FILE -s $SEQ_FILE --outname "${OUTNAME}" \
    --outdir . --clean >> $RUNLOG

printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitDb group"
$RUN SplitDb.py group -d "${OUTNAME}_db-pass.tab" -f FUNCTIONAL >> $RUNLOG
mv "${OUTNAME}_db-pass_F.tab" "${OUTNAME}_non-functional.tab"
mv "${OUTNAME}_db-pass_T.tab" "${OUTNAME}_functional.tab"

# Assign clones
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "DefineClones bygroup"
$RUN DefineClones.py bygroup -d "${OUTNAME}_functional.tab" --model $DC_MODEL --dist $DC_DIST \
    --mode gene --act $DC_ACT --nproc $NPROC --outname "${OUTNAME}" --log CloneLog.log >> $RUNLOG

# Create germlines
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CreateGermlines"
$RUN CreateGermlines.py -d "${OUTNAME}_clone-pass.tab" -r $GERM_DIR -g $CG_GERM \
    --vfield $CG_FIELD --cloned --outname "${OUTNAME}" --log GermLog.log >> $RUNLOG

if $ZIP_FILES; then
    tar -cf LogFiles.tar *Log.log
    gzip LogFiles.tar
    rm *Log.log
    
    tar -cf TempFiles.tar *fail.tab *db-pass.tab *clone-pass.tab
    gzip TempFiles.tar
    rm *fail.tab *db-pass.tab *clone-pass.tab
fi

# End
echo -e "DONE\n" 
cd ..