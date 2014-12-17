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
DEFINE_CLONES=true

# DefineClones parameters
DC_MODEL=hs5f
DC_DIST=0.01
DC_ACT=first

# Create germlines parameters
CG_GERM=dmask
CG_FIELD=V_CALL

# Define log files
RUNLOG="Pipeline.log"
TIMELOG="Time.log"

# Make output directory and empty log files
mkdir -p $OUTDIR; cd $OUTDIR
echo '' > $RUNLOG
if $LOG_RUNTIMES; then
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
echo -e "\nSTART"
STEP=0

# Parse IMGT output
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MakeDb imgt"
#echo $RUN MakeDb.py imgt -z $IMGT_FILE -s $SEQ_FILE --outdir . --clean ">>" $RUNLOG
$RUN MakeDb.py imgt -z $IMGT_FILE -s $SEQ_FILE --outname "${OUTNAME}" \
    --outdir . --clean >> $RUNLOG

printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitDb group"
$RUN SplitDb.py group -d "${OUTNAME}_db-pass.tab" -f FUNCTIONAL >> $RUNLOG
mv "${OUTNAME}_db-pass_F.tab" "${OUTNAME}_non-functional.tab"
mv "${OUTNAME}_db-pass_T.tab" "${OUTNAME}_functional.tab"

Assign clones
if $DEFINE_CLONES; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "DefineClones bygroup"
    $RUN DefineClones.py bygroup -d "${OUTNAME}_functional.tab" --model $DC_MODEL --dist $DC_DIST \
        --mode gene --act $DC_ACT --nproc $NPROC --outname "${OUTNAME}" --log CloneLog.log >> $RUNLOG
    CG_FILE="${OUTNAME}_clone-pass.tab"
else
    CG_FILE="${OUTNAME}_functional.tab"
fi

# Create germlines
if $DEFINE_CLONES; then
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CreateGermlines"
    $RUN CreateGermlines.py -d $CG_FILE -r $GERM_DIR -g $CG_GERM \
        --vfield $CG_FIELD --cloned --outname "${OUTNAME}" --log GermLog.log >> $RUNLOG
else
    printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CreateGermlines"
    $RUN CreateGermlines.py -d $CG_FILE -r $GERM_DIR -g $CG_GERM \
        --vfield $CG_FIELD --outname "${OUTNAME}" --log GermLog.log >> $RUNLOG
fi

# Zip intermediate and log files

if $ZIP_FILES; then
    tar -cf LogFiles.tar *Log.log
    rm *Log.log
    gzip LogFiles.tar

    if $DEFINE_CLONES; then
        tar -cf TempFiles.tar *fail.tab *db-pass.tab *clone-pass.tab
        rm *fail.tab *db-pass.tab *clone-pass.tab
    else
        tar -cf TempFiles.tar *fail.tab *db-pass.tab
        rm *fail.tab *db-pass.tab
    fi
    gzip TempFiles.tar

fi

# End
echo -e "DONE\n" 
cd ..