#!/bin/bash
#$ -l h_rt=12:00:00
#$ -P cui-buchem
#$ -N SEDJOBNAME
#$ -j y
#$ -V
#$ -o /projectnb/cui-buchem/rosadche/dnap/qm_mm/logfiles/
#$ -pe omp 1 #single core multithread line

#---------------------------------------------------------------
#                     Required Modules
#---------------------------------------------------------------
module load intel/2021.1
module load openmpi/3.1.4_intel-2021
module use /projectnb/cui-buchem/rosadche/development/plumed_make/lib/plumed/
module load plumed_2.7.4_reilly

export OMP_NUM_THREADS=$NSLOTS

CHARMMEXE=/projectnb/cui-buchem/rosadche/development/charmm_make_rmax/bin/charmm

#---------------------------------------------------------------
#                    Other Variables
#---------------------------------------------------------------
INDI=SEDRUN
JOBBASE=SEDJOBNAME
WRKDIR=SEDWRKDIR
TEMPLATE=SEDTEMPLATE

start=$SECONDS

#Check current round
last=${WRKDIR}/last.txt
if [ -f "${last}" ]; then
    cnt=$(cat ${last})
    echo "round file exists, reading in current round to do: ${cnt}"
else
    echo "round file does NOT exist, round set to 0"
    cnt=0
fi


#Check total rounds
totalrunsfile=${WRKDIR}/totalruns.txt
if [ -f "${totalrunsfile}" ]; then
    TOTALRUNS=$(cat ${totalrunsfile})
    echo "total rounds file exists, reading in current value: ${TOTALRUNS}"
else
    echo "total rounds file does NOT exist, value set to SEDTOTALROUNDS"
    TOTALRUNS=SEDTOTALROUNDS
    touch ${totalrunsfile}
    echo ${TOTALRUNS} > ${totalrunsfile}
fi

JOB=${JOBBASE}_indirun_${INDI}_round_${cnt}

#---------------------------------------------------------------
#                    Set up this round
#---------------------------------------------------------------
cd ${WRKDIR}

rm ${JOB}.inp
cp ${TEMPLATE} ${JOB}.inp

sed -i "s/SEDROUND/${cnt}/g"      ${JOB}.inp

#---------------------------------------------------------------
#                    Run this round
#---------------------------------------------------------------
# To prevent CHARMM-DFTB segmentation fault
ulimit -s unlimited

${CHARMMEXE} -i ${JOB}.inp > ${JOB}.out

#---------------------------------------------------------------
#           Check completion status
#---------------------------------------------------------------
if [ $? -eq 0 ]; then
  duration=$(( SECONDS - start ))
  if (( $duration > 3600 )); then
    let "hours=duration/3600"
    let "minutes=(duration%3600)/60"
    let "seconds=(duration%3600)%60"
  elif (( $duration > 60 )); then
    let "hours=0"
    let "minutes=(duration%3600)/60"
    let "seconds=(duration%3600)%60"
  else
    let "hours=0"
    let "minutes=0"
    let "seconds=duration"
  fi
  echo "=========================================================="
  echo "CHARMM $JOBBASE $TEMPLATE independent run ${INDI} round ${cnt} completed successfully."
  echo "End date: $(date)"
  echo "Run time: $hours hour(s), $minutes minute(s) and $seconds second(s)"
  echo "=========================================================="
else
  echo "=========================================================="
  echo "CHARMM $JOBBASE $TEMPLATE independent run ${INDI} round ${cnt} did not complete..."
  echo "End date: $(date)"
  echo "Failed on round: ${cnt}"
  echo "=========================================================="
  exit 1
fi
wait

#---------------------------------------------------------------
#       Set up next in case of desired restart
#---------------------------------------------------------------
#Increase the round
cnt=$((cnt+1))
echo $cnt > $last

#---------------------------------------------------------------
#       Resubmit
#---------------------------------------------------------------
RUNMINONE=$(( TOTALRUNS-1 ))
if [ "$cnt" -le "${RUNMINONE}" ]; then
    qsub submit_production.sh
    echo "New CHARMM $JOBBASE $TEMPLATE independent run ${INDI} round ${cnt} was submitted"
fi
