#!/bin/bash
IN=$1

# get the filename without the extension
JOB=${IN%.*}

SUBMIT=${JOB}.sub
PWD=`pwd`

cat > $SUBMIT <<!EOF
#!/bin/bash -l
#$ -l h_rt=150:00:00
#$ -P cui-buchem
#$ -N $JOB
#$ -e ${JOB}.err
#$ -j y
#$ -V
#$ -pe omp 1 #Change cores here

export OMP_NUM_THREADS=\$NSLOTS

cd $PWD

module purge
module load intel/2021.1
module load openmpi/3.1.4_intel-2021
module use /projectnb/cui-buchem/rosadche/development/plumed_make/lib/plumed/
module load plumed_2.7.4_reilly

#This is to prevent segfault with CHARMM-DFTB
ulimit -s unlimited

CHARMMEXE=/projectnb/cui-buchem/rosadche/development/charmm_make_sugar/bin/charmm
\$CHARMMEXE -i ${JOB}.inp > ${JOB}.out

!EOF

qsub $SUBMIT
rm $SUBMIT