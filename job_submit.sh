#!/bin/bash

cores=4
case=1
HowMany=338
PTO=$SCRATCH/Set2

export OMP_NUM_THREADS=$cores
export OMP_PROC_BIND=true

while [ $case -le $HowMany ]
do
#	export GMON_OUT_PREFIX=$case
#	echo $GMON_OUT_PREFIX
	bsub -n $cores -oo $PTO/Epidemics_$case/output.txt ./CEvo.o $PTO/Epidemics_$case/parameters_cluster.dat
	((case = case + 1))
done
