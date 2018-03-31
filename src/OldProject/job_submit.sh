#!/bin/bash

cores=2
case=1
HowMany=160
PTO=$SCRATCH/1000bp_10SNP_doubleMalus/
export OMP_NUM_THREADS=$cores
export OMP_PROC_BIND=true

while [ $case -le $HowMany ]
do
#	export GMON_OUT_PREFIX=$case
#	echo $GMON_OUT_PREFIX
	bsub -n $cores -W 18:00 -oo $PTO/Epidemics_$case/output.txt ./CEvo.o $PTO/Epidemics_$case/parameters_cluster.dat
	((case = case + 1))
done
