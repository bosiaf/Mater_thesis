#!/bin/bash

cores=2
case=1
HowMany=544

export OMP_NUM_THREADS=$cores
export OMP_PROC_BIND=true

while [ $case -le $HowMany ]
do
#	export GMON_OUT_PREFIX=$case
#	echo $GMON_OUT_PREFIX
	bsub -n $cores -oo Output/Epidemics_$case/output ./CEvo.o Output/Epidemics_$case/parameters_cluster.dat
	((case = case + 1))
done
