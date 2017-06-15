#!/bin/bash

INF_RATE_CONST=(1e-9 3e-9 5e-9 1e-8 4e-8 6.0e-8 9.0e-8 1.5e-7 4.0e-7 6.0e-7 9.0e-7 1.5e-6 2.5e-6)
FITNESS=(0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0)
FIT_NON_SNP=(-0.1 -0.15 -0.2 -0.25 -0.35 -0.4 -0.45 -0.5 -0.55 -0.65 -0.7 -0.75 -0.8 -0.85 -0.95 -1 -1.05 -1.1 -1.2 -1.3)
case=1
SEED=93110
PTO=$SCRATCH

cp parameter_search.sh $PTO/

for i in ${INF_RATE_CONST[*]}
do
	for j in {0..19}
	do
		for k in 0 1
		do
			echo "Case: $case"
			cp -r ./Output/Epidemics_Tmpl $PTO/Epidemics_$case
			
			awk 'NR==6 {$0="'$PTO'/Epidemics_'$case'/dyn/"}1' $PTO/Epidemics_$case/parameters_cluster.dat > $PTO/Epidemics_$case/par.dat
			rm $PTO/Epidemics_$case/parameters_cluster.dat
			mv $PTO/Epidemics_$case/par.dat $PTO/Epidemics_$case/parameters_cluster.dat
			awk 'NR==8 {$0="'$PTO'/Epidemics_'$case'/seq/"}1' $PTO/Epidemics_$case/parameters_cluster.dat > $PTO/Epidemics_$case/par.dat
			rm $PTO/Epidemics_$case/parameters_cluster.dat
			mv $PTO/Epidemics_$case/par.dat $PTO/Epidemics_$case/parameters_cluster.dat
			awk 'NR==29 {$0='$i'}1' $PTO/Epidemics_$case/parameters_cluster.dat > $PTO/Epidemics_$case/par.dat
			rm $PTO/Epidemics_$case/parameters_cluster.dat
			mv $PTO/Epidemics_$case/par.dat $PTO/Epidemics_$case/parameters_cluster.dat
			awk 'NR==39 {$0='${FITNESS[$j]}'}1' $PTO/Epidemics_$case/parameters_cluster.dat > $PTO/Epidemics_$case/par.dat
			rm $PTO/Epidemics_$case/parameters_cluster.dat
			mv $PTO/Epidemics_$case/par.dat $PTO/Epidemics_$case/parameters_cluster.dat
			awk 'NR==41 {$2='${FIT_NON_SNP[$j]}'}1' $PTO/Epidemics_$case/parameters_cluster.dat > $PTO/Epidemics_$case/par.dat
			rm $PTO/Epidemics_$case/parameters_cluster.dat
			mv $PTO/Epidemics_$case/par.dat $PTO/Epidemics_$case/parameters_cluster.dat
			awk 'NR==53 {$0='$k'}1' $PTO/Epidemics_$case/parameters_cluster.dat > $PTO/Epidemics_$case/par.dat
	        rm $PTO/Epidemics_$case/parameters_cluster.dat
	        mv $PTO/Epidemics_$case/par.dat $PTO/Epidemics_$case/parameters_cluster.dat
			
			awk 'NR==59 {$0='$SEED'}1' $PTO/Epidemics_$case/parameters_cluster.dat > $PTO/Epidemics_$case/par.dat
            rm $PTO/Epidemics_$case/parameters_cluster.dat
            mv $PTO/Epidemics_$case/par.dat $PTO/Epidemics_$case/parameters_cluster.dat


			echo ${FITNESS[$j]}
			echo ${FITNESS_NON_SNP[$j]}
			((case = case + 1))
		done
	done
done
