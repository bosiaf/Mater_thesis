#!/bin/bash

INF_RATE_CONST=(1e-5 3e-5 5e-5 1e-4 4e-4 6.0e-4 9.0e-4 1.5e-3 4.0e-3 6.0e-3 9.0e-3 1.5e-2 2.5e-2)
FITNESS=(0.05 0.1 0.15 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1)
FIT_NON_SNP=(-0.5 -0.55 -0.7 -0.8 -0.9 -1 -1.2 -1.3 -1.5 -1.7 -1.9 -2.1 -2.3)
case=1
SEED=93110
PTO=$SCRATCH/Set2/

cp parameter_search.sh $PTO/

for i in ${INF_RATE_CONST[*]}
do
	for j in {0..12}
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
			awk 'NR==31 {$0='$i'}1' $PTO/Epidemics_$case/parameters_cluster.dat > $PTO/Epidemics_$case/par.dat
			rm $PTO/Epidemics_$case/parameters_cluster.dat
			mv $PTO/Epidemics_$case/par.dat $PTO/Epidemics_$case/parameters_cluster.dat
			awk 'NR==41 {$0='${FITNESS[$j]}'}1' $PTO/Epidemics_$case/parameters_cluster.dat > $PTO/Epidemics_$case/par.dat
			rm $PTO/Epidemics_$case/parameters_cluster.dat
			mv $PTO/Epidemics_$case/par.dat $PTO/Epidemics_$case/parameters_cluster.dat
			awk 'NR==43 {$2='${FIT_NON_SNP[$j]}'}1' $PTO/Epidemics_$case/parameters_cluster.dat > $PTO/Epidemics_$case/par.dat
			rm $PTO/Epidemics_$case/parameters_cluster.dat
			mv $PTO/Epidemics_$case/par.dat $PTO/Epidemics_$case/parameters_cluster.dat
			awk 'NR==55 {$0='$k'}1' $PTO/Epidemics_$case/parameters_cluster.dat > $PTO/Epidemics_$case/par.dat
	        rm $PTO/Epidemics_$case/parameters_cluster.dat
	        mv $PTO/Epidemics_$case/par.dat $PTO/Epidemics_$case/parameters_cluster.dat
			
			awk 'NR==61 {$0='$SEED'}1' $PTO/Epidemics_$case/parameters_cluster.dat > $PTO/Epidemics_$case/par.dat
            rm $PTO/Epidemics_$case/parameters_cluster.dat
            mv $PTO/Epidemics_$case/par.dat $PTO/Epidemics_$case/parameters_cluster.dat


			echo ${FITNESS[$j]}
			echo ${FITNESS_NON_SNP[$j]}
			((case = case + 1))
		done
	done
done
