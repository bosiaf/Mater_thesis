#!/bin/bash

INF_RATE_CONST=(3.0e-8 6.0e-8 9.0e-8 1.5e-7 4.0e-7 6.0e-7 9.0e-7 1.5e-6 4.0e-6 6.0e-6 9.0e-6)
FITNESS=(0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.7 0.8)
FIT_NON_SNP=(-0.3 -0.35 -0.4 -0.45 -0.5 -0.55 -0.6 -0.65 -0.7 -0.8 -0.9)
case=1


for i in ${INF_RATE_CONST[*]}
do
	for j in {0..10}
	do
		echo "Case: $case"
		cp -r ./Output/Epidemics_Tmpl ./Output/Epidemics_$case
		awk 'NR==6 {$0="./Output/Epidemics_'$case'/dyn/"}1' ./Output/Epidemics_$case/parameters_cluster.dat > ./Output/Epidemics_$case/parameters_cluster.dat
		awk 'NR==8 {$0="./Output/Epidemics_'$case'/seq/"}1' ./Output/Epidemics_$case/parameters_cluster.dat > ./Output/Epidemics_$case/parameters_cluster.dat
		awk 'NR==29 {$0='$i'}1' ./Output/Epidemics_$case/parameters_cluster.dat > ./Output/Epidemics_$case/parameters_cluster.dat
		awk 'NR==39 {$0='${FITNESS[$j]}'}1' ./Output/Epidemics_$case/parameters_cluster.dat > ./Output/Epidemics_$case/parameters_cluster.dat
		awk 'NR==41 {$2='${FIT_NON_SNP[$j]}'}1' ./Output/Epidemics_$case/parameters_cluster.dat > ./Output/Epidemics_$case/parameters_cluster.dat
		echo ${FITNESS[$j]}
		echo ${FITNESS_NON_SNP[$j]}
		((case = case + 1))
	done
done
