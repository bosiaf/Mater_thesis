#!/bin/bash

INF_RATE_CONST=(1e-4 2e-4 3e-4 4e-4 5e-4 6e-4 8e-4 1e-3)
FITNESS=(0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5)
FIT_NON_SNP=(-0.1 -0.2 -0.3 -0.4 -0.5 -0.6 -0.7 -0.8 -0.9 -1.0)
case=1
SEED=42
SEQ_LEN=500
FRACTION_SNP=20
PTO=$SCRATCH/500bp_5SNP_doubleMalus/



rm seq_${SEQ_LEN}.dat

for ((i = 0 ; i < SEQ_LEN ; i++))
do
	printf "A" >> seq_${SEQ_LEN}.dat
done

mkdir $PTO

cp parameter_search.sh $PTO/

for i in ${INF_RATE_CONST[*]}
do
	for j in {0..9}
	do
		for h in $j
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
			
				awk 'NR==10 {$0="seq_'$SEQ_LEN'.dat"}1' $PTO/Epidemics_$case/parameters_cluster.dat > $PTO/Epidemics_$case/par.dat
				rm $PTO/Epidemics_$case/parameters_cluster.dat
				mv $PTO/Epidemics_$case/par.dat $PTO/Epidemics_$case/parameters_cluster.dat
				
				awk 'NR==12 {$0="1-'$((SEQ_LEN / FRACTION_SNP))'"}1' $PTO/Epidemics_$case/parameters_cluster.dat > $PTO/Epidemics_$case/par.dat
				rm $PTO/Epidemics_$case/parameters_cluster.dat
				mv $PTO/Epidemics_$case/par.dat $PTO/Epidemics_$case/parameters_cluster.dat
				
				
				awk 'NR==31 {$0='$i'}1' $PTO/Epidemics_$case/parameters_cluster.dat > $PTO/Epidemics_$case/par.dat
				rm $PTO/Epidemics_$case/parameters_cluster.dat
				mv $PTO/Epidemics_$case/par.dat $PTO/Epidemics_$case/parameters_cluster.dat
				awk 'NR==41 {$0='${FITNESS[$j]}'}1' $PTO/Epidemics_$case/parameters_cluster.dat > $PTO/Epidemics_$case/par.dat
				rm $PTO/Epidemics_$case/parameters_cluster.dat
				mv $PTO/Epidemics_$case/par.dat $PTO/Epidemics_$case/parameters_cluster.dat
				awk 'NR==43 {$2='${FIT_NON_SNP[$h]}'}1' $PTO/Epidemics_$case/parameters_cluster.dat > $PTO/Epidemics_$case/par.dat
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
done
