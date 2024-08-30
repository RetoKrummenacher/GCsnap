#!/bin/bash

## script to run GCsnap on sciCORE with mpi4py

path=/scicore/home/schwede/kruret00/MT/experiments/
gcsnap_path=/scicore/home/schwede/kruret00/MT/GCsnap2/
exp_path=${path}/gcsnap2_cluster/
target_files=${path}targets/

## load Python module
ml Python

# how many repetitions
repetitions=1
cpus_per_task=2

## targets loop
#for n_targets in 10000 50000 100000 250000 500000 
for n_targets in 10000
do		
	## nodes
	for nodes in 8
	do
		## ranks on each node
		for ranks_per_node in 8 
		do		

			## repetition loop
			for (( rep=1; rep<=${repetitions}; rep++ ))
			do		
				ident=${n_targets}_${nodes}_${ranks_per_node}_${cpus_per_task}_${rep}
		
				## exported are all the arguments needed for the Python script
				## slurm configuration (e.g., --nodes) are set here when changing during iterations
				sbatch 	--export=ALL,exp_path=${exp_path},gcsnap_path=${gcsnap_path},target_files=${target_files},n_targets=${n_targets},nodes=${nodes},ranks_per_node=${ranks_per_node},cpus_per_task=${cpus_per_task},rep=${rep},ident=${ident} \
						--job-name=${ident} \
						--nodes=${nodes}  \
						--ntasks-per-node=${ranks_per_node} \
						--cpus-per-task=${cpus_per_task} \
						--output=${exp_path}run_${ident}.out \
						${exp_path}gcsnap.job
									
				sleep 0.2	
			done
		done 
	done
done 