#!/bin/bash

## Load modlues
#ml Python/3.10.4-GCCcore-11.3.0
export PYTHONPATH=$PYTHONPATH:/scicore/home/schwede/kruret00/MT/GCsnap2

## script to to call run.job which calls db_update_sequences.py tu add sequence to db for GCsnap2.0
## Author: Reto Krummenacher

# paths
path=/scicore/home/schwede/GROUP/gcsnap_db/
exp_path=${path}create_db_scripts/
gcsnap_path=/scicore/home/schwede/kruret00/MT/GCsnap2/gcsnap/

# arguments
nodes=1
n_processes=10

ident=update_seq_${n_processes}

sbatch --export=ALL,gcsnap_path=${gcsnap_path},n_processes=${n_processes} --job-name=${ident} --nodes=${nodes} --ntasks=${n_processes} --output=${exp_path}${ident}.out ${exp_path}sequence_update.job


