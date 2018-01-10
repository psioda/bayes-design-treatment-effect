#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem 4000
#SBATCH --output=./../cluster-out/05-%a.out
#SBATCH --error=./../cluster-err/05-%a.err
#SBATCH --array=1-600

## add SAS
module add sas/9.4


## run SAS command
sas -work /dev/shm -noterminal ./../05_calculate_power.sas -log "./../cluster-logs/05-$SLURM_ARRAY_TASK_ID.log" -sysparm "$SLURM_ARRAY_TASK_ID"