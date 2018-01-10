#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem 4000
#SBATCH --output=./../cluster-out/03-%a.out
#SBATCH --error=./../cluster-err/03-%a.err
#SBATCH --array=1-600

## add SAS
module add sas/9.4


## run SAS command
sas -work /dev/shm -noterminal ./../03_calibrate_a0.sas -log "./../cluster-logs/03-$SLURM_ARRAY_TASK_ID.log" -sysparm "$SLURM_ARRAY_TASK_ID"