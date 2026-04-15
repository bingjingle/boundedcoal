#!/bin/bash
#SBATCH --job-name=DIS2_100
#SBATCH --array=1-30 
#SBATCH --time=2-00:00:00 
#SBATCH --mem=20G 
#SBATCH -p hns

ml R/4.2.0
Rscript dis_syn2_ntip100.R $SLURM_ARRAY_TASK_ID 
