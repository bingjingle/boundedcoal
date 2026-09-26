#!/bin/bash
#SBATCH --job-name=DIS1_100
#SBATCH --array=1-30 
#SBATCH --time=2-00:00:00 
#SBATCH --mem=20G 
#SBATCH -p hns

# Site-specific template: adjust partition, resources, and module names.
ml R/4.2.0
Rscript dis_syn1_ntip100.R "$SLURM_ARRAY_TASK_ID"
