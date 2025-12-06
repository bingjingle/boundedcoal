#!/bin/bash
#SBATCH --job-name=sny1BM50
#SBATCH --array=0-29
#SBATCH --time=0-2:00:00
#SBATCH --mem=20G
#SBATCH -p hns

ml python/3.9.0
python3 boundedcoal_BM_std_syn1_50_jitter7.py $SLURM_ARRAY_TASK_ID
