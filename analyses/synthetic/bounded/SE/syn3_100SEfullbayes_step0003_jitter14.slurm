#!/bin/bash
#SBATCH --job-name=sny3BM100
#SBATCH --array=0-29
#SBATCH --time=0-20:00:00
#SBATCH --mem=20G
#SBATCH -p hns

# Site-specific template: adjust partition, resources, and module names.
ml python/3.9.0
python3 boundedcoal_SEsyn3_100fullbayes_step0003_ite100000_jitter14.py "$SLURM_ARRAY_TASK_ID"
