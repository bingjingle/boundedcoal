#!/bin/bash
#SBATCH --job-name=real100_DIS_BM
#SBATCH --array=0-15
#SBATCH --time=2-00:00:00
#SBATCH --mem=20G
#SBATCH -p stat
#SBATCH --output=slurm-%A_%a.out

ml python/3.9.0

# 4 const values x 4 common nsim1=nsim2 settings = 16 tasks
const_values=(0.1 0.5 1 5)
nsim_values=(100000 200000 300000 500000)

task_id=$SLURM_ARRAY_TASK_ID
const_idx=$(( task_id % 4 ))
nsim_idx=$(( task_id / 4 ))

const_value=${const_values[$const_idx]}
nsim_value=${nsim_values[$nsim_idx]}

echo "task=$task_id const=$const_value nsim1=$nsim_value nsim2=$nsim_value"

python3 real100_DIS_BM_const_nsim.py \
  --const "$const_value" \
  --nsim1 "$nsim_value" \
  --nsim2 "$nsim_value" \
  --output-dir "/scratch/groups/juliapr/output_Bingjing/boundedcoaldata/real/"
