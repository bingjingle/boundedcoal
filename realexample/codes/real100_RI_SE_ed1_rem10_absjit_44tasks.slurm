#!/bin/bash
#SBATCH --job-name=real100_RI_SE_ed1_rem10_absjit
#SBATCH --array=0-43
#SBATCH --time=7-00:00:00
#SBATCH --mem=20G
#SBATCH -p hns
#SBATCH --output=slurm-%A_%a.out

ml python/3.9.0

# 11 proposal step sizes x 4 common nsim values = 44 tasks
proposal_values=(5e-6 5e-5 5e-4 1e-3 5e-3 5e-2 1e-2 0.1 0.15 0.2 0.3)
nsim_values=(100000 200000 300000 500000)

task_id=$SLURM_ARRAY_TASK_ID
proposal_idx=$(( task_id % 11 ))
nsim_idx=$(( task_id / 11 ))

proposal_sigma_l=${proposal_values[$proposal_idx]}
nsim=${nsim_values[$nsim_idx]}

echo "task=$task_id proposal=$proposal_sigma_l nsim1=$nsim nsim2=$nsim"

python3 real100_RI_SE_ed1_rem10_absjit.py \
  --proposal-sigma-l "$proposal_sigma_l" \
  --jitter 1e-14 \
  --nsim1 "$nsim" \
  --nsim2 "$nsim" \
  --rem 10 \
  --output-dir "/scratch/groups/juliapr/output_Bingjing/boundedcoaldata/real/"
