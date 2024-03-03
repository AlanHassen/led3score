#!/bin/bash
#SBATCH --job-name=led3_caspyrus20k_hist_xgboost
#SBATCH --partition="cpu-short"
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40GB
#SBATCH --mail-user="<mail>@liacs.leidenuniv.nl"
#SBATCH --mail-type="ALL"
#SBATCH --array=0-999

CURRENT_ARRAY=$((SLURM_ARRAY_TASK_ID))

ENTITY=<user>
PROJECT=led3_caspyrus20k_hist_xgboost
WANDB_ID=<wand-db-id>

echo "Current array index: $CURRENT_ARRAY"

wandb agent --count 1 $ENTITY/$PROJECT/$WANDB_ID