#!/bin/bash
#SBATCH --job-name=<PROJECT_NAME>
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
PROJECT=<PROJECT_NAME>
WANDB_ID=<SWEEP_ID>

echo "Current array index: $CURRENT_ARRAY"

wandb agent --count 1 $ENTITY/$PROJECT/$WANDB_ID