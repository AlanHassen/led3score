#!/bin/bash
#SBATCH --job-name=led3_chembl_200k
#SBATCH --partition="cpu-short"
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10GB
#SBATCH --mail-user="mail"
#SBATCH --mail-type="ALL"
#SBATCH --array=0-19999

ROW_INDEX=$((SLURM_ARRAY_TASK_ID))
echo "Running example $ROW_INDEX"

smiles_csv="<PATH>/paper/led3_score/1_route_planning/led3/chembl_200k/splitted/chembl29_200k_random_sampled_molecules_part_${ROW_INDEX}.csv"
config="<PATH>/paper/led3_score/1_route_planning/led3/led3_azf_config.yml"
output="<PATH>/paper/led3_score/1_route_planning/led3/chembl_200k/results/chembl29_200k_random_sampled_molecules_part_${ROW_INDEX}_result.hdf"
stocks="led3"
policy="uspto"

aizynthcli --smiles ${smiles_csv} --config ${config} --output ${output} --stocks ${stocks} --policy ${policy}