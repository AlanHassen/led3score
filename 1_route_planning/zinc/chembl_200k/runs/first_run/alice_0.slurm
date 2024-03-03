#!/bin/bash
#SBATCH --job-name=zinc_chembl_200k
#SBATCH --partition="cpu-short"
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10GB
#SBATCH --mail-user="mail@liacs.leidenuniv.nl"
#SBATCH --mail-type="ALL"
#SBATCH --array=0-19999

ROW_INDEX=$((SLURM_ARRAY_TASK_ID))

echo "Running example $ROW_INDEX"

smiles_csv="<PATH>/led3_score/1_route_planning/zinc/chembl_200k/splitted/chembl29_200k_random_sampled_molecules_part_${ROW_INDEX}.csv"
config="<PATH>/led3_score/1_route_planning/zinc/zinc_azf_config.yml"
output="<PATH>/led3_score/1_route_planning/zinc/chembl_200k/results/chembl29_200k_random_sampled_molecules_part_${ROW_INDEX}_result.hdf"
stocks="zinc"
policy="uspto"

aizynthcli --smiles ${smiles_csv} --config ${config} --output ${output} --stocks ${stocks} --policy ${policy}