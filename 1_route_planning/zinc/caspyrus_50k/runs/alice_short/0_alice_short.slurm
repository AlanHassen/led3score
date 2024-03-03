#!/bin/bash
#SBATCH --job-name=zinc_caspyrus_50k
#SBATCH --partition="cpu-short"
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10GB
#SBATCH --mail-user=mail@liacs.leidenuniv.nl"
#SBATCH --mail-type="ALL"
#SBATCH --array=0-4999

ROW_INDEX=$((SLURM_ARRAY_TASK_ID))

smiles_csv="<PATH>/led3_score/1_route_planning/zinc/caspyrus_50k/splitted/CASPyrus50k_part_${ROW_INDEX}.csv"
config="<PATH>/led3_score/1_route_planning/zinc/zinc_azf_config.yml"
output="<PATH>/led3_score/1_route_planning/zinc/caspyrus_50k/results/CASPyrus50k_part_${ROW_INDEX}_result.hdf"
stocks="zinc"
policy="uspto"

aizynthcli --smiles ${smiles_csv} --config ${config} --output ${output} --stocks ${stocks} --policy ${policy}