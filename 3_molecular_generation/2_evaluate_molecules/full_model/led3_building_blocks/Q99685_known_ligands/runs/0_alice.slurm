#!/bin/bash
#SBATCH --job-name=LED_BB_known
#SBATCH --partition="amd-short"
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=30GB
#SBATCH --array=0-99

CURRENT_ARRAY=$((SLURM_ARRAY_TASK_ID))
echo "Current array index: $CURRENT_ARRAY"

smiles_csv="<PATH>/led3_score/3_molecular_generation/2_evaluate_molecules/full_model/led3_building_blocks/Q99685_known_ligands/splitted/Q99685_papyrus_known_ligands_cleaned_part_${CURRENT_ARRAY}.csv"
config="<PATH>/led3_score/3_molecular_generation/2_evaluate_molecules/full_model/led3_building_blocks/led3_azf_config.yml"
output="<PATH>/led3_score/3_molecular_generation/2_evaluate_molecules/full_model/led3_building_blocks/Q99685_known_ligands/results/Q99685_papyrus_known_ligands_cleaned_part_${CURRENT_ARRAY}_result.hdf"
stocks="led3"
policy="uspto"

aizynthcli --smiles ${smiles_csv} --config ${config} --output ${output} --stocks ${stocks} --policy ${policy}