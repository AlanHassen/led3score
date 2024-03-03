#!/bin/bash
#SBATCH --job-name=LED_BB_LED3SCORE_CHEMBL_200K
#SBATCH --partition="amd-short"
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=20GB
#SBATCH --array=0-9999

ROW_INDEX=$((SLURM_ARRAY_TASK_ID))

echo "Running example $ROW_INDEX"

smiles_csv="<PATH>/led3_score/3_molecular_generation/2_evaluate_molecules/full_model/led3_building_blocks/Q99685_Led3score_chembl_200k/splitted/Q99685_chembl200k_led3score_400ep_final_cleaned_part_${ROW_INDEX}.csv"
config="<PATH>/led3_score/3_molecular_generation/2_evaluate_molecules/full_model/led3_building_blocks/led3_azf_config.yml"
output="<PATH>/led3_score/3_molecular_generation/2_evaluate_molecules/full_model/led3_building_blocks/Q99685_Led3score_chembl_200k/results/batches/Q99685_chembl200k_led3score_400ep_final_cleaned_part_${ROW_INDEX}_result.hdf"
stocks="led3"
policy="uspto"

aizynthcli --smiles ${smiles_csv} --config ${config} --output ${output} --stocks ${stocks} --policy ${policy}