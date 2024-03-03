#!/bin/bash
#SBATCH --job-name=LED_BB_chembl200k
#SBATCH --partition="cpu_short"
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=20GB
#SBATCH --array=0-9999

source activate led3_score

ROW_INDEX=$((SLURM_ARRAY_TASK_ID))
echo "Current array index: $CURRENT_ARRAY"

smiles_csv="<PATH>/led3_score/3_molecular_generation/2_evaluate_molecules/full_model/led3_building_blocks/Q99685_chembl_200k_zincscore/splitted/Q99685_chembl200k_zincscore_400ep_final_cleaned_part_${ROW_INDEX}.csv"
config="<PATH>/led3_score/3_molecular_generation/2_evaluate_molecules/full_model/led3_building_blocks/led3_azf_config.yml"
output="<PATH>/led3_score/3_molecular_generation/2_evaluate_molecules/full_model/led3_building_blocks/Q99685_chembl_200k_zincscore/results/Q99685_chembl200k_zincscore_400ep_final_cleaned_part_${ROW_INDEX}_result.hdf"
stocks="led3"
policy="uspto"

aizynthcli --smiles ${smiles_csv} --config ${config} --output ${output} --stocks ${stocks} --policy ${policy}