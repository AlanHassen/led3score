#!/bin/bash
#SBATCH --job-name=ZINC_BB_Led3_caspyrus_10k
#SBATCH --partition="gpu-short"
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=40GB
#SBATCH --array=0-9999

ROW_INDEX=$((SLURM_ARRAY_TASK_ID))

echo "Running example $ROW_INDEX"

smiles_csv="<PATH>/led3_score/3_molecular_generation/2_evaluate_molecules/full_model/led3_building_blocks/Q99685_Led3_caspyrus_10k/splitted/Q99685_casp10k_led3score_400ep_final_cleaned_part_${ROW_INDEX}.csv"
config="<PATH>/led3_score/3_molecular_generation/2_evaluate_molecules/full_model/zinc_building_blocks/zinc_azf_config.yml"
output="<PATH>/led3_score/3_molecular_generation/2_evaluate_molecules/full_model/zinc_building_blocks/Q99685_Led3_caspyrus_10k/results/Q99685_casp10k_led3score_400ep_final_cleaned_part_${ROW_INDEX}_result.hdf"
stocks="zinc"
policy="uspto"

aizynthcli --smiles ${smiles_csv} --config ${config} --output ${output} --stocks ${stocks} --policy ${policy}