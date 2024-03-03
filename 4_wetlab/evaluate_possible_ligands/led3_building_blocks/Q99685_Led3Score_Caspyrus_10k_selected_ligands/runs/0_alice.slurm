#!/bin/bash
#SBATCH --job-name=LED_BB_Led3_caspyrus_10k
#SBATCH --partition="cpu-short"
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=100GB
#SBATCH --array=0

smiles_csv="<PATH>/led3_score/4_wetlab/evaluate_possible_ligands/led3_building_blocks/Q99685_Led3Score_Caspyrus_10k_selected_ligands/splitted/new_ligand.csv"
config="<PATH>/led3_score/4_wetlab/evaluate_possible_ligands/led3_building_blocks/led3_azf_config.yml"
output="<PATH>/led3_score/4_wetlab/evaluate_possible_ligands/led3_building_blocks/Q99685_Led3Score_Caspyrus_10k_selected_ligands/results/new_ligands.hdf"
stocks="led3"
policy="uspto"

aizynthcli --smiles ${smiles_csv} --config ${config} --output ${output} --stocks ${stocks} --policy ${policy}