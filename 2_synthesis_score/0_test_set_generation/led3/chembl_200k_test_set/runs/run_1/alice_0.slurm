#!/bin/bash
#SBATCH --job-name=test_led3_chembl_200k
#SBATCH --partition="amd-short"
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10GB
#SBATCH --mail-user="mail@liacs.leidenuniv.nl"
#SBATCH --mail-type="ALL"
#SBATCH --array=0-19999

ROW_INDEX=$((SLURM_ARRAY_TASK_ID))
echo "Running example $ROW_INDEX"

smiles_csv="<PATH>/led3_score/2_synthesis_score/0_test_set_generation/led3/chembl_200k_test_set/splitted/chembl29_random_200k_test_set_part_${ROW_INDEX}.csv"
config="<PATH>/led3_score/2_synthesis_score/0_test_set_generation/led3/led3_azf_config.yml"
output="<PATH>/led3_score/2_synthesis_score/0_test_set_generation/led3/chembl_200k_test_set/results/chembl29_random_200k_test_set_part_${ROW_INDEX}_result.hdf"
stocks="led3"
policy="uspto"

aizynthcli --smiles ${smiles_csv} --config ${config} --output ${output} --stocks ${stocks} --policy ${policy}