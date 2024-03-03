import pandas

# from csv_save_path, read in the data but only keep the target and number of solved routes columns

def extract_results(azf_result_csv_path):
    print("Extracting results from {}".format(azf_result_csv_path))
    results_df = pandas.read_csv(azf_result_csv_path, usecols=["target", "number_of_solved_routes"])
    return results_df

led3_azf_result_csv_path = "<PATH>/led3_score/2_synthesis_score/0_test_set_generation/led3/chembl_200k_test_set/led_chembl200k_test_set_results.csv"
led3_azf_results_df = extract_results(led3_azf_result_csv_path)

zinc_azf_result_csv_path = "<PATH>/led3_score/2_synthesis_score/0_test_set_generation/zinc/chembl_200k_test_set/zinc_chembl200k_test_set_results.csv"
zinc_azf_result_df = extract_results(zinc_azf_result_csv_path)

original_chembl200k_path = "<PATH>/led3_score/data/chemical_datasets/test_sets/chembl29_random_200k_test_set.csv"
original_chembl200k_df = pandas.read_csv(original_chembl200k_path)

# create a column in original_chembl200k_df that is the number of solved routes for led3
original_chembl200k_df["led3_result"] = led3_azf_results_df["number_of_solved_routes"]
original_chembl200k_df["zinc_result"] = zinc_azf_result_df["number_of_solved_routes"]

# assert that the led3_result and zinc_result columns are the same as the number_of_solved_routes columns in the led3 and zinc results dataframes
assert (original_chembl200k_df["led3_result"] == led3_azf_results_df["number_of_solved_routes"]).all()
assert (original_chembl200k_df["zinc_result"] == zinc_azf_result_df["number_of_solved_routes"]).all()

print("Careful, AZF doesnt produce the correct stereochemistry as output, so the smiles are not the same")

# find the differences
led3_diff = original_chembl200k_df["clean_smiles"] != led3_azf_results_df["target"]
zinc_diff = original_chembl200k_df["clean_smiles"] != zinc_azf_result_df["target"]

# print index of the differences
print("led3_diff: {}".format(led3_diff[led3_diff == True].index))
print("zinc_diff: {}".format(zinc_diff[zinc_diff == True].index))

print("led3 differences")
# subset the original_chembl200k_df to only the rows where the smiles are different
orginal_diff_to_led3 = original_chembl200k_df[led3_diff]
led3_diff_to_original = led3_azf_results_df[led3_diff]

# calculate the inchis for the original and led3
from rdkit import Chem
orginal_diff_to_led3["original_inchi"] = orginal_diff_to_led3["clean_smiles"].apply(lambda x: Chem.MolToInchi(Chem.MolFromSmiles(x)))
led3_diff_to_original["led3_inchi"] = led3_diff_to_original["target"].apply(lambda x: Chem.MolToInchi(Chem.MolFromSmiles(x)))

#find the differences in the inchis
inchi_diff = orginal_diff_to_led3["original_inchi"] != led3_diff_to_original["led3_inchi"]
print("inchi_diff: {}".format(inchi_diff[inchi_diff == True].index))

print("zinc differences")
# subset the original_chembl200k_df to only the rows where the smiles are different
orginal_diff_to_zinc = original_chembl200k_df[zinc_diff]
zinc_diff_to_original = zinc_azf_result_df[zinc_diff]

# calculate the inchis for the original and zinc
orginal_diff_to_zinc["original_inchi"] = orginal_diff_to_zinc["clean_smiles"].apply(lambda x: Chem.MolToInchi(Chem.MolFromSmiles(x)))
zinc_diff_to_original["zinc_inchi"] = zinc_diff_to_original["target"].apply(lambda x: Chem.MolToInchi(Chem.MolFromSmiles(x)))

#find the differences in the inchis
inchi_diff = orginal_diff_to_zinc["original_inchi"] != zinc_diff_to_original["zinc_inchi"]
print("inchi_diff: {}".format(inchi_diff[inchi_diff == True].index))

# assert that the inchis are the same
#assert (orginal_diff_to_led3["original_inchi"] == led3_diff_to_original["led3_inchi"]).all()
#assert (orginal_diff_to_zinc["original_inchi"] == zinc_diff_to_original["zinc_inchi"]).all()
print("Note: the differences are checked with 03_check_differences.py")

print("Create Labels")

from led3_score.fingerprints.fingerprints import DrugExFingerprint
from led3_score.fingerprints.fingerprint_generator import FingerprintGenerator

drugex_fingerprint = DrugExFingerprint()
fingerprint_generator = FingerprintGenerator(fingerprint = drugex_fingerprint)

# at least one route found
original_chembl200k_df["Y_led3"] = original_chembl200k_df["led3_result"] > 0
original_chembl200k_df["Y_zinc"] = original_chembl200k_df["zinc_result"] > 0

original_chembl200k_df = fingerprint_generator.create_fingerprints_pandas_dataframe(data_frame = original_chembl200k_df, smiles_column = "clean_smiles", fingerprint_column = "fingerprint")
original_chembl200k_df = fingerprint_generator.flatten_fingerprints_in_dataframe(dataframe = original_chembl200k_df, fingerprint_column = "fingerprint")

result_path = "<PATH>/led3_score/data/results/2_synthesis_score/test_data/test_chembl200k_random/test_chembl200k_results.csv"
print("Saving results to {}".format(result_path))
#original_chembl200k_df.to_csv(result_path, index=False)