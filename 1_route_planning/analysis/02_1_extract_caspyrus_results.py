import pandas

# from csv_save_path, read in the data but only keep the target and number of solved routes columns

def extract_results(azf_result_csv_path):
    print("Extracting results from {}".format(azf_result_csv_path))
    results_df = pandas.read_csv(azf_result_csv_path, usecols=["target", "number_of_solved_routes"])
    return results_df

led3_azf_result_csv_path = "./led3_score/1_route_planning/led3/caspyrus_50k/led3_caspyrus50k_results.csv"
led3_azf_results_df = extract_results(led3_azf_result_csv_path)

zinc_azf_result_csv_path = "./led3_score/1_route_planning/zinc/caspyrus_50k/zinc_caspyrus50k_results.csv"
zinc_azf_results_df = extract_results(zinc_azf_result_csv_path)

original_caspyrus_path = "./led3_score/data/chemical_datasets/casp/CASPyrus50k.csv"
original_caspyrus_df = pandas.read_csv(original_caspyrus_path)

# create a column in original_caspyrus_df that is the number of solved routes for led3
original_caspyrus_df["led3_result"] = led3_azf_results_df["number_of_solved_routes"]
original_caspyrus_df["zinc_result"] = zinc_azf_results_df["number_of_solved_routes"]

result_path = "./led3_score/data/results/1_route_planning/caspyrus_50k_results.csv"
original_caspyrus_df.to_csv(result_path, index=False)

led3_result_path = "./led3_score/data/results/1_route_planning/single_results/led3_caspyrus50k_results.csv"
led3_azf_results_df.to_csv(led3_result_path, index=False)

zinc_result_path = "./led3_score/data/results/1_route_planning/single_results/zinc_caspyrus50k_results.csv"
zinc_azf_results_df.to_csv(zinc_result_path, index=False)

# assert that the led3_result and zinc_result columns are the same as the number_of_solved_routes columns in the led3 and zinc results dataframes
assert (original_caspyrus_df["led3_result"] == led3_azf_results_df["number_of_solved_routes"]).all()
assert (original_caspyrus_df["zinc_result"] == zinc_azf_results_df["number_of_solved_routes"]).all()

print("Careful, AZF doesnt produce the correct stereochemistry as output, so the smiles are not the same")

# find the differences
led3_diff = original_caspyrus_df["clean_smiles"] != led3_azf_results_df["target"]
zinc_diff = original_caspyrus_df["clean_smiles"] != zinc_azf_results_df["target"]

# print index of the differences
print("led3_diff: {}".format(led3_diff[led3_diff == True].index))
print("zinc_diff: {}".format(zinc_diff[zinc_diff == True].index))

# subset the original_caspyrus_df to only the rows where the smiles are different
orginal_diff_to_led3 = original_caspyrus_df[led3_diff]
led3_diff_to_original = led3_azf_results_df[led3_diff]

# calculate the inchis for the original and led3
from rdkit import Chem
orginal_diff_to_led3["original_inchi"] = orginal_diff_to_led3["clean_smiles"].apply(lambda x: Chem.MolToInchi(Chem.MolFromSmiles(x)))
led3_diff_to_original["led3_inchi"] = led3_diff_to_original["target"].apply(lambda x: Chem.MolToInchi(Chem.MolFromSmiles(x)))

assert (orginal_diff_to_led3["original_inchi"] == led3_diff_to_original["led3_inchi"]).all()
print("All inchis are the same for led3")

orginal_diff_to_zinc = original_caspyrus_df[zinc_diff]
zinc_diff_to_original = zinc_azf_results_df[zinc_diff]

# calculate the inchis for the original and zinc
orginal_diff_to_zinc["original_inchi"] = orginal_diff_to_zinc["clean_smiles"].apply(lambda x: Chem.MolToInchi(Chem.MolFromSmiles(x)))
zinc_diff_to_original["zinc_inchi"] = zinc_diff_to_original["target"].apply(lambda x: Chem.MolToInchi(Chem.MolFromSmiles(x)))

assert (orginal_diff_to_zinc["original_inchi"] == zinc_diff_to_original["zinc_inchi"]).all()
print("All inchis are the same for zinc")