## We need a test set of chembl molecules that are not in the Chembl200k random sample nor in the caspyrus dataset

chembl29_path = "<PATH>/led3_score/data/chemical_datasets/full_datasets/cleanChembl29.csv"
chembl200k_random_train_sample_path = "<PATH>/led3_score/data/chemical_datasets/casp/chembl29_200k_random_sampled_molecules.csv"
caspyrus50k_path = "<PATH>/led3_score/data/chemical_datasets/casp/CASPyrus50k.csv"

import pandas

# Read in the datasets
chembl29 = pandas.read_csv(chembl29_path)
chembl200k_random_train_sample = pandas.read_csv(chembl200k_random_train_sample_path)
caspyrus50k = pandas.read_csv(caspyrus50k_path)

from rdkit import Chem
import swifter

chembl29["inchi"] = chembl29["clean_smiles"].swifter.apply(lambda x: Chem.MolToInchi(Chem.MolFromSmiles(x)))
chembl200k_random_train_sample["inchi"] = chembl200k_random_train_sample["clean_smiles"].swifter.apply(lambda x: Chem.MolToInchi(Chem.MolFromSmiles(x)))
caspyrus50k["inchi"] = caspyrus50k["clean_smiles"].swifter.apply(lambda x: Chem.MolToInchi(Chem.MolFromSmiles(x)))

# sample 200k molecules that are not in the chembl200k random sample nor in the caspyrus dataset
chembl200k_test_set = chembl29[~chembl29["inchi"].isin(chembl200k_random_train_sample["inchi"])]
chembl200k_test_set = chembl200k_test_set[~chembl200k_test_set["inchi"].isin(caspyrus50k["inchi"])]
chembl200k_test_set = chembl200k_test_set.sample(n=200000, random_state=42)

# assert that the test set does not contain any of the molecules in the chembl200k random sample nor in the caspyrus dataset
assert not chembl200k_test_set["clean_smiles"].isin(chembl200k_random_train_sample["clean_smiles"]).any()
assert not chembl200k_test_set["clean_smiles"].isin(caspyrus50k["clean_smiles"]).any()

# assert that the test set does not contain any of the molecules in the chembl200k random sample nor in the caspyrus dataset
assert not chembl200k_test_set["inchi"].isin(chembl200k_random_train_sample["inchi"]).any()
assert not chembl200k_test_set["inchi"].isin(caspyrus50k["inchi"]).any()

# save the test set
#chembl200k_test_set.to_csv("<PATH>/led3_score/data/chemical_datasets/test_sets/chembl29_random_200k_test_set.csv", index=False) 