PROJECT_NAME = zinc_chembl200k_random_hist_xgboost
DATASET_PATH = <PATH>/led3_score/data/results/2_synthesis_score/train_data/chembl200k_random/chembl_200k_random_train.csv

create_folder:
	# copy template folder to new folder with the name PROJECT_NAME
	cp -r template $(PROJECT_NAME)

	# replace <PROJECT_NAME> in all files of the new folder with the name PROJECT_NAME
	find $(PROJECT_NAME) -type f -exec sed -i 's%<PROJECT_NAME>%$(PROJECT_NAME)%g' {} +
	
	# replace <DATASET_PATH> in all files of the new folder with the name DATASET_PATH
	find $(PROJECT_NAME) -type f -exec sed -i 's%<DATASET_PATH>%$(DATASET_PATH)%g' {} +

	echo "Created folder $(PROJECT_NAME) with the dataset $(DATASET_PATH)"