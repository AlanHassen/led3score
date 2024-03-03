PROJECT_NAME = led3_score
AZF_FOLDER = data/aizynthfinder/

create_conda_env:
	conda env create -f condaenvironment.yaml
	###############################################################################
	# Please activate the conda environment first with: conda activate $(PROJECT_NAME) #
	###############################################################################

update_poetry:
	# update poetry lock file
	poetry config virtualenvs.create false
	poetry lock -n && poetry export --without-hashes

install_packages:
	# Please activate the conda environment first with: conda activate $(PROJECT_NAME)
	poetry config virtualenvs.create false
	poetry install -n

download_data: download_aizynthfinder_data

download_aizynthfinder_data:
	download_public_data $(AZF_FOLDER)