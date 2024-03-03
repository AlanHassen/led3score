# DrugEx Training and Molecule Generation

This folder contains all code needed to train a DrugEx model, its scorers, generate new molecules and analyze data. The required environment can be obtained by installing packages from the `requirements.txt` file. In order to reproduce any experiments performed for the MAGL target, download and extract the following files with data and models:

The models and data associated with this paper will be made available upon its publication or upon request.

## Fine-Tuning

In the first step, a fine-tuned model is created to reproduce compounds from the chemical space of the target with [`FT.ipynb`](FT.ipynb). The generated molecules from this model and training can be analyzed in [`FT_analysis.ipynb`](FT_analysis.ipynb).

## Reinforcement Learning

Scorers are available in the [`./models`](./models) folder where QSAR models are also trained. The [`RL_initial.ipynb`](RL_initial.ipynb) notebook contains some code to further analyze FT data to help decide thresholds for different scores, but is not strictly necessary. [`RL.ipynb`](RL.ipynb) then implements reinforcement learning via DrugEx itself, but the [`RL_run.py`](RL_run.py) script can be used as well for distributed computing setup. [`RL_analysis_overall.ipynb`](RL_analysis_overall.ipynb) then contains code to analyze the RL training data and [`RL_generate.ipynb`](RL_generate.ipynb) generates molecules and filters them to create an initial shortlist for retrosynthetic route planning and further filtering in [`../2_evaluate_molecules`](../2_evaluate_molecules).