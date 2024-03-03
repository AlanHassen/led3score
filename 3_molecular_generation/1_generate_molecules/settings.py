import os
import sys
import numpy as np
import pandas as pd
sys.path.append("../../src/") # load led3score module directory

TARGET_ID = 'Q99685'

ROOT_DIR = os.path.dirname(__file__)
DATA_FOLDER = f'{ROOT_DIR}/data'
MODELS_FOLDER = f'{ROOT_DIR}/models'
DATA_GENERATED_FOLDER = f'{DATA_FOLDER}/generated'

TRAIN_DF = pd.read_csv(os.path.join(DATA_FOLDER, f'{TARGET_ID}_papyrus_low.tsv'), sep='\t', header=0)

N_PROC = 64

# data preprocessing
DATASETS_ENCODED_PATH = f"{DATA_FOLDER}/encoded"
os.makedirs(DATASETS_ENCODED_PATH, exist_ok=True)

# pretrained models
MODELS_PR_PATH = f'{MODELS_FOLDER}/pretrained/'
MODEL_FILE_PR = f'{MODELS_PR_PATH}/Papyrus05.5_graph_trans_PT/Papyrus05.5_graph_trans_PT.pkg'
VOCAB_FILE_PR = f'{MODELS_PR_PATH}/Papyrus05.5_graph_trans_PT/Papyrus05.5_graph_trans_PT.vocab'

# finetuning
N_EPOCHS_FT = 500
GPUS_FT = [0,1,2]
MODELS_FT_PATH = f'models/finetuned/{TARGET_ID}'
## path to hold all compounds up to fine-tuning:
GENERATED_FT_PATH_ALL = f'{DATA_GENERATED_FOLDER}/finetuning/{TARGET_ID}_FT_all.tsv'
os.makedirs(os.path.dirname(GENERATED_FT_PATH_ALL), exist_ok=True)

# reinforcement learning
GPUS_RL = [0,1,2]
N_EPOCHS_RL = 500
MODEL_FILE_FT = f'{MODELS_FT_PATH}/{TARGET_ID}_FT.pkg'
MODELS_RL_PATH = f"models/RL/{os.environ['DRUGEX_JOB']}" if "DRUGEX_JOB" in os.environ else f"models/RL/{TARGET_ID}"
RL_SCORERS = os.environ['DRUGEX_SCORERS'].strip('|').split('|') if 'DRUGEX_SCORERS' in os.environ else ['led3', 'qsar_cls']
print("RL_SCORERS =", RL_SCORERS)
RL_OUT_PREFIX = "-".join(RL_SCORERS)
LED3_MODEL="zinc_chembl200k_random_hist_xgboost.model" if not ('LED3_MODEL' in os.environ) else os.environ['LED3_MODEL']
DRUGEX_EPSILON=0.2 if not ('DRUGEX_EPSILON' in os.environ) else float(os.environ['DRUGEX_EPSILON'])
