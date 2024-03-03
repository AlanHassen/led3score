from drugex.logs import logger
from settings import *

import numpy as np
import pandas as pd
import logging
import os
import copy

from rdkit import Chem
from rdkit.Chem import RDConfig
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdDepictor, rdMolDraw2D
opts = Draw.DrawingOptions()
Draw.SetComicMode(opts)

from drugex.data.corpus.vocabulary import VocGraph
from drugex.data.datasets import GraphFragDataSet
from drugex.training.models.transform import GraphModel
from drugex.training.models.explorer import GraphExplorer
from drugex.training.monitors import FileMonitor

def initLogger(filename, dir_name='data/logs/'):
    """
    Initializes a logging directory if necessary and places all DrugEx outputs in the specified file.
    
    Args:
        filename: name of the log file for DrugEx outputs
        dir_name: directory where the log file will be placed
    """
    
    filename = os.path.join(dir_name, filename)
    
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
        
    if os.path.exists(filename):
        os.remove(filename)
    
    formatter = logging.Formatter('%(asctime)s | %(levelname)s | %(message)s', 
                              '%m-%d-%Y %H:%M:%S')
    fh = logging.FileHandler(filename)
    fh.setFormatter(formatter)
    
    logger.addHandler(fh)

# different grid visualizations
standard_grid = Chem.Draw.MolsToGridImage
def interactive_grid(mols, *args, molsPerRow=5, **kwargs):
    """
    install mols2grid with pip to use
    """
    
    import mols2grid
    
    return mols2grid.display(mols, *args, n_cols=molsPerRow, **kwargs)

# show molecules as grid
make_grid = interactive_grid # change this to 'standard_grid' if you do not have the mols2grid package
def smilesToGrid(smiles, *args, molsPerRow=5, **kwargs):
    mols = []
    for smile in smiles:
        try:
            m = Chem.MolFromSmiles(smile)
            if m:
                AllChem.Compute2DCoords(m)
                mols.append(m)
            else:
                raise Exception(f'Molecule empty for SMILES: {smile}')
        except Exception as exp:
            pass
        
    return make_grid(mols, *args, molsPerRow=molsPerRow, **kwargs)

def smiles_to_df(smiles_list, names_list):
    dc = {
        "SMILES" : [],
        "Set" : []
    }
    for idx,ls in enumerate(smiles_list):
        dc["SMILES"].extend(ls)
        dc["Set"].extend([names_list[idx]] * len(ls))
    
    return pd.DataFrame(dc)

def get_timestamp(prefix=None):
    from datetime import datetime

    now = datetime.now()
    prefix = f"{prefix}_" if prefix else ""
    return f"{prefix}{now.date()}_{str(now.time()).replace(':', '-').replace('.', '_')}"

def fit_rl(environment, ft_model_path_now, pr_model_path_now=MODEL_FILE_PR, epsilon=0.1):
    vocabulary = VocGraph.fromFile(f"{DATASETS_ENCODED_PATH}/{TARGET_ID}_train.tsv.vocab")

    train = GraphFragDataSet(f"{DATASETS_ENCODED_PATH}/{TARGET_ID}_train.tsv", rewrite=False)
    test = GraphFragDataSet(f"{DATASETS_ENCODED_PATH}/{TARGET_ID}_test.tsv", rewrite=False)
    
    # pretrained model
    pretrained = GraphModel(voc_trg=vocabulary, use_gpus=GPUS_RL)
    pretrained.loadStatesFromFile(pr_model_path_now)
    
    # finetuned model
    finetuned = GraphModel(voc_trg=vocabulary, use_gpus=GPUS_RL)
    finetuned.loadStatesFromFile(f'{ft_model_path_now}')
    
    # fit
    path = f'{MODELS_RL_PATH}/{get_timestamp(prefix=RL_OUT_PREFIX)}'
    explorer = GraphExplorer(agent=pretrained, env=environment, mutate=finetuned, epsilon=epsilon, use_gpus=GPUS_RL)
    monitor = FileMonitor(f"{path}/{TARGET_ID}_RL", verbose=True)
    explorer.fit(train.asDataLoader(batch_size=512), test.asDataLoader(batch_size=512), monitor=monitor, epochs=N_EPOCHS_RL)

    return os.path.dirname(path), f"{path}/{TARGET_ID}_RL.pkg"

def show_mols_in_epoch(df_smiles, epoch=1):
    return smilesToGrid(df_smiles[df_smiles.Epoch == epoch].Smiles)

def fetch_files_from_folder(folder, extensions):
    contents = [os.path.join(folder, x) for x in os.listdir(folder)]
    next_dirs = []
    found_files = []
    for item in contents:
        if os.path.isdir(item):
            next_dirs.append(item)
            continue
        for extension in extensions:
            if item.endswith(extension):
                found_files.append(item)
            
    if next_dirs:
        for dr in next_dirs:
            found_files.extend(fetch_files_from_folder(dr, extensions))
    return found_files
