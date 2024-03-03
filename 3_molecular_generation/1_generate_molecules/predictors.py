from settings import *
from led3_score.fingerprints.fingerprints import DrugExFingerprint
from drugex.training.interfaces import Scorer
from drugex.training.scorers.modifiers import ClippedScore, SmoothClippedScore
from rdkit import Chem
import numpy as np
from xgboost import XGBClassifier
from drugex.logs import logger
import pandas as pd
from drugex.training.scorers.properties import Property
from qsprpred.scorers.predictor import Predictor
from drugex.training.scorers.modifiers import ClippedScore, SmoothClippedScore       
from tqdm.auto import tqdm

class LED3Scorer(Scorer):
    
    def __init__(self, model_path):
        super().__init__()
        self.model = XGBClassifier()
        self.model.load_model(f'models/sa_predictors/{LED3_MODEL}') # led3_caspyrus10k_hist_xgboost.model, led3_chembl200k_random_hist_xgboost.model, zinc_caspyrus10k_hist_xgboost.model or zinc_chembl200k_random_hist_xgboost.model (RAScore equivalent)
        self.fp = DrugExFingerprint()
        self.path = model_path
    
    def getScores(self, mols, frags=None):
        """
        Processes molecules and returns a score for each (i.e. a QSAR model prediction).
        """
        
        scores = []
        for mol in mols:
            if type(mol) == str:
                mol = Chem.MolFromSmiles(mol)
                
            if not mol:
                scores.append(0.0)
                continue
            
            try:
                fp = self.fp.create_rdkit_based_fingerprint(mol)
                scores.append(self.model.predict_proba([fp])[0,1]) # class one means synthesizable
            except Exception as exp:
                logger.error(f'Unknown exception while calculating fingerprints and scores in {self.getKey()}')
                logger.exception(exp)
                
        return np.array(scores)
    
    def getKey(self):
        """
        Unique Identifier among all the scoring functions used in a single environment.
        """
        
        return "LED3Scorer"
        # return f"LED3Scorer_{os.path.basename(self.path).split('.')[0]}"

def get_qsar():
    score = Predictor.fromFile(
        'models/qsar/scorers/', 
        'XGBClassifier', 
        'pchembl_value_Median_class', 
        type='CLS', 
        name='XGBClassifier',
        scale=False
    )
    score.setModifier(ClippedScore(lower_x=0.2, upper_x=0.8))
    return score

def get_qsar_reg():
    return Predictor.fromFile(
        'models/qsar/scorers', 
        'XGBRegressor', 
        'pchembl_value_Median', 
        type='REG', 
        name='XGBRegressor',
        scale=True
    )
    
def get_led3():
    score = LED3Scorer()
    score.setModifier(ClippedScore(lower_x=0.2, upper_x=0.8))
    return score

def get_sascore():
    from drugex.training.scorers.properties import Property
    from drugex.training.scorers.modifiers import SmoothClippedScore

    sascore = Property(
        "SA"
    )
    sascore.setModifier(SmoothClippedScore(lower_x=7, upper_x=4))
    return sascore

def run_predictors(mols, predictors):
    scores = {
        x.getKey() : [] for x in predictors
    }
    for predictor in tqdm(predictors):
        scores[predictor.getKey()].extend(predictor.getScores(mols))
    ret = pd.DataFrame(scores)
    ret['SMILES'] = mols
    return ret

SCORERS = {
    "led3" : get_led3(),
    "qsar_cls" : get_qsar(),
    "sascore" : get_sascore()
}

THRESHOLDS = {
    "led3" : 0.5,
    "qsar_cls" : 0.5,
    "sascore" : 0.5
}
