from settings import *
from predictors import SCORERS, THRESHOLDS
from utils import fit_rl

scorers = [SCORERS[scorer] for scorer in RL_SCORERS]
thresholds = [THRESHOLDS[scorer] for scorer in RL_SCORERS]

environment = DrugExEnvironment(scorers, thresholds, reward_scheme=ParetoCrowdingDistance())

models = []
for env in [environment]:
    # add more instances of environment for more replicates
    models.append(
        fit_rl(
            env, 
            f'{MODELS_FT_PATH}/2023-02-03_12-26-31_027350/Q99685_FT.pkg', # exploration network
            # f'{MODELS_RL_PATH}/2023-02-03_14-58-57_756346/Q99685_RL.pkg', # agent (exploitation network), comment out to use MODEL_FILE_PR (the default pretrained network)
            epsilon=DRUGEX_EPSILON # exploration rate from settings
        )
    )
