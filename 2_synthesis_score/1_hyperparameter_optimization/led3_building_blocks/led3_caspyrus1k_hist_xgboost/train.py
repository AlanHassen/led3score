from sklearn.model_selection import train_test_split, StratifiedKFold
from led3_score.fingerprints.fingerprint_generator import FingerprintGenerator
import numpy
import xgboost as xgb
import sklearn.metrics
import pandas
import wandb
import yaml

def k_fold(dataset, learning_rate = 0.3, max_depth = 6, n_estimators = 100, min_split_loss = 0, y_column="Y_led3", number_of_folds=5):
    cv = StratifiedKFold(n_splits=number_of_folds, shuffle=True, random_state=42)
    X = FingerprintGenerator.get_fingerprint_columns(dataset, "fp_").to_numpy()
    Y = numpy.vstack(dataset[y_column].to_numpy()).ravel()
    
    X_test = []
    Y_test = []
    Y_hat = []

    for i, (train, test) in enumerate(cv.split(X, Y)):
        print(f"Fold {i}")
        model = xgb.XGBClassifier(objective="binary:logistic", tree_method="hist", random_state=42, n_jobs = 1, learning_rate = learning_rate, max_depth = max_depth, n_estimators = n_estimators, min_split_loss = min_split_loss)
        model.fit(X[train], Y[train])
        y_hat = model.predict(X[test])

        X_test.append(X[test])
        Y_test.append(Y[test])
        Y_hat.append(y_hat)

    # assert that Y_test and Y have the same entries
    assert numpy.all(numpy.isin(numpy.concatenate(X_test), X))
    assert numpy.all(numpy.isin(numpy.concatenate(Y_test), Y))

    # flatten the arrays
    Y_test = numpy.vstack(Y_test).ravel()
    Y_hat = numpy.vstack(Y_hat).ravel()

    results = {}

    results["accuracy"] = sklearn.metrics.accuracy_score(Y_test, Y_hat)
    results["precision"] = sklearn.metrics.precision_score(Y_test, Y_hat)
    results["recall"] = sklearn.metrics.recall_score(Y_test, Y_hat)
    results["f1"] = sklearn.metrics.f1_score(Y_test, Y_hat)
    results["mcc"] = sklearn.metrics.matthews_corrcoef(Y_test, Y_hat)

    return results


def main():
    
    with open('./wandb_config.yaml') as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
    
    run = wandb.init(config=config)

    # note that we define values from `wandb.config` instead of 
    # defining hard values
    learning_rate  =  wandb.config.learning_rate
    max_depth = wandb.config.max_depth
    n_estimators = wandb.config.n_estimators
    min_split_loss = wandb.config.min_split_loss

    dataset_path = wandb.config.dataset_path

    # load data
    training_csv = pandas.read_csv(dataset_path)

    results = k_fold(dataset= training_csv, learning_rate = learning_rate, max_depth = max_depth, n_estimators = n_estimators, min_split_loss = min_split_loss)
    wandb.log({
      'accuracy': results["accuracy"],
      'precision': results["precision"], 
      'recall': results["recall"], 
      'f1': results["f1"],
      'mcc': results["mcc"]
    })

main()