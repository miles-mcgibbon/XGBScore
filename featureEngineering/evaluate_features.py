'''
Uses sklearns cross-validated recursive feature elimination function
with a random forest classifier estimator to determine optimum
number of features to use based on the highest acheived ROC AUC
score.
'''

import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
import sys
import pandas as pd
import numpy as np

def parse_args(args): # parse CLI user inputs

    h5_path = args[args.index('-h5') + 1]

    estimator_jobs =  int(args[args.index('-est_jobs') + 1])

    estimator_trees = int(args[args.index('-est_trees') + 1])

    selector_jobs =  int(args[args.index('-sel_jobs') + 1])

    steps =  int(args[args.index('-steps') + 1])

    cvs =  int(args[args.index('-cvs') + 1])

    output_dir =  args[args.index('-out') + 1]

    return h5_path, estimator_jobs, estimator_trees, selector_jobs, steps, cvs, output_dir

def main(): # run script using CLI

    h5_path, estimator_jobs, estimator_trees, selector_jobs, steps, cvs, output_dir = parse_args(sys.argv)

    # load dataset
    df = pd.read_hdf(h5_path, key='df', mode='r')
    label_headers = ['PDBCode','BindingDataType', 'BindingValue','BindingUnits','Label', 'Database', 'Binding_Data_in_uM']
    X = df.drop(label_headers , axis=1)
    headers = X.columns
    y = df['Label'].copy()

    # define model
    model = RandomForestClassifier(n_estimators=estimator_trees, n_jobs=estimator_jobs)

    # define recursive eliminator
    selector = RFECV(model, step=steps, cv=cvs, n_jobs=selector_jobs, scoring='roc_auc', verbose=10)

    # run elimination
    selector = selector.fit(X, y)

    # fetch boolean list of features to include
    print(selector.support_)
    choices = list(selector.support_)

    # fetch feature rankings
    print(selector.ranking_)
    ranks = list(selector.ranking_)

    # return optimal feature counts
    print("Optimal number of features : %d" % selector.n_features_)

    # save importance, boolean and ranking results to dataframe
    results = pd.DataFrame({'Feature':headers, 'includeBool':choices, 'Rank':ranks})
    results.to_csv(f'{output_dir}RFECVFeatureRanks.csv', index=False)

    # save plot of aurocs to number of features
    plt.figure()
    plt.xlabel("Number of features selected")
    plt.ylabel("Cross validation score (nb of correct classifications)")
    plt.plot(np.linspace(1, len(headers), len(selector.grid_scores_)),
             selector.grid_scores_)
    plt.savefig(f'{output_dir}featureROCAUCPlot.png')

    # save aurocs of feature counts to dataframe
    feature_numbers = np.linspace(1, len(headers), len(selector.grid_scores_))
    auc_scores = selector.grid_scores_
    auc_results = pd.DataFrame({'Num Features':feature_numbers, 'ROC_AUC':auc_scores})
    auc_results.to_csv(f'{output_dir}Ideal{selector.n_features_}ROCAUCResults.csv', index=False)

if __name__ == '__main__':
    main()
