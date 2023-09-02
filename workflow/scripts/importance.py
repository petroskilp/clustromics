import sys

# logging
sys.stderr = open(snakemake.log[0], "w")
sys.stdout = open(snakemake.log[0], "w")

import time

import pandas as pd
import numpy as np

from sklearn import preprocessing
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import confusion_matrix
from sklearn.metrics import RocCurveDisplay
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
import joblib
import matplotlib.pyplot as plt
from xgboost import XGBClassifier
from sklearn.neighbors import KNeighborsClassifier

###############################################################################
# Parse dataset
###############################################################################
df = pd.read_csv(snakemake.input.parsedcounts, sep="\t", index_col=0)

conditions = list(snakemake.config["diffexp"]["variables_of_interest"].keys())
X = df.drop(columns=conditions)

cond = snakemake.wildcards["variable"]
le = preprocessing.LabelEncoder()
le.fit(df[cond])
np.save(snakemake.output.classes, le.classes_)
y = le.transform(df[cond])


###############################################################################
# Feature selection
###############################################################################

# import XGBClassifier
from xgboost import XGBClassifier

# declare parameters
params = {
            'objective':'binary:logistic',
            'booster':'gbtree',
            'max_depth': 4,
            'alpha': 10,
            'learning_rate': 0.1,
            'n_estimators':100
        }         
               
# instantiate the classifier 
xgb_clf = XGBClassifier(**params)


print("Start feature importance for {}".format(cond))
X_train, X_test, y_train, y_test = train_test_split(X,y ,
                                random_state=104, 
                                test_size=0.25, 
                                shuffle=True)

xgb_clf.fit(X=X_train, y=y_train,
                    eval_set=[(X_train, y_train), (X_test, y_test)],
                    verbose=200)

fold_importance_df = pd.DataFrame()
fold_importance_df["feature"] = X.columns
fold_importance_df["importance"] = xgb_clf.feature_importances_

fold_importance_df.sort_values(by=["importance"], inplace=True)
fold_importance_df.to_csv(snakemake.output.importance, sep="\t", index=False)