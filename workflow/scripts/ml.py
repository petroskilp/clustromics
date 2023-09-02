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
le.classes_ = np.load(snakemake.input.classes, allow_pickle=True)
y = le.transform(df[cond])


###############################################################################
# Models
###############################################################################

models = {
    'knn': KNeighborsClassifier(),
    'svc': SVC(),
    'lr' : LogisticRegression(), 
    'dt' : DecisionTreeClassifier(),
    'gb' : GaussianNB(),
    'rf' : RandomForestClassifier(),
    'mlp': MLPClassifier()
    }

###############################################################################
# Training pipeline
###############################################################################
print("Start process for condition: {}".format(cond))
fig = plt.figure()
ax = plt.gca()
model_name = cond = snakemake.wildcards["model_name"]
    
print(100*'*')
print("Start process for model: {}".format(model_name))

model = models[model_name]

reports = []

    
start = time.time()
X_train, X_test, y_train, y_test = train_test_split(X,y ,
                            random_state=104, 
                            test_size=0.25, 
                            shuffle=True)

# Train classifier
train_start = time.time()
model.fit(X_train, y_train)
train_finish = time.time() - train_start

print("Trained model: {} in {}".format(model_name, train_finish))

# Predict on test set
y_pred = model.predict(X_test)

print("Predicted model: {} in {}".format(model_name, time.time() - start))
score = accuracy_score(y_test, y_pred)
print("Accuracy {}: {}".format(model_name, score))
cm=confusion_matrix(y_test, y_pred)
print(cm)
#plot_roc_curve(model, X_test, y_test, ax=ax)
RocCurveDisplay.from_estimator(model, X_test, y_test, ax=ax)
joblib.dump(model, snakemake.output.model)
plt.savefig(snakemake.output.roc)