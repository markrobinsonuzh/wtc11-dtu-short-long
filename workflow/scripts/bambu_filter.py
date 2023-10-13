import numpy as np
import pandas as pd
from xgboost import XGBClassifier
import os 

from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score
from sklearn.model_selection import GridSearchCV

os.chdir(r"/home/asamant/pacbio_benchmark/isoquant/feature_tables/all_samples_50%_anno")

data_file =  "feature_table.tsv" # snakemake.input[["feature_matrix"]]
data = pd.read_csv(data_file, delimiter='\t')
columns = data.columns

X = data[list(set(columns) - set(["label","transcript_id"]))]
y = data["label"]

scaler = StandardScaler()
scaler.fit(X)

X = scaler.transform(X)

parameters_xgboost = {"n_estimators":[50,100,250,500],"max_depth":[5,10,15],"learning_rate":[0.001,0.01,0.1]}

# model = XGBClassifier(n_estimators = 50)
# model.fit(X,y)

clf = GridSearchCV(XGBClassifier(), parameters_xgboost, scoring="f1_macro", n_jobs = 16)
clf.fit(X,y)

model = clf.best_estimator_

y_pred_probs = model.predict_proba(X)
y_pred_probs

y_pred = model.predict(X)
acc = accuracy_score(y,y_pred)

