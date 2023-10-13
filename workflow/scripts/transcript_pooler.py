import numpy as np
import pandas as pd
from xgboost import XGBClassifier
import os 

from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score
from sklearn.model_selection import GridSearchCV


os.chdir(r"/home/asamant/pacbio_benchmark")

def label_assigner(x):
    novel = None
    for label in x:
        if label != None:
            novel = label
        if label.startswith("SIRV"):
            return x
    
    return novel
    


sample_paths = snakemake.input["samples"]
sample_names = [sample.split("/")[3] for sample in sample_paths]

gffcmp_file_path = "gffcmp.tracking"

all_samples_probability_table = pd.read_csv(gffcmp_file_path,delimiter='\t',header=None)
all_samples_probability_table = all_samples_probability_table.iloc[:,4:]
all_samples_probability_table.columns = sample_names

for name in sample_names:
    all_samples_probability_table[name] = all_samples_probability_table[name].apply(lambda x: None if (x == None or x == "-") else x.split("|")[1])

for sample in sample_paths:

    sample_name = sample.split("/")[3]
    data = pd.read_csv(sample, delimiter='\t')
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

    probs_table = pd.DataFrame({f"{sample_name}": data["transcript_id"], f"{sample_name}_prob": y_pred_probs})
    all_samples_probability_table = all_samples_probability_table.merge(probs_table, on=f"{sample_name}")

all_samples_probability_table["transcript_id"] = all_samples_probability_table[sample_names].apply(label_assigner) 
all_samples_probability_table["max_TPS"] = all_samples_probability_table[[sample_name + "_prob" for sample_name in sample_names]].apply(lambda x: max(x))






