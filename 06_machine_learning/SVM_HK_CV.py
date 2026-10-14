#!/usr/bin/env python3

"""
Python script to train an SVM model applying 5-fold stratified cross validation
and get the feature importance values
"""

import pandas as pd
import numpy as np
from sklearn.metrics import confusion_matrix
from sklearn.svm import LinearSVC
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

###Load sensor domain composition dataframe
df = pd.read_csv("/home/roman.sanz/data/diamond_pblast/australia/HK_abundance_matrix_prot_50cov_50score.csv", index_col=0,header=0)

#Filter prevalent sensor domains
with open('/home/roman.sanz/scripts/data/hk_aus_prev_filtered2.txt') as hk_prev:
    hk_prev_list=[]
    for line in hk_prev:
        hk_prev_list.append(line[:-1])

df = df[df.columns.intersection(hk_prev_list)]
drop_cols = df.columns[(df != 0).sum() < 10] #Drop low-occurrent sensor domains
df.drop(drop_cols, axis = 1, inplace = True) 


hk_annotation=pd.read_csv('/home/roman.sanz/data/HK_related/outer_hk_fam_ligand_annot_v2.csv',header=0,index_col=0)
df = df[df.columns.intersection(hk_annotation.index.to_list())]

df = df.transform(lambda x: np.log2((x*1e6)+1))

#Load metadata
metadata=pd.read_csv("/home/roman.sanz/scripts/data/Australia_merged_metadata.csv", index_col=0,header=0)

df=df.reindex(metadata.index)

x=df
habitat=metadata['general_env_feature'].to_list()
aridity=[]
for env in habitat:
    if env=='Arid':
        aridity.append(1)
    else:
        aridity.append(0)
y=aridity

y = np.asarray(y)

skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

coefs = []

for train_idx, val_idx in skf.split(x, y):
    #print(train_idx)
    X_train, X_val = x.iloc[train_idx], x.iloc[val_idx]
    y_train, y_val = y[train_idx], y[val_idx]

    model = Pipeline([
        ("scaler", StandardScaler()),
        ("svm", LinearSVC(C=100, max_iter=5000))
    ])

    model.fit(X_train, y_train)
    print(confusion_matrix(y_val,model.predict(X_val)))
    # LinearSVC coef_: (n_classes, n_features)
    coef = model.named_steps["svm"].coef_

    # Binary classification → take the single row
    if coef.shape[0] == 1:
        coef = coef.flatten()
    else:
        # Multiclass: average absolute importance across classes
        coef = np.mean(np.abs(coef), axis=0)

    coefs.append(coef)


coefs = np.array(coefs)  # shape: (n_folds, n_features)

mean_importance = np.mean(coefs, axis=0)
std_importance = np.std(coefs, axis=0)
feature_list=df.columns.to_list()


with open('/home/roman.sanz/scripts/machine_learning/svm_lin_cv_feat_imp_outerhks.txt','w') as outfile:
    outfile.write('feature,value\n')
    for i in range(len(mean_importance)):
        outfile.write(f'{feature_list[i]},{mean_importance[i]}\n')

