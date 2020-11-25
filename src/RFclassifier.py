#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 10:09:03 2019

@author: wuyuanqi
"""

import numpy as np
import pandas as pd
import random as rd
import os
import copy
import re
import matplotlib.pyplot as plt
from sklearn import svm
from scipy import interp
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import SelectFromModel
from sklearn.datasets import make_classification
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedShuffleSplit, StratifiedKFold
from sklearn.svm import SVC
import cloudpickle as pickle
from scipy import stats
from statsmodels.stats.multitest import multipletests
from scipy.spatial.distance import braycurtis
import scipy.stats as st

 
#modeling
def get_kfold_auc(data, meta_group, label, max_features=0.1, n_estimators=501, k=10):#, min_samples_leaf=80
    aucs = []
    tprs = []
    mean_fpr = np.linspace(0, 1, 100)
    plot_data = []
    i = 0
    splitor = StratifiedKFold(n_splits=k, shuffle=True,random_state=SEED) 
#    sample_leaf_options = [1,5,10,50,100,200,500]
    clf = RandomForestClassifier(n_estimators = n_estimators, oob_score = True, random_state =SEED,
                                max_features = max_features)#, min_samples_leaf = min_samples_leaf
    
    for train_index, test_index in splitor.split(data, meta_group):
        y_train, y_test = label[train_index], label[test_index]
        X_train, X_test = np.array(data)[train_index], np.array(data)[test_index]
        probas = clf.fit(X_train, y_train).predict_proba(X_test)
        fpr, tpr, thresholds = roc_curve(y_test, probas[:, 1])
        roc_auc = auc(fpr, tpr)
        aucs.append(roc_auc)
        ### plot data
        tprs.append(interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0
        plot_data.append([fpr, tpr, 'ROC Fold %d(AUC = %0.2f)' %(i+1, roc_auc)])
        i += 1
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    return mean_auc, clf, (plot_data, mean_fpr, mean_tpr, tprs, np.std(aucs))


##########################
def dataset_reader(ds='an'): # ds: an, ca ,cn
    study_ids = ['6070', '290926', '389927'] if ds!='ca' else ['6070', '290926', '362366', '389927']
    rm_state = 'C' if ds=='an' else ('A' if ds=='cn' else 'N')
    file_name = 'AN' if ds == 'an' else ('CN' if ds == 'cn' else 'CA')
    control_state = 'N' if ds!='ca' else 'A'
    ### read dataset
    data = pd.read_csv('OTU_del_all_shannon2.csv', index_col=0)
    label = np.array([i[0] for i in data.index])
    data_sub = data.loc[label!=rm_state, :]
    data = data_sub
    data = data.iloc[:, 1:]
    ### read features wilcoxon differential ASVs
    diff = pd.read_csv(file_name+ '/p.val_'+ds+'_block_otus.csv', index_col=0)
    select_features = list(diff.loc[diff['all']<0.05, :].index)
    print('### Wilcoxon test :', len(select_features))    
    meta_cols = ['Age', 'Gender', 'BMI', 
                'Shannon', 'Simpson', 'InSimpson', 'spe_number']
    select_features.extend(meta_cols)
    data = data.loc[:, select_features]
    return study_ids, control_state, data


###########
###########################
## IFE
select = list(data.columns)
best_auc = 0
best_plot_data = []
best_features = []
feature_rank = []
while(len(select)>1):
    aucs = []
    for ni in select:
        temp = copy.deepcopy(select)
        temp.remove(ni)
        roc_auc, _, plot_data = get_kfold_auc(data.loc[:, temp], meta_group, label,max_features=0.1, n_estimators=501, k=10)#,min_samples_leaf=80
        aucs.append([temp, roc_auc, plot_data])
        #print(temp, roc_auc)
    select, roc_auc, plot_data = sorted(aucs, key=lambda x:x[1], reverse = True)[0]
    if roc_auc >= best_auc:
        best_auc = roc_auc
        best_features = select
        best_plot_data = plot_data
    feature_rank.append([select, roc_auc])
    print('### Best AUC :', len(select), round(best_auc, 3), round(roc_auc, 3))


#pickle.dump([SEED, best_auc, best_features, best_plot_data, feature_rank], open('/Test_outs_'+ds+'501.pkl, 'wb'))

