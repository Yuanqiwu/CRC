#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 25 15:58:15 2019

@author: wuyuanqi
"""

import numpy as np
import pandas as pd
import random as rd
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
from sklearn.model_selection import KFold


#####study-to-study
def get_kfold_auc(data, id1, id2, n_estimators=501, max_features=0.1):
    aucs = []
    tprs = []
    mean_fpr = np.linspace(0, 1, 100)
    plot_data = []
    i = 0
    clf = RandomForestClassifier(n_estimators=n_estimators, oob_score=True, max_features=max_features, random_state=SEED) #, random_state=12
    X_train = data.loc[(np.array([i.strip().split("-")[1] for i in list(data.index)]) == id1),:]
    X_test = data.loc[(np.array([i.strip().split("-")[1] for i in list(data.index)]) == id2),:]
    #print(sum(np.array([i.strip().split("-")[1] for i in list(data.index)]) == id1))
    #print(sum(np.array([i.strip().split("-")[1] for i in list(data.index)]) == id2))    
    y_train_label = np.array([i[0] for i in X_train.index])
    y_train = np.array([0 if i== control_state else 1 for i in y_train_label])
    y_test_label = np.array([i[0] for i in X_test.index])
    y_test = np.array([0 if i== control_state else 1 for i in y_test_label])
    probas = clf.fit(X_train, y_train).predict_proba(X_test)
    fpr, tpr, thresholds = roc_curve(y_test, probas[:, 1])
    roc_auc = auc(fpr, tpr)
    print(roc_auc)
    return roc_auc

def get_kfold_auc3(data, id1, n_estimators=501, max_features=0.1):
    aucs = []
    tprs = []
    mean_fpr = np.linspace(0, 1, 100)
    plot_data = []
    i = 0
    clf = RandomForestClassifier(n_estimators=n_estimators, oob_score=True, max_features=max_features, random_state=SEED) #, random_state=12
    sub_data = data.loc[(np.array([i.strip().split("-")[1] for i in list(data.index)]) == id1),:]
    sub_label = np.array([i[0] for i in sub_data.index])
    sub_label = label = np.array([0 if i== control_state else 1 for i in sub_label])
    ss=StratifiedKFold(n_splits=5,random_state = SEED)
    for train_index, test_index in ss.split(sub_data, sub_label):
       y_train, y_test = np.array(sub_label)[train_index], np.array(sub_label)[test_index]
       X_train, X_test = sub_data.iloc[train_index,:], sub_data.iloc[test_index,:]
       probas = clf.fit(X_train, y_train).predict_proba(X_test)
       fpr, tpr, thresholds = roc_curve(y_test, probas[:, 1])
       roc_auc = auc(fpr, tpr)
       aucs.append(roc_auc)
    print(np.mean(aucs))
    return np.mean(aucs)

def get_kfold_auc2(data, id1, n_estimators=501, max_features=0.1):
    aucs = []
    tprs = []
    mean_fpr = np.linspace(0, 1, 100)
    plot_data = []
    i = 0
    clf = RandomForestClassifier(n_estimators=n_estimators, oob_score=True, max_features=max_features, random_state=SEED) #, random_state=12
    X_train = data.loc[[i for i in list(data.index) if i.strip().split("-")[1] in (set(study_ids)-set([id1]))],:]
    X_test = data.loc[(np.array([i.strip().split("-")[1] for i in list(data.index)]) == id1),:]
    #print(sum(np.array([i.strip().split("-")[1] for i in list(data.index)]) == id1))
    #print(sum(np.array([i.strip().split("-")[1] for i in list(data.index)]) == id2))
    y_train_label = np.array([i[0] for i in X_train.index])
    y_train = np.array([0 if i== control_state else 1 for i in y_train_label])
    y_test_label = np.array([i[0] for i in X_test.index])
    y_test = np.array([0 if i== control_state else 1 for i in y_test_label])
    probas = clf.fit(X_train, y_train).predict_proba(X_test)
    fpr, tpr, thresholds = roc_curve(y_test, probas[:, 1])
    roc_auc = auc(fpr, tpr)
    print(roc_auc)
    return roc_auc

def dataset_reader(ds='an'): # ds: an, ca ,cn
    study_ids = ['6070', '290926', '389927'] if ds!='ca' else ['6070', '290926', '362366', '389927']
    rm_state = 'C' if ds=='an' else ('A' if ds=='cn' else 'N')
    control_state = 'N' if ds!='ca' else 'A'
    ### read dataset
    data = pd.read_csv('OTU_del_all_shannon2.csv', index_col=0)
    label = np.array([i[0] for i in data.index])
    data_sub = data.loc[label!=rm_state, :]
    data = data_sub
    data = data.iloc[:, 1:]
    ### read features
    SEED, best_auc, best_features, best_plot_data, feature_rank = pickle.load(open('Test_outs_'+ds+'501.pkl', 'rb'))
    data = data.loc[:, best_features]
    return study_ids, control_state, data
