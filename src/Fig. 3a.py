#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 15:51:15 2020

@author: wuyuanqi
"""


import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.lines as mlines
from matplotlib.font_manager import FontProperties
import seaborn as sns
from scipy.stats import norm, pearsonr, spearmanr
import scipy.stats as stats
from scipy.spatial import distance
import cloudpickle as pickle

from sklearn import svm, tree, linear_model, neighbors, naive_bayes, ensemble, discriminant_analysis, gaussian_process
#Common Model Helpers
from sklearn.preprocessing import OneHotEncoder, LabelEncoder
from sklearn import feature_selection
from sklearn import model_selection
from sklearn import metrics
from sklearn import preprocessing
import warnings
warnings.filterwarnings("ignore")

### 设置模型优化的基础模型和超参数池
def set_tune_params(max_features=[1.0], max_samples=[1.0], cpu=2):
    tune_params = {
        'Bagging_kn':['ensemble.BaggingClassifier()', 'roc_auc', 
                        {'base_estimator':[neighbors.KNeighborsClassifier(algorithm='auto', metric='braycurtis', n_neighbors=3, weights='distance'),],
                         'n_estimators': [501],
                         'max_features':max_features,
                         'max_samples':max_samples,
                         'bootstrap_features':[True], # 随机置换feature
                         'bootstrap':[True], # 随机置换sample
                         'oob_score': [True],
                         'random_state': [0],
                         'n_jobs':[cpu], # CPU核数，可在大机器上添加核数加速运算
                        }],

        'RandomForest':['ensemble.RandomForestClassifier()', 'roc_auc', 
                        {'n_estimators': [501], #default=10
                         'criterion': ['gini'],# entropy
                         'max_features':max_features,
                         'max_samples':max_samples,
                         #'min_samples_leaf':[1, 2, 3],
                         'max_depth': [1, 2, 3], # 防止过拟合，减少max_depth
                         'oob_score': [True],
                         'random_state': [0],
                         'n_jobs':[cpu],
                        }],
        }
    return tune_params

### 模型优化，选择最优超参数
def tune_model(X, y, cv_split, model, param_grid, scoring='roc_auc'):
    #basic model 基础模型训练（默认超参数）
    basic_model = eval(model)
    basic_results = model_selection.cross_validate(basic_model, X, y, cv=cv_split, scoring = scoring, return_train_score=True)
    #tune model 模型优化 (超参数网格遍历)
    tune_model = model_selection.GridSearchCV(eval(model), param_grid=param_grid, 
                                              scoring = scoring, cv=cv_split, return_train_score=True)
    _ = tune_model.fit(X, y)
    ### 获取最优超参数，并建模
    best_param = tune_model.best_params_
    final_model = eval(model).set_params(**best_param)
    final_results = model_selection.cross_validate(final_model, X, y, cv=cv_split, scoring = scoring, return_train_score=True)
    return [final_model, best_param, final_results['test_score'],
            basic_results['train_score'].mean(), basic_results['test_score'].mean(), 
            tune_model.cv_results_['mean_train_score'][tune_model.best_index_], 
            tune_model.cv_results_['mean_test_score'][tune_model.best_index_]]

### 导入数据
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

### self建模
def model_self(data, study, control_state, model, scoring, param_grid):
    index = np.array([i.split('-')[1] for i in data.index])==study
    X = data.loc[index, :].values
    y = np.array([0 if i[0]==control_state else 1 for i in data.loc[index, :].index])
    nor = preprocessing.MinMaxScaler()
    X[:, -3:] = nor.fit_transform(X[:, -3:])
    ### cross validate 设置
    cv_split = list(model_selection.StratifiedKFold(n_splits=5, random_state = RANDOM_SEED).split(X, y))
    ### 模型优化
    tune_results = tune_model(X, y, cv_split, model, param_grid, scoring)
    return tune_results, 0.0

### Study-study
def model_cross_study(model, data, study_train, study_test, control_state):
    ### Train
    train_index = np.array([i.split('-')[1] for i in data.index])==study_train
    X_train = data.loc[train_index, :].values
    y_train = np.array([0 if i[0]==control_state else 1 for i in data.loc[train_index, :].index])
    ### Test
    test_index = np.array([i.split('-')[1] for i in data.index])==study_test
    X_test = data.loc[test_index, :].values
    y_test = np.array([0 if i[0]==control_state else 1 for i in data.loc[test_index, :].index])
    nor = preprocessing.MinMaxScaler()
    X_train[:, -3:] = nor.fit_transform(X_train[:, -3:])
    X_test[:, -3:] = nor.transform(X_test[:, -3:])
    # 使用最优模型重新训练全部训练数据，并在测试数据中验证
    model.fit(X_train, y_train)
    probas = model.predict_proba(X_test)
    fpr, tpr, thresholds = metrics.roc_curve(y_test, probas[:, 1])
    score = metrics.auc(fpr, tpr)
    return score

### LODO
def model_LODO(data, study_ids, study, control_state, model, scoring, param_grid, cv_per_study=5, cv_ratio=0.8):
    train_index = np.array([i.split('-')[1] for i in data.index])!=study
    X_train = data.loc[train_index, :].values
    y_train = np.array([0 if i[0]==control_state else 1 for i in data.loc[train_index, :].index])
    test_index = np.array([i.split('-')[1] for i in data.index])==study
    X_test = data.loc[test_index, :].values
    y_test = np.array([0 if i[0]==control_state else 1 for i in data.loc[test_index, :].index])
    nor = preprocessing.MinMaxScaler()
    X_train[:, -3:] = nor.fit_transform(X_train[:, -3:])
    X_test[:, -3:] = nor.transform(X_test[:, -3:]) 
    # 使用训练数据进行cross-validation, 并挑选最优模型
    ### cross validate 设置 【原来的训练数据分为train和valid】
    ### Cross-validation between studies 【在LODO中，为了使模型捕获study之间的批次效应，验证设置采用不同stud数据，即训练和验证从不同study中抽取】
    train_ids = data.index[train_index] # cross-validation只在训练数据中进行
    cv_split = []
    for valid_study in set(study_ids)-set([study]): # 设置验证study
        train_index = np.arange(len(train_ids))[np.array([i.split('-')[1] for i in train_ids])!=valid_study] # 训练数据来源
        valid_index = np.arange(len(train_ids))[np.array([i.split('-')[1] for i in train_ids])==valid_study] # 验证数据来源
        for rt in range(cv_per_study): # 每个study抽取N次随机样本，构建训练验证数据
            cv_split.append([np.random.choice(train_index, int(len(train_index)*cv_ratio)),  # 每次仅使用0.8的样本
                             np.random.choice(valid_index, int(len(valid_index)*cv_ratio))])
    #cv_split = model_selection.StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_SEED) # 随机的cv方式，结果差一些
    ### 模型优化
    tune_results = tune_model(X_train, y_train, cv_split, model, param_grid, scoring)
    # 使用最优模型重新训练全部训练数据，并在测试数据中验证
    probas = tune_results[0].fit(X_train, y_train).predict_proba(X_test)
    fpr, tpr, thresholds = metrics.roc_curve(y_test, probas[:, 1])
    score = metrics.auc(fpr, tpr)
    return tune_results, score

### main function
def model_one_dataset(ds, methods, max_features, max_samples, cpu, outfile):
    study_ids, control_state, data = dataset_reader(ds)
    tune_params = set_tune_params(max_features=max_features, max_samples=max_samples, cpu=cpu)
    # self and study to study
    print("$$$$$$$$$$$$ Start model of study to study $$$$$$$$$$$$")
    outfile.write("### Start model of study to study...\n")
    for study_i in study_ids:
        for model_name in methods:
            model, scoring, param_grid = tune_params[model_name]
            # 首先5fold cross-validation测试自身效果，并确定最优模型超参数
            [final_model, params, valid_scores, 
             basic_train_score, basic_valid_score, 
             tune_train_score, tune_valid_score], test_score = model_self(data, study_i, control_state, model, scoring, param_grid)
            # 使用最优模型测试其他study数据
            scores = []
            for study_j in study_ids:
                score = model_cross_study(final_model, data, study_i, study_j, control_state)
                scores.append([study_j, score])
            print("### Train study:{}, Model:{}, Basic model[Train:{:.3f}, Valid:{:.3f}], Tune model[Train:{:.3f}, Valild:{:.3f}], Test:{:.3f}". format(study_i, model_name, basic_train_score, basic_valid_score, tune_train_score, tune_valid_score, test_score))
            outfile.write("### Train study:{}, Model:{}, Basic model[Train:{:.3f}, Valid:{:.3f}], Tune model[Train:{:.3f}, Valild:{:.3f}], Test:{:.3f}". format(study_i, model_name, basic_train_score, basic_valid_score, tune_train_score, tune_valid_score, test_score))
            print('Valid Scores:', valid_scores)
            print('Study-Study Scores:', scores)
            outfile.write('Valid Scores:'+str(valid_scores)+'\n')
            outfile.write('Study-Study Scores:'+str(scores)+'\n')
            print(params)
            outfile.write(str(params)+'\n')
            outfile.write('\n')
            outfile.flush()
        print('$$$$$$$$$$$$\n\n')
    
    # LODO
    print("$$$$$$$$$$$$ Start model of LODO $$$$$$$$$$$$")
    for study in study_ids:
        for model_name in methods:
            model, scoring, param_grid = tune_params[model_name]
            # 在训练数据测试模型，确定最优模型超参数，最后在测试数据中测试
            [final_model, params, valid_scores, 
             basic_train_score, basic_valid_score, 
             tune_train_score, tune_valid_score], test_score = model_LODO(data, study_ids, study, control_state, model, scoring, param_grid)
            print("### Test study:{}, Model:{}, Basic model[Train:{:.3f}, Valid:{:.3f}], Tune model[Train:{:.3f}, Valild:{:.3f}], Test:{:.3f}". format(study, model_name, basic_train_score, basic_valid_score, tune_train_score, tune_valid_score, test_score))
            outfile.write("### Test study:{}, Model:{}, Basic model[Train:{:.3f}, Valid:{:.3f}], Tune model[Train:{:.3f}, Valild:{:.3f}], Test:{:.3f}". format(study, model_name, basic_train_score, basic_valid_score, tune_train_score, tune_valid_score, test_score))
            print('Valid Scores:', valid_scores)
            outfile.write('Valid Scores:'+str(valid_scores)+'\n')
            print(params)
            outfile.write(str(params)+'\n')
            print('')
            outfile.write('\n\n')
            outfile.flush()
        print('$$$$$$$$$$$$\n\n')

if __name__ == '__main__':
    pass
    
