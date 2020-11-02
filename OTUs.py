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
SEED=2

 
#建模
def get_kfold_auc(data, meta_group, label, max_features=0.1, n_estimators=501, k=10):#, min_samples_leaf=80
    aucs = []
    tprs = []
    mean_fpr = np.linspace(0, 1, 100)
    plot_data = []
    i = 0
    splitor = StratifiedKFold(n_splits=k, shuffle=True,random_state=SEED) 
#    splitor = StratifiedShuffleSplit(n_splits=k,train_size=0.7,random_state=SEED)
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



#def get_2fold_auc(data, meta_group, label, n_estimators=201, k=2, max_features=0.1):
#    splitor = StratifiedShuffleSplit(n_splits=2, train_size=0.7)
#    clf = RandomForestClassifier(n_estimators=n_estimators, oob_score=True, max_features=max_features)
#    for train_index, test_index in splitor.split(data, meta_group):
#        y_train, y_test = label[train_index], label[test_index]
#        X_train, X_test = np.array(data)[train_index], np.array(data)[test_index]
#        probas = clf.fit(X_train, y_train).predict_proba(X_test)
#        fpr, tpr, thresholds = roc_curve(y_test, probas[:, 1])
#        roc_auc = auc(fpr, tpr)
#        return roc_auc, clf

##########################
data = pd.read_csv('../new_run/OTUs/new/OTU_del_all_shannon2.csv', index_col=0)
label = np.array([i[0] for i in data.index])
data_sub = data.loc[label!='C', :]
label_sub = label[label!='C']
data = data_sub
label = label_sub
label = np.array([0 if i=='N' else 1 for i in label])
meta_group = np.array(data['Label'])
data = data.iloc[:, 1:]

#样本选择
#meta_data = pd.read_csv('../new_run/metadata_4study.csv', index_col=1)
#sub_label = []
#for i in data.index:
#    sample_id = i.split('-')[1]+'-'+i.split('-')[0]
#    sub_label.append(meta_data.loc[sample_id, 'Group_2'])
#
#for i in set(sub_label):
#    print(i, sub_label.count(i))
##select_sample = [True if i in ['Large Adenoma', 'Small Adenoma', 'adv_Adenoma', 'Adenoma', 'Normal'] else False for i in sub_label]#select1
##select_sample = [True if i in ['Large Adenoma', 'adv_Adenoma', 'Adenoma', 'Normal','High Risk Normal'] else False for i in sub_label]#select2
##select_sample = [True if i in ['Large Adenoma', 'adv_Adenoma', 'Adenoma', 'Normal'] else False for i in sub_label]
#select_sample = [True if i in ['Large Adenoma', 'Adenoma', 'Normal'] else False for i in sub_label]
##select_sample = [True if i in ['adv_Adenoma', 'Adenoma', 'Normal'] else False for i in sub_label]
##select_sample = [True if i in ['Adenoma', 'Normal'] else False for i in sub_label]
##
#data = data.iloc[select_sample, :]
#sub_label = np.array(sub_label)[select_sample]
#label = label[select_sample]
#meta_group = meta_group[select_sample]
################### wilcoxon 差异菌种选择
diff = pd.read_csv('../new_run/OTUs/new/AN/p.val_an_block_otus.csv', index_col=0)
select_features = list(diff.loc[diff['all']<0.05, :].index)
print('### Wilcoxon test :', len(select_features))
print('StratifiedKFlod,k=10,max_features=0.1,seed=2, n_eatimates =501, cn_cafeature')

meta_cols = ['Age', 'Gender', 'BMI', 
             'Shannon', 'Simpson', 'InSimpson', 'spe_number']
select_features.extend(meta_cols)
data = data.loc[:, select_features]
###########
###########################
##遍历
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


    
#pickle.dump([SEED, best_auc, best_features, best_plot_data, feature_rank], open('../new_run/OTUs/new/CN/Test_outs_cn_anfeature.pkl', 'wb'))
#SEED, best_auc, best_features, best_plot_data, feature_rank = pickle.load(open('../new_run/OTUs/new/CA/Test_outs_ca501.pkl', 'rb'))


###########################
#
#import numpy as np
#import pandas as pd
#tax_df = pd.read_csv('../PRJNA290926/PRJNA290926-taxonomy.tsv',sep='\t',index_col =0)
#tax_species = []
#name = []
#for i in tax_df.loc[:,'Taxon']:
#    if 'D_6__' not in i  :
#        if 'D_5__' not in i:
#            if 'D_4__' in i:
#                name.append(i.split('D_4__')[1])
#            else:
#                name.append('NA')
#        elif 'D_5__' in i:
#            name.append(i.split('D_5__')[1])
#    elif 'D_6__' in i:
#        if (i.split('D_6__')[1] == 'uncultured bacterium') or (i.split('D_6__')[1] == 'uncultured organism') or (i.split('D_6__')[1] == 'unidentified'):
#            name.append(i.split('__')[6].split(';')[0]+' sp.')
#        else:
#            name.append(i.split('D_6__')[1])
#           
#            
#        
#tax_df['species'] = name
#
###########################
# 画图
### plot AUC
#lines, mean_fpr, mean_tpr, tprs, std_auc = best_plot_data
#plt.figure(1, (6, 5.5), dpi=100)
#for fpr, tpr, label in lines:
#    plt.plot(fpr, tpr, lw=1, alpha=0.3, label=label)
#plt.plot([0, 1], [0, 1], ls='--', lw=2, color='r')
#plt.plot(mean_fpr, mean_tpr, color='b', label='Mean ROC(AUC=%0.2f $\pm$ %0.2f)'%(best_auc, std_auc))
#plt.legend(loc='lower right',fontsize= 8)
#plt.xlim([-0.05, 1.05])
#plt.ylim([-0.05, 1.05])
#plt.xlabel('False Positive Rate', fontsize=13)
#plt.ylabel('True Positive Rate', fontsize=13)
#plt.show()

## 画图 （有阴影）
## plot AUC

### #############################################################################lines, mean_fpr, mean_tpr, tprs, std_auc = best_plot_data
#lines, mean_fpr, mean_tpr, tprs, std_auc = best_plot_data
#font1 = {'family' : 'Arial','weight' : 'normal','size': 12}
#font2 = {'family' : 'Arial','weight' : 'normal','size': 16}
#
##plt.figure(1, (6, 5.5), dpi=300)
#fig, ax = plt.subplots(figsize=(6,5.5))
#for fpr, tpr, label in lines:
#    plt.plot(fpr, tpr, lw=1, alpha=0.3)
#plt.plot([0, 1], [0, 1], ls='--', lw=2, color='r')
#plt.plot(mean_fpr, mean_tpr, color='b', lw=2.5,label='Mean ROC(AUC=%0.2f $\pm$ %0.2f)'%(best_auc, std_auc))
#
### 画阴影
#std_tpr = np.std(tprs, axis=0)
#tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
#tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
#plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
#                label=r'$\pm$ 1 std. dev.')
#
####s= [0.68,0.81,0.78,0.89,0.79,0.82,0.84,0.78,0.69,0.88] ## an
####s= [0.87,0.93,0.85,0.92,0.89,0.88,0.91,0.92,0.90,0.85] ## ca
##s= [0.86,0.97,0.97,0.92,0.93,0.92,0.94,0.95,0.89,0.98] ## cn
##
##interval=stats.norm.interval(0.95,loc = np.mean(s),scale = st.sem(s))
##interval
##plt.plot([0.2,0.1],[0.4,0.1], color = 'w',label = '95%% CI: %0.2f-%0.2f'%(interval[0], interval[1]))
#
#plt.legend(loc='lower right',prop=font1)
#plt.xlim([-0.05, 1.05])
#plt.ylim([-0.05, 1.05])
#plt.xlabel('False Positive Rate', font2)
#plt.ylabel('True Positive Rate', font2)
#plt.yticks(fontproperties = 'Arial', size = 13)
#plt.xticks(fontproperties = 'Arial', size = 13)
#plt.show()
#fig.savefig('../new_run/OTUs/new/RF_501_ca_11.1.pdf',dpi=300,format='pdf')
#
#### plot Feature
#def autolabel(rects, value):
#    i = 0
#    for rect in rects:
#        width = rect.get_width()
#        left = value[i]
#        plt.text(1.01*(width+left), rect.get_y()+rect.get_height()/2.-0.3, '%s' % round(float(width), 3))
#        i+=1
#
##def get_species_name(feature):
##    if 'D_5__' not in feature:
##        return 'NA'
##    name = feature.split('D_5__')[1]
##    name = name.split(';')[0]
##    return name
#
#def get_species_name(feature):
#    if (feature == 'Age') or (feature == 'BMI') or (feature == 'Gender'):
#        return feature
#    else:
#        return tax_df.loc[feature,'species']
#
#
#def plot_feature_importance(feature_rank, top=30):
#    feature_rank = sorted(feature_rank, key=lambda x:len(x[0]))
#    current = set(['Intercept'])
#    current_auc = 0.5
#    result = [[0, 'Intercept', current_auc, 0.5]]
#    index = 1
#    for sel_features, aucscore in feature_rank[:top]:
#        aucscore = round(aucscore, 3)
#        add_features = list(set(sel_features)-current)[0]
#        #add_auc = aucscore - current_auc
##        add_feature_name = add_features if (add_features == 'BMI' or 'Age') else tax_df.loc[add_features,'species']
##        add_feature_name = add_features if 'D_0' not in add_features else get_species_name(add_features)
#        add_feature_name = get_species_name(add_features)
#        result.append([index, add_feature_name, current_auc, aucscore])
#        current = set(sel_features)
#        current_auc = aucscore
#        #print(index, add_features, add_auc)
#        index += 1
#    result = pd.DataFrame(result)
#    return result
############################
#
############################
##Top = 11
#Top = 30
#result = plot_feature_importance(feature_rank, top=Top)
#plt.figure(figsize=(6, 9), dpi=100)
#color = ['#48d1cc' if (result.loc[i, 3]-result.loc[i, 2])>0 else '#ffdead' for i in result.index]
#rects = plt.barh(range(len(result)), width=result.loc[:, 3]-result.loc[:, 2], left=result.loc[:, 2],
#         align='center', height=1, color=color, linewidth=1, edgecolor='black')
#autolabel(rects, list(result.loc[:, 2]))
#plt.yticks(range(len(result)), result.loc[:, 1])
#plt.xlim([0.5, best_auc+0.1])
#plt.ylim([-0.5, min([len(feature_rank), Top])+0.5])
#plt.ylabel('Features',fontsize=13)
#plt.xlabel('AUC of 10-fold cross validation',fontsize=13)
#plt.show()

############################
#
#
#
#
#
