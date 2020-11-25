#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 20:27:24 2020

@author: wuyuanqi
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cloudpickle as pickle

## plot ROC
SEED, best_auc, best_features, best_plot_data, feature_rank = pickle.load(open('data/Test_asvs_'+ds+'501.pkl', 'rb'))
lines, mean_fpr, mean_tpr, tprs, std_auc = best_plot_data
font1 = {'family' : 'Arial','weight' : 'normal','size': 12}
font2 = {'family' : 'Arial','weight' : 'normal','size': 16}

#plt.figure(1, (6, 5.5), dpi=300)
fig, ax = plt.subplots(figsize=(6,5.5))
for fpr, tpr, label in lines:
    plt.plot(fpr, tpr, lw=1, alpha=0.3)
plt.plot([0, 1], [0, 1], ls='--', lw=2, color='r')
plt.plot(mean_fpr, mean_tpr, color='b', lw=2.5,label='Mean ROC(AUC=%0.2f $\pm$ %0.2f)'%(best_auc, std_auc))


std_tpr = np.std(tprs, axis=0)
tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                label=r'$\pm$ 1 std. dev.')

plt.legend(loc='lower right',prop=font1)
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate', font2)
plt.ylabel('True Positive Rate', font2)
plt.yticks(fontproperties = 'Arial', size = 13)
plt.xticks(fontproperties = 'Arial', size = 13)
plt.show()
fig.savefig('figure/Figure2c.pdf',dpi=300,format='pdf')

