#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 21:44:39 2020

@author: wuyuanqi
"""
import numpy as np
import pandas as pd
import random
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib import axes
import seaborn as sns

####calculate intra validation results

plot_data = []
for id1 in study_ids:
    temp = []
    for id2 in study_ids:
        try:
            a = get_kfold_auc(data, id1, id2, n_estimators=501, max_features=0.1)
            temp.append(round(a,2))
        except:
            pass
    plot_data.append(temp)


plot_diagonal = []
for id1 in study_ids:
    k = get_kfold_auc3(data, id1, n_estimators=501, max_features=0.1)
    plot_diagonal.append(round(k,2))


LODO = []
for id1 in study_ids:
    try:
        n = get_kfold_auc2(data, id1, n_estimators=501, max_features=0.1)
    except:
        pass  
    LODO.append(round(n,2))

LODO.append(np.mean(LODO))

#LODO = np.array(LODO).reshape(1,4)
LODO = np.array(LODO).reshape(1,5)

######prepare plot_heatmap_data
matrix = np.array(plot_data)
main = np.array(plot_diagonal)
np.fill_diagonal(matrix, main)#study-to-study4*4 matrix
matrix = np.c_[matrix,np.mean(matrix, axis=1)]
#a = np.mean(matrix, axis =0).reshape(1,4)
a = np.mean(matrix, axis =0).reshape(1,5)
plot_matrix = np.around(np.r_[matrix,a],2)
#plot_matrix[3,3] = np.mean(plot_matrix[0:3,0:3])
plot_matrix[4,4] = np.mean(plot_matrix[0:4,0:4])
#print(plot_matrix)
plot_matrix = np.around(np.r_[plot_matrix,LODO],2)

####plot
#
xLabel = ['FR', 'US2', 'US1', 'CA', 'Average']
yLabel = ['FR', 'US2', 'US1', 'CA', 'Average','LODO']

fig = plt.figure(1, (6, 5.5), dpi=300)
font1 = {'family' : 'Arial','weight' : 'normal','size': 12}
font2 = {'family' : 'Arial','weight' : 'normal','size': 16}
##plot heatmap

ax = sns.heatmap(plot_matrix,vmin=0.5, vmax=1, cmap = 'YlGnBu_r', linewidths= 0.8, annot=True,annot_kws={'size':15},
                 center=0.8, xticklabels = xLabel, yticklabels = yLabel)  
#y tick
ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

#set lable
label_x = ax.get_xticklabels()
plt.setp(label_x, rotation= -30, horizontalalignment='right')
label_y = ax.get_yticklabels()
plt.setp(label_y, rotation= 0)
#add title
plt.title("Adenoma vs Cancer Cross-prediction",font2)
#plt.title("Control vs Adenoma Cross-prediction",font2)
plt.ylabel('Training Cohorts',font2)
plt.xlabel('Testing Cohorts',font2)
plt.yticks(fontproperties = 'Arial', size = 13)
plt.xticks(fontproperties = 'Arial', size = 13)

fig.tight_layout()
plt.show()
fig.savefig('figure/Figure3b.pdf',dpi=300,format='pdf')
