#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 00:12:38 2020

@author: wuyuanqi
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

## control versus adenoma
data = pd.read_csv('data/Figure3c.csv', index_col=0)

fig = plt.figure(figsize=(6, 5), dpi=300)

font1 = {'family' : 'Arial','weight' : 'normal','size': 14}
font2 = {'family' : 'Arial','weight' : 'normal','size': 18}

plt.ylabel('Average AUC',font1)
plt.xlabel("Used features",font1)
plt.yticks(fontproperties = 'Arial', size = 13)
plt.xticks(fontproperties = 'Arial', size = 13)

plt.plot(list(reversed(data.index)), list(reversed(data.iloc[:,0])), color='#EF9B9E', marker='o',label='FR',linestyle='dashed')
plt.plot(list(reversed(data.index)), list(reversed(data.iloc[:,1])), color='#5DA0ED', marker='o',label='US2',linestyle='dashed')
plt.plot(list(reversed(data.index)), list(reversed(data.iloc[:,2])), color='#37A2A9', marker='o',label='CA',linestyle='dashed')


plt.ylim(0.40,0.75)##AN_study
#plt.ylim(0.40,0.90)##AN_LODO
plt.legend(loc='lower left',prop=font1) # 
fig.tight_layout()
plt.show()
fig.savefig('figure/Figure3c.pdf',dpi=300,format='pdf')


## adenoma versus cancer
data = pd.read_csv('data/Figure3d.csv', index_col=0)

fig = plt.figure(figsize=(6.5, 5), dpi=300)

font1 = {'family' : 'Arial','weight' : 'normal','size': 14}
font2 = {'family' : 'Arial','weight' : 'normal','size': 18}

plt.ylabel('Average AUC',font1)
plt.xlabel("Used features",font1)
plt.xticks(fontproperties = 'Arial',size=13)
plt.yticks(fontproperties = 'Arial',size=13)

plt.plot(list(reversed(data.index)), list(reversed(data.iloc[:,0])), color='#EF9B9E', marker='o',label='FR',linestyle='dashed')
plt.plot(list(reversed(data.index)), list(reversed(data.iloc[:,1])), color='#5DA0ED', marker='o',label='US2',linestyle='dashed')
plt.plot(list(reversed(data.index)), list(reversed(data.iloc[:,2])), color='#E2BA6B', marker='o',label='US1',linestyle='dashed')
plt.plot(list(reversed(data.index)), list(reversed(data.iloc[:,3])), color='#37A2A9', marker='o',label='CA',linestyle='dashed')

plt.ylim(0.48,0.85)##CA_study
#plt.ylim(0.4,1.00)##CA_LODO
plt.legend(loc='lower right',prop=font1)
fig.tight_layout()
plt.show()
fig.savefig('figure/Figure3d.pdf',dpi=300,format='pdf')

