#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot erosion of core stratification and onset of dynamo as a function of 60Fe
Uses data in Paper_runfe which is from extra 60fe runs and the main model runs
"""
#%% load modules
import matplotlib.pyplot as plt
import pandas as pd

#%% import data
data = pd.read_csv('../Results_combined/Paper_runfe/all_sucess_info.csv',skiprows=[1],delimiter=',')
#change fe0= 0 point
data.loc[data['Fe0']==0,'Fe0']=3e-10 #in Inkscape will relabel this point
#%% plot data
plt.figure()
plt.scatter(data['Fe0'],data['terode'],label='Erosion of stratification',marker='x',color='black')
plt.scatter(data['Fe0'],data['magon_1'],label='Dynamo onset',marker='+',color='#A00143',s=100)
plt.xlabel('$^{60}Fe/^{56}Fe$')
plt.ylabel('Time after CAI formation')
plt.legend(fontsize=9,loc='center right')
plt.xscale('log')
plt.xlim([2e-10,1e-6])
plt.fill_between([0,1e-6],3.8,5,alpha=0.3,color='grey')
plt.fill_betweenx([0,12],1.6e-9,color='seagreen',alpha=0.2)
plt.fill_between([0,1e-6],11,11.1,alpha=0.3,color='#A00143')
plt.text(1e-8,4.2,'Volcanic angrite and NWA 7325 ages',size=9)
plt.text(2e-8,11.1,'Angra dos Reis',size=9)
plt.text(5e-10,8,'Possible \nrange of \n$^{60}Fe/^{56}Fe$',size=9)
#plt.vlines(1.75e-9,0,12,linestyle='dashed',color='black')
plt.ylim([1,12])
plt.savefig('../Plots/EPSL_paper/angrite_fe.pdf',dpi=450,bbox_inches='tight')