#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plots for peak mantle temperature across all parameter variations
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

#%% Set up folders etc.
folder = 'Paper_run2/'
subfolders = {'rcmf':1,'eta0':2,'beta':3,'etal':4,'Xs_0':5,'Fe0':6,'alpha_n':7,'r':8}
labels = {'rcmf':'$\\phi_{{RCMF}}$','eta0':'$\\eta_0$','beta':'$\\beta$','etal':'$\\eta_l$ ','Xs_0':'$X_{{s,0}}$','Fe0':'$^{{60}}Fe/^{{56}}Fe$','alpha_n':'$\\alpha_n$','r':'radius'}
units = {'rcmf':'','eta0':'Pas','beta':'$K^{-1}$','etal':'Pas','Xs_0':'wt %','Fe0':'','alpha_n':'','r':'km'}
logs =[False,True,False,True,False,True,False,False]
Myr = 365*24*3600*1e6 #number of s in Myr

#plot unchanged B
bcol = 'royalblue'
rcol = 'forestgreen'

save = True

#%% Load a given variable
variables = ['rcmf','eta0','beta','etal','Xs_0','Fe0','alpha_n','r']

for i, var in enumerate(variables):
    unit = units[var]
    varlab = labels[var]
    logvar = logs[i]
    path = '../Results_combined/'+folder+f"params_{subfolders[var]}/"
    
    #find run numbers
    var_data = pd.read_csv(path+'auto_params.csv',skiprows=[1])
    var_results = pd.read_csv(f'../Results_combined/{folder}/all_sucess_info.csv',skiprows=[1])
    minrun = min(var_data['run'])
    maxrun = max(var_data['run'])
    nrun = len(var_data)
    data = var_results[(var_results['run']>=minrun)&(var_results['run']<=maxrun)].copy(deep=True)
    data.reset_index(inplace=True,drop=True)
    
    if var == 'r':
        data[var] = data[var]/1e3 #convert to km
    if var == 'Fe0':
        data.loc[data['Fe0']==0,'Fe0']=1e-10

#%% Make the plot       
    #find upper and lower bounds for time colour maps
    combt = np.concatenate([data['tmax'],data['tcoremax']])
    tmin = min(combt)
    tmax = max(combt)
    #make the figure
    
    plt.figure()
    plt.scatter(data[var],data['peakT'],label='Peak mantle temperature',marker='o',c=data['tmax'],vmin=tmin,vmax=tmax)
    plt.scatter(data[var],data['peak_coreT'],label='Peak core temperature',marker='v',c=data['tcoremax'],vmin=tmin,vmax=tmax)
    plt.xlabel(f'{varlab} {unit}')
    plt.legend()
    plt.ylabel('Temperature /K')
    plt.colorbar(label='Time of maxima/Myr')
    if logvar == True:
        plt.xscale('log')
    if var == 'Fe0':
        plt.xticks(data[var],var_data['Fe0'])
    if save == True:
        plt.savefig(f'../Plots/{folder}/peakT_{var}.png')
        
        
    