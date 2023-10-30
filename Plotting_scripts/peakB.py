#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plots for peak field strength and Rem across all parameter variations
"""
import matplotlib.pyplot as plt
import pandas as pd

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
    fig = plt.figure()
    ax1 = plt.axes()
    ax2 = ax1.twinx()
    ln1 = ax1.scatter(data[var],data['max_B']/1e-6,label='B',marker='v',color=bcol)
    ln2 = ax2.scatter(data[var],data['max_R'],label='Rem',marker='^',color=rcol)
    ax1.set_xlabel(varlab)
    if logvar == True:
        ax1.set_xscale('log')
    if var == 'Fe0':
        plt.xticks(data[var],var_data['Fe0'])
    fig.legend(['B','Rem'],bbox_to_anchor=[0.4,0.8])
    ax1.set_ylabel('Peak field strength /$\mu$T')
    ax2.set_ylabel('Peak magnetic Reynolds number')
    
    #changing colours
    ax1.tick_params(axis='y',colors=bcol)
    ax2.tick_params(axis='y',colors=rcol)
    #ax1.yaxis.label.set_color(bcol) 
    #ax2.yaxis.label.set_color(rcol) 
    #ax2.spines['left'].set_color(bcol) #ax2 plots over ax1 so change colour of both
    #ax2.spines['right'].set_color(rcol)
    if save == True:
        plt.savefig(f'../Plots/{folder}/peakB_all_{var}.png',bbox_inches='tight')
        
    