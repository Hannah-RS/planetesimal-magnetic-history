#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plots for dynamo duration across all parameter variations
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

#%% Load a given variable
variables = ['rcmf','eta0','beta','etal','Xs_0','Fe0','alpha_n','r']

for i, var in enumerate(variables):
    unit = units[var]
    varlab = labels[var]
    logvar = logs[i]
    save = False
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
    
    var1=var_data.loc[var_data['run']==minrun,var].values[0]
    var2=var_data.loc[var_data['run']==maxrun,var].values[0]

#%% Make the plot
    if logvar == False:
        width = (data.loc[1,var]-data.loc[0,var])/3 #width is 1/4 gap between variables
    else:
        width = np.diff(data[var])/10
        width = np.append(width,width[-1]*10)
        
    dur1 = data['magoff_1']-data['magon_1']
    dur2 = data['magoff_2']-data['magon_2']
    gap = data['magon_2']-data['magoff_1']
    fig = plt.figure()
    ax = plt.axes()
    ax2 = ax.twinx()
    ln1 = ax.bar(data[var],dur1,width=width,label='First dynamo',color='#4477AA')
    ln2 = ax.bar(data[var],dur2,bottom=dur1,width=width,label='Second dynamo',color='#66CCEE')
    if logvar == True:
        ln3 = ax2.bar((data[var]+width)[gap>0],gap[gap>0],width=width[gap>0],color='#AA3377',label='gap in generation')
    else:
        ln3 = ax2.bar((data[var]+width)[gap>0],gap[gap>0],width=width/2,color='#AA3377',label='gap in generation')
    ax.legend([ln1,ln2,ln3],['First dynamo','Second dynamo','Gap between dynamos'],framealpha=1,bbox_to_anchor=[0.7,-0.2])
    ax.set_ylim([0,max(dur1+dur2)+10])
    ax.set_ylabel('Dynamo duration /Myr')
    ax2.set_ylabel('Gap in dynamo generation /Myr')
    ax.tick_params(axis='y',colors='#4477AA')
    ax2.tick_params(axis='y',colors='#AA3377')
    ax.set_xlabel(f'{varlab}{unit}')
    ax.set_xticks(data[var],data[var])
    
    if logvar == True:
        ax.set_xscale('log')
    
    if save == True:
        plt.savefig(f'../Plots/{folder}/dynamobar_{var}.png',bbox_inches='tight')