#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plots for average field strength and supercritical Rem across all parameter variations
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
split = True
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
    
    if split == True: #split up by dynamo generation periods
        Bav = np.zeros([2,nrun])
        Remav = np.zeros([2,nrun])
    else:
        Bav = np.zeros([nrun])
        Remav = np.zeros([nrun])
    
    for i in range(nrun):
        run = int(minrun+i)
        varval = var_data.loc[var_data['run']==run,var].values[0]
        #import data
        npzfile = np.load(f'{path}run_{run}.npz')
        B = npzfile['B']
        Rem = npzfile['Rem']
        t = npzfile['t']/Myr
        threshold = 10
        #find averages
        if len(Rem[Rem>threshold])>0:
            if split == True:
                on1 =data.loc[data['run']==run,'magon_1'].values[0]
                off1 = data.loc[data['run']==run,'magoff_1'].values[0]
                on2 = data.loc[data['run']==run,'magon_2'].values[0]
                off2 = data.loc[data['run']==run,'magoff_2'].values[0]
                if on1 > 0:
                    Bav[0,i] = np.average(B[(Rem>threshold)&(t>on1)&(t<off1)])
                    Remav[0,i] = np.average(Rem[(Rem>threshold)&(t>on1)&(t<off1)])
                if on2 > 0:
                    Bav[1,i] = np.average(B[(Rem>threshold)&(t>on2)&(t<off2)])
                    Remav[1,i] = np.average(Rem[(Rem>threshold)&(t>on2)&(t<off2)])
            else:
                Bav[i] = np.average(B[Rem>threshold])
                Remav[i] = np.average(Rem[Rem>threshold])
        
#%% Make the plot       
    fig = plt.figure()
    ax1 = plt.axes()
    ax2 = ax1.twinx()
    if split == True:
        ln1 = ax1.scatter(data[var],Bav[0,:]/1e-6,label='B - first',marker='v',color=bcol)
        ln2 = ax2.scatter(data[var],Remav[0,:],label='Rem - first',marker='v',color=rcol)
        ln3 = ax1.scatter(data[Bav[1,:]>0][var],Bav[1,Bav[1,:]>0]/1e-6,label='B - second',marker='^',color='navy')
        ln4 = ax2.scatter(data[Bav[1,:]>0][var],Remav[1,Bav[1,:]>0],label='Rem - second',marker='^',color='limegreen')
        fig.legend(bbox_to_anchor=[0.6,-0.01])
    else:
        ln1 = ax1.scatter(data[var],Bav/1e-6,label='B',marker='v',color=bcol)
        ln2 = ax2.scatter(data[var],Remav,label='Rem',marker='^',color=rcol)
        fig.legend(['B','Rem'],bbox_to_anchor=[0.4,0.8])
    ax1.set_xlabel(f'{varlab}/{unit}')
    if logvar == True:
        ax1.set_xscale('log')
    if var == 'Fe0':
        plt.xticks(data[var],var_data['Fe0'])
    
    ax1.set_ylabel('Average field strength /$\mu$T')
    ax2.set_ylabel('Average supercritical Rem')
    
    #changing colours
    ax1.tick_params(axis='y',colors=bcol)
    ax2.tick_params(axis='y',colors=rcol)
    #ax1.yaxis.label.set_color(bcol) 
    #ax2.yaxis.label.set_color(rcol) 
    #ax2.spines['left'].set_color(bcol) #ax2 plots over ax1 so change colour of both
    #ax2.spines['right'].set_color(rcol)
    if save == True:
        plt.savefig(f'../Plots/{folder}/avB_split_{var}.png',bbox_inches='tight')
        
    