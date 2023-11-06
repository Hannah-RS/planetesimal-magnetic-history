#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make stacked timings heatmap for all body sizes
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%% Load data for chosen variable
folders = ['Paper_run100km/','Paper_run200km/','Paper_run4/','Paper_run400km/','Paper_run500km/']
radius = [100,200,300,400,500]
subfolders = {'rcmf':1,'eta0':2,'beta':3,'etal':4,'Xs_0':5,'Fe0':6,'alpha_n':7}
labels = {'rcmf':'$\\phi_{{RCMF}}$','eta0':'$\\eta_0$','beta':'$\\beta$','etal':'$\\eta_l$ ','Xs_0':'$X_{{s,0}}$','Fe0':'$^{{60}}Fe/^{{56}}Fe$','alpha_n':'$\\alpha_n$'}
units = {'rcmf':'','eta0':'Pas','beta':'$K^{-1}$','etal':'Pas','Xs_0':'wt %','Fe0':'','alpha_n':''}
logs =[False,True,False,True,False,True,False]
variables = ['rcmf','eta0','beta','etal','Xs_0','Fe0','alpha_n']
y = np.arange(len(variables)+1)
varlabels = []
Myr = 365*24*3600*1e6 #number of s in Myr
bcol = 'white'
ecol = 'gray'
xmax = 400 #max xlimit in Myr
save = True

for k, folder in enumerate(folders):
    r = radius[k]
#%% Start plot
    plt.figure(figsize=[12,5])
    #%% Load a given variable
    
    for i, var in enumerate(variables):
        unit = units[var]
        if k ==0: #create labels for y axis
            varlab = labels[var]
            varlabels.append(varlab)
        logvar = logs[i]
        path = '../Results_combined/'+folder+f"params_{subfolders[var]}/"
        
        #find run numbers
        var_data = pd.read_csv(path+'auto_params.csv',skiprows=[1])
        var_results = pd.read_csv(f'../Results_combined/{folder}/all_sucess_info.csv',skiprows=[1])
        minrun = min(var_data['run'])
        maxrun = max(var_data['run'])
        var1=var_data.loc[var_data['run']==minrun,var].values[0]
        var2=var_data.loc[var_data['run']==maxrun,var].values[0]
        min_start = var_data.loc[var_data['run']==minrun,'t_acc_m'].values[0] #min start
        max_end = max(var_results.loc[var_results['r']==r*1e3,'tsolid']) #max end
        data = var_results[(var_results['run']>=minrun)&(var_results['run']<=maxrun)].copy(deep=True)
        data.reset_index(inplace=True,drop=True)
        nrun = len(data)
        
        #%% Create time array
        tdyn = np.concatenate([data[data['magon_1']>0]['magon_1'],data[data['magon_1']>0]['magoff_1'],data[data['magon_2']>0]['magon_2'],data[data['magon_2']>0]['magoff_2']])
        tdyn = np.sort(tdyn) #range from low to high
        tmid = (tdyn[:-1]+tdyn[1:])/2 #find midpoints
        ton = np.concatenate([data['magon_1'],data['magon_2']]) #on times
        toff = np.concatenate([data['magoff_1'],data['magoff_2']]) #off times
        weight=np.zeros([len(tmid),1])
        
        for j, t in enumerate(tmid):
            weight[j] = 100*(len(ton[ton<t])-len(toff[toff<t]))/nrun
        
            
        #%% Make the plot
        if len(tdyn)>0: #check the dynamo is on at all
            plt.pcolormesh(tdyn,y[i:i+2],np.transpose(weight),shading='flat',vmin=0,vmax=100)
        
    #%% Customise the plot
    plt.colorbar(label='% dynamos on')
    plt.yticks((y[:-1]+y[1:])/2,varlabels)
    plt.ylabel('Variable')
    plt.xlabel('Time/Myr')
    if r==100: #force all variables to appear
        plt.ylim([0,y[-1]])
    plt.xlim([0.8,xmax])
    if xmax > 200:
        plt.title(f'Chance of a dynamo being on across parameters for a {r}km body')
    else:
        plt.title(f'Onset times for a {r}km body')
    #plt.text(210,3.5,'Core solidified',bbox=dict(edgecolor=ecol,facecolor=bcol))
    #plt.text(1.4,3,'Core \n thermally \n stratified',bbox=dict(edgecolor=ecol,facecolor=bcol))
    if save == True:
        if xmax > 200:
            plt.savefig(f'../Plots/CoS/timing_heatmap_{r}.png',dpi=450,bbox_inches='tight')
        else:
            plt.savefig(f'../Plots/CoS/onset_heatmap_{r}.png',dpi=450,bbox_inches='tight')