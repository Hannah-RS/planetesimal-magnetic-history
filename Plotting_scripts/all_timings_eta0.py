#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make stacked timings heatmap for all body sizes for eta0 grouped by above or below gap in dynamo generation
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

#%% Load data for chosen variable
folders = ['Paper_run100km/','Paper_run200km/','Paper_run300km/','Paper_run400km/','Paper_run500km/']

from plot_params import subfolders, labels, units, logs
variables = ['eta0']
y = np.arange(len(variables)+1)
radius = [100,200,300,400,500]
varlabels = []
bcol = 'white'
ecol = 'gray'
xmax = 650 #max xlimit in Myr
save = True

#make colormap
cmap = mpl.cm.viridis
bounds = [0,100]
norm = mpl.colors.Normalize(vmin=0, vmax=100)


fig, ax = plt.subplots(5,1,sharey='row',sharex='col',tight_layout=True,figsize=[5,5])


for k, folder in enumerate(folders):
    r = radius[k]
#%% Start plot
    
    #%% Load a given variable
    for i, var in enumerate(variables):
        unit = units[var]
        if k ==0: #create labels for y axis
            varlab = labels[var]
            varlabels.append(varlab)
        logvar = logs[var]
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
        tdyn = np.concatenate([data[data['magon_1']>0]['magon_1'],data[data['magon_1']>0]['magoff_1'],data[data['magon_2']>0]['magon_2'],data[data['magon_2']>0]['magoff_2'],data[data['magon_3']>0]['magon_3'],data[data['magon_3']>0]['magoff_3']])
        tdyn = np.sort(tdyn) #range from low to high
        tmid = (tdyn[:-1]+tdyn[1:])/2 #find midpoints
        ton = np.concatenate([data['magon_1'],data['magon_2'],data['magon_3']]) #on times
        toff = np.concatenate([data['magoff_1'],data['magoff_2'],data['magoff_3']]) #off times
        weight=np.zeros([len(tmid),1])
        
        for j, t in enumerate(tmid):
            weight[j] = 100*(len(ton[ton<t])-len(toff[toff<t]))/nrun
        
            
        #%% Make the plot
        if len(tdyn)>0: #check the dynamo is on at all
            ax[k].pcolormesh(tdyn,y[0:2],np.transpose(weight),shading='flat',vmin=0,vmax=100)
        
    #%% Customise the plot
    #ax[k].set_yticks([],[])
    ax[k].set_yticks([y[1]/2],['$\\eta_0$'])
    
    
    
    if r==100: #force all variables to appear
        ax[k].set_ylim([0,y[-1]])
    ax[k].set_xlim([0.8,xmax])
    ax[k].set_title(f' Planetesimal radius = {r}km')

ax[4].set_xlabel('Time after CAI formation/Ma')

cax = fig.add_axes([1, 0.36, 0.025, 0.3])
fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm),cax=cax, orientation='vertical', label='% dynamos on')

if save == True:
    plt.savefig(f'../Plots/EPSL_paper/timings_heatmap_eta0.pdf',dpi=500,bbox_inches='tight')