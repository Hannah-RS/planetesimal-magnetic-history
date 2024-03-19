#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Heatmap of dynamo on times as a function of variable
The bottom graph without 60Fe needs fixing as does the creation of the array of desired tick locations
Also can't figure out how to add a colorbar with the same colour scale'
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%% Load data for chosen variable
folder = 'Paper_run2/'
subfolders = {'rcmf':1,'eta0':2,'beta':3,'etal':4,'Xs_0':5,'Fe0':6,'alpha_n':7,'r':8}
labels = {'rcmf':'$\\phi_{{RCMF}}$','eta0':'$\\eta_0$','beta':'$\\beta$','etal':'$\\eta_l$ ','Xs_0':'$X_{{s,0}}$','Fe0':'$^{{60}}Fe/^{{56}}Fe$','alpha_n':'$\\alpha_n$','r':'radius'}
units = {'rcmf':'','eta0':'Pas','beta':'$K^{-1}$','etal':'Pas','Xs_0':'wt %','Fe0':'','alpha_n':'','r':'km'}
logs =[False,True,False,True,False,True,False,False]
variables = ['rcmf','eta0','beta','etal','Xs_0','Fe0','alpha_n','r']
y = np.arange(len(variables)+1)
varlabels = []
ytick_lab = []
ytick_val = np.zeros([len(variables)*2])
Myr = 365*24*3600*1e6 #number of s in Myr

#yplot = np.arange(43)*2
#plot unchanged B
bcol = 'white'
ecol = 'gray'

save = False

#%% Start plot
plt.figure(figsize=[12,5])
#%% Load a given variable

for i, var in enumerate(variables):

    unit = units[var]
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
    data = var_results[(var_results['run']>=minrun)&(var_results['run']<=maxrun)].copy(deep=True)
    data.reset_index(inplace=True,drop=True)
    nrun = len(data)
    if logvar == True:
        ytick_lab.append(f'{varlab}={var1:.1e}{unit}')
        ytick_lab.append(f'{varlab}={var2:.1e}{unit}')
    else:
        ytick_lab.append(f'{varlab}={var1:.1g}{unit}')
        ytick_lab.append(f'{varlab}={var2:.1g}{unit}')
    #this doesn't work because runs 29-31 are missing
    # if i==0:
    #     ytick_val[2*i] = 0
    #     ytick_val[2*i+1] = ytick_val[2*i]+nrun-1
    # elif i==5:
    #     ytick_val[2*i] = ytick_val[2*i-1]+1
    #     ytick_val[2*i+1] = ytick_val[2*i]+nrun-1 
    # else:
    #     ytick_val[2*i] = ytick_val[2*i-1]+1 
    #     ytick_val[2*i+1] = ytick_val[2*i]+nrun-1   
 
    
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
    
    plt.pcolormesh(tdyn,y[i:i+2],np.transpose(weight),shading='flat',vmin=0,vmax=100)

#%% Customise the plot
plt.colorbar(label='% dynamos on')
plt.yticks((y[:-1]+y[1:])/2,varlabels)
plt.xscale('log')
plt.ylabel('Variable')
plt.xlabel('Time/Myr')
plt.text(210,3.5,'Core solidified',bbox=dict(edgecolor=ecol,facecolor=bcol))
plt.text(1.4,3,'Core \n thermally \n stratified',bbox=dict(edgecolor=ecol,facecolor=bcol))
if save == True:
    plt.savefig(f'../Plots/{folder}timing_heatmap.png',dpi=450,bbox_inches='tight') 
    
#%% Massive stacked barchart
nruns = len(var_results)
yplot = np.arange(nruns)*2
ytick_val = np.array([0,5,6,16,17,21,22,24,25,27,28,32,33,37,38,42])*2
maxB_norm = (var_results['max_B']-min(var_results['max_B']))/(max(var_results['max_B'])-min(var_results['max_B']))
plt.figure(figsize=[10,10])
for i in range(nruns):
    plt.barh(yplot[i],var_results.loc[i,'magoff_1']-var_results.loc[i,'magon_1'],left=var_results.loc[i,'magon_1'],color='navy',alpha=maxB_norm[i])
    plt.barh(yplot[i],var_results.loc[i,'magoff_2']-var_results.loc[i,'magon_2'],left=var_results.loc[i,'magon_2'],color='cornflowerblue',alpha=maxB_norm[i])
plt.xscale('log')  
plt.xlabel('Time/Myr') 
plt.ylabel('Variable') 
plt.yticks(ytick_val,ytick_lab)
if save == True:
    plt.savefig(f'../Plots/{folder}timing_bars.png',dpi=450,bbox_inches='tight') 

#%% Same plot without 60Fe
var_results2 = var_results[(var_results['run']<31)|(var_results['run']>35)]
var_results2.reset_index(inplace=True,drop=True)
nruns2 = len(var_results2)
yplot2 = np.arange(nruns2)*2
ytick_val2 = np.array([0,5,6,16,17,21,22,24,25,27,28,32,33,37])*2
ytick_lab2 = ytick_lab.copy()
#ytick_lab2 = ytick_lab2.remove('$^{{60}}Fe/^{{56}}Fe$=0.0e+00')
#ytick_lab2 = ytick_lab2.remove('$^{{60}}Fe/^{{56}}Fe$=6.0e-07',inplace=True)
#the code above here isn't working
ytick_lab2=['$\\phi_{{RCMF}}$=0.2',
 '$\\phi_{{RCMF}}$=0.5',
 '$\\eta_0$=1.0e+14Pas',
 '$\\eta_0$=1.0e+24Pas',
 '$\\beta$=0.01$K^{-1}$',
 '$\\beta$=0.04$K^{-1}$',
 '$\\eta_l$ =1.0e+00Pas',
 '$\\eta_l$ =1.0e+02Pas',
 '$X_{{s,0}}$=3e+01wt %',
 '$X_{{s,0}}$=3e+01wt %',
 '$\\alpha_n$=2e+01',
 '$\\alpha_n$=4e+01',
 'radius=1e+05km',
 'radius=5e+05km']
maxB_norm2 = (var_results2['max_B']-min(var_results2['max_B']))/(max(var_results2['max_B'])-min(var_results2['max_B']))
plt.figure(figsize=[10,10])
for i in range(nruns2):
    plt.barh(yplot[i],var_results2.loc[i,'magoff_1']-var_results2.loc[i,'magon_1'],left=var_results2.loc[i,'magon_1'],color='navy',alpha=maxB_norm2[i])
    plt.barh(yplot[i],var_results2.loc[i,'magoff_2']-var_results2.loc[i,'magon_2'],left=var_results2.loc[i,'magon_2'],color='cornflowerblue',alpha=maxB_norm2[i])
#plt.xscale('log')  
plt.xlabel('Time/Myr') 
plt.ylabel('Variable') 
plt.yticks(ytick_val2,ytick_lab2)
if save == True:
    plt.savefig(f'../Plots/{folder}timing_bars_nofe.png',dpi=450,bbox_inches='tight') 