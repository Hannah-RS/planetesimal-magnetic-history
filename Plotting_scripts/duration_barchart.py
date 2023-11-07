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
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap, ListedColormap


#%% Load data for chosen variable
folder = 'Paper_run4/'
subfolders = {'rcmf':1,'eta0':2,'beta':3,'etal':4,'Xs_0':5,'Fe0':6,'alpha_n':7,'r':8}
labels = {'rcmf':'$\\phi_{{RCMF}}$','eta0':'$\\eta_0$','beta':'$\\beta$','etal':'$\\eta_l$ ','Xs_0':'$X_{{s,0}}$','Fe0':'$^{{60}}Fe/^{{56}}Fe$','alpha_n':'$\\alpha_n$','r':'radius'}
units = {'rcmf':'','eta0':'Pas','beta':'$K^{-1}$','etal':'Pas','Xs_0':'wt %','Fe0':'','alpha_n':'','r':'km'}
logs =[False,True,False,True,False,True,False,False]
variables = ['rcmf','eta0','beta','etal','Xs_0','Fe0','alpha_n','r']

varlabels = []
ytick_lab = []
Myr = 365*24*3600*1e6 #number of s in Myr

#plot unchanged B
bcol = 'white'
ecol = 'gray'
barcol = '#0202c4'
barcol2= '#04acc9'
save = True

#%% Start plot
plt.figure(figsize=[12,10])
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
        ytick_lab.extend(['']*(nrun-2))
        ytick_lab.append(f'{varlab}={var2:.1e}{unit}')
    else:
        ytick_lab.append(f'{varlab}={var1:.1g}{unit}')
        ytick_lab.extend(['']*(nrun-2))
        ytick_lab.append(f'{varlab}={var2:.1g}{unit}')   

#%% Massive stacked barchart
nruns = len(var_results)
yplot = np.arange(nruns)*2
#ytick_val = np.array([0,5,6,16,17,21,22,24,25,27,28,32,33,37,38,42])*2
maxB_norm = (var_results['max_B']-min(var_results['max_B']))/(max(var_results['max_B'])-min(var_results['max_B']))
fig, ax = plt.subplots(ncols=1)
for i in range(nruns):
    ax.barh(yplot[i],var_results.loc[i,'magoff_1']-var_results.loc[i,'magon_1'],left=var_results.loc[i,'magon_1'],color=barcol,alpha=maxB_norm[i])
    ax.barh(yplot[i],var_results.loc[i,'magoff_2']-var_results.loc[i,'magon_2'],left=var_results.loc[i,'magon_2'],color=barcol2,alpha=maxB_norm[i])
    ax.vlines(var_results.loc[i,'tsolid_start'],yplot[i]-0.4,yplot[i]+0.4,color='black')
    if (ytick_lab[i]!='')&(ytick_lab[i-1]!='')&(i>0):
        ax.hlines(yplot[i]-1,0.8,500,color='grey',linestyle='dashed',alpha=0.5,linewidth=0.5) 
#plt.xscale('log') 

ax.set_xlabel('Time/Myr') 
ax.set_ylabel('Variable') 
ax.set_yticks(yplot,ytick_lab)

transparency_ticks = 50
color_out =[]
for trans in range(transparency_ticks):
    color_out.append(np.ndarray.tolist(mpl.colors.hsv_to_rgb((0.6667, trans/transparency_ticks, 0.77))))
cmap = mpl.colors.ListedColormap(color_out)
norm = mpl.colors.Normalize(vmin=min(var_results['max_B'])/1e-6, vmax=max(var_results['max_B'])/1e-6)
cax = fig.add_axes([0.73, 0.17, 0.01, 0.3])

mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='vertical',label='Max B/$\\mu T$')
transparency_ticks = 50
color_out =[]
for trans in range(transparency_ticks):
    color_out.append(np.ndarray.tolist(mpl.colors.hsv_to_rgb((0.525, trans/transparency_ticks, 0.79))))
cmap = mpl.colors.ListedColormap(color_out)
norm = mpl.colors.Normalize(vmin=min(var_results['max_B'])/1e-6, vmax=max(var_results['max_B'])/1e-6)
cax = fig.add_axes([0.83, 0.17, 0.01, 0.3])

mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='vertical')

if save == True:
    plt.savefig(f'../Plots/{folder}timing_bars.png',dpi=450,bbox_inches='tight') 

#%% Same plot without 60Fe
# var_results2 = var_results[(var_results['run']<31)|(var_results['run']>35)]
# var_results2.reset_index(inplace=True,drop=True)
# nruns2 = len(var_results2)
# yplot2 = np.arange(nruns2)*2
# ytick_val2 = np.array([0,5,6,16,17,21,22,24,25,27,28,32,33,37])*2
# ytick_lab2 = ytick_lab.copy()
# #ytick_lab2 = ytick_lab2.remove('$^{{60}}Fe/^{{56}}Fe$=0.0e+00')
# #ytick_lab2 = ytick_lab2.remove('$^{{60}}Fe/^{{56}}Fe$=6.0e-07',inplace=True)
# #the code above here isn't working
# ytick_lab2=['$\\phi_{{RCMF}}$=0.2',
#  '$\\phi_{{RCMF}}$=0.5',
#  '$\\eta_0$=1.0e+14Pas',
#  '$\\eta_0$=1.0e+24Pas',
#  '$\\beta$=0.01$K^{-1}$',
#  '$\\beta$=0.04$K^{-1}$',
#  '$\\eta_l$ =1.0e+00Pas',
#  '$\\eta_l$ =1.0e+02Pas',
#  '$X_{{s,0}}$=3e+01wt %',
#  '$X_{{s,0}}$=3e+01wt %',
#  '$\\alpha_n$=2e+01',
#  '$\\alpha_n$=4e+01',
#  'radius=1e+05km',
#  'radius=5e+05km']
# maxB_norm2 = (var_results2['max_B']-min(var_results2['max_B']))/(max(var_results2['max_B'])-min(var_results2['max_B']))
# plt.figure(figsize=[10,10])
# for i in range(nruns2):
#     plt.barh(yplot[i],var_results2.loc[i,'magoff_1']-var_results2.loc[i,'magon_1'],left=var_results2.loc[i,'magon_1'],color='navy',alpha=maxB_norm2[i])
#     plt.barh(yplot[i],var_results2.loc[i,'magoff_2']-var_results2.loc[i,'magon_2'],left=var_results2.loc[i,'magon_2'],color='cornflowerblue',alpha=maxB_norm2[i])
# #plt.xscale('log')  
# plt.xlabel('Time/Myr') 
# plt.ylabel('Variable') 
# plt.yticks(ytick_val2,ytick_lab2)
# if save == True:
#     plt.savefig(f'../Plots/{folder}timing_bars_nofe.png',dpi=450,bbox_inches='tight') 