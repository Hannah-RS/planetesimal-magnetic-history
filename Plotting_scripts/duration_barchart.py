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
from plot_params import subfolders, labels, units, logs, variables, Myr
folder = 'Paper_run4/'
varlabels = []
ytick_lab = []
ytick_lab2 = []
ytick_val2 = []
paths = []
runval = np.array([])
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
    logvar = logs[var]
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
    paths.extend([path]*nrun) #save paths for opening files later
    runval = np.concatenate([runval,var_data['run'].to_numpy()])
    #for labelling min and max values
    if logvar == True:
        ytick_lab.append(f'{varlab}={var1:.1e}{unit}')
        ytick_lab.extend(['']*(nrun-2))
        ytick_lab.append(f'{varlab}={var2:.1e}{unit}')
    else:
        ytick_lab.append(f'{varlab}={var1:.1g}{unit}')
        ytick_lab.extend(['']*(nrun-2))
        ytick_lab.append(f'{varlab}={var2:.1g}{unit}')   
    #for labelling midpoints
    ytick_lab2.append(varlab)   
    ytick_val2.append(len(ytick_lab)*2-nrun) #len(ytick_lab)-nrun/2 but yplot gets doubled later, so x2
    
#%% Find max B values
nruns = len(var_results)
B1max = np.zeros([nruns])
B2max = np.zeros([nruns])
for i in range(nruns):
    #load B data
    run = int(runval[i])
    #import data
    npzfile = np.load(f'{paths[i]}run_{run}.npz')
    t = npzfile['t']/Myr
    B = npzfile['B']
    #find max value in first interval
    toff1 = var_results.loc[var_results['run']==run,'magoff_1'].values[0]
    if toff1 > 0:
        B1max[i] = np.max(B[t<toff1])
        #find max value in second interval
        B2max[i] = np.max(B[t>toff1])
    
#normalise
maxB1_norm = (B1max-min(B1max))/(max(B1max)-min(B1max))
maxB2_norm = (B2max-min(B2max))/(max(B2max)-min(B2max))
#%% Massive stacked barchart
yplot = np.arange(nruns)*2

fig, ax = plt.subplots(ncols=1)
for i in range(nruns):
    run = int(runval[i])
    #find max B for first and second period
    ax.barh(yplot[i],var_results.loc[var_results['run']==run,'magoff_1']-var_results.loc[var_results['run']==run,'magon_1'],left=var_results.loc[var_results['run']==run,'magon_1'],color=barcol,alpha=maxB1_norm[i])
    ax.barh(yplot[i],var_results.loc[var_results['run']==run,'magoff_2']-var_results.loc[var_results['run']==run,'magon_2'],left=var_results.loc[var_results['run']==run,'magon_2'],color=barcol2,alpha=maxB2_norm[i])
    ax.vlines(var_results.loc[var_results['run']==run,'tsolid_start'],yplot[i]-0.4,yplot[i]+0.4,color='black')
    if (ytick_lab[i]!='')&(ytick_lab[i-1]!='')&(i>0):
        ax.hlines(yplot[i]-1,0.8,500,color='grey',linestyle='dashed',alpha=0.5,linewidth=0.5) 
        
#plt.xscale('log') 

ax.set_xlabel('Time/Myr') 
ax.set_ylabel('Variable') 
#ax.set_yticks(yplot,ytick_lab)
ax.set_yticks(ytick_val2,ytick_lab2)
#first period colourbar
transparency_ticks = 50
color_out =[]
for trans in range(transparency_ticks):
    color_out.append(np.ndarray.tolist(mpl.colors.hsv_to_rgb((0.6667, trans/transparency_ticks, 0.77))))
cmap = mpl.colors.ListedColormap(color_out)
norm = mpl.colors.Normalize(vmin=min(B1max)/1e-6,vmax=max(B1max)/1e-6)
cax = fig.add_axes([0.73, 0.17, 0.01, 0.3])
mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='vertical',label='Max B/$\\mu T$')

#second period colourbar
transparency_ticks = 50
color_out =[]
for trans in range(transparency_ticks):
    color_out.append(np.ndarray.tolist(mpl.colors.hsv_to_rgb((0.525, trans/transparency_ticks, 0.79))))
cmap2 = mpl.colors.ListedColormap(color_out)
norm2 = mpl.colors.Normalize(vmin=min(B2max)/1e-6,vmax=max(B2max)/1e-6)
cax2 = fig.add_axes([0.83, 0.17, 0.01, 0.3])
mpl.colorbar.ColorbarBase(cax2, cmap=cmap2, norm=norm2, orientation='vertical')

if save == True:
    plt.savefig(f'../Plots/{folder}timing_bars.pdf',dpi=450,bbox_inches='tight') 