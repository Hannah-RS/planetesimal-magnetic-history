#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Heatmap of dynamo on times as a function of variable
Can optionally shade a reference run over the top
This version has a zoom in of early times
This is customised for just plotting fe
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

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
xlim = 12 #limit of zoomed region
save = True
ref = False #do you want reference run over the top
ron = False #do you want r on the plot?

variables = ['Fe0']

#if ron == False:
#    variables = variables[:-1] #exclude r
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
nruns = len(runval)
B1max = np.zeros([nruns])
B2max = np.zeros([nruns])
B3max = np.zeros([nruns])
for i in range(nruns):
    #load B data
    run = int(runval[i])
    #import data
    npzfile = np.load(f'{paths[i]}run_{run}.npz')
    t = npzfile['t']/Myr
    B = npzfile['B']
    #find max value in first interval
    toff1 = var_results.loc[var_results['run']==run,'magoff_1'].values[0]
    toff2 = var_results.loc[var_results['run']==run,'magoff_2'].values[0]
    if toff1 > 0:
        B1max[i] = np.max(B[t<toff1])
    if toff2> 0:
        #find max value in second interval
        B2max[i] = np.max(B[(t>toff1)&(t<toff2)])
        B3max[i] = np.max(B[(t>toff2)])
    #ytick_lab.append(f"{var_results.loc[var_results['run']==run,'Fe0'].values[0]:.0e}")
#normalise, use min value of 0 so everything appears and max value of all B values
maxB = np.max([B1max,B2max,B3max])
maxB1_norm = (B1max-0)/(maxB-0)
maxB2_norm = (B2max-0)/(maxB-0)
maxB3_norm = (B3max-0)/(maxB-0)
#%% Massive stacked barchart
yplot = np.arange(nruns)*2

fig, axes = plt.subplots(1,2,sharey='row',figsize=[7.5,5],gridspec_kw={'width_ratios': [1,3]},tight_layout=True)
for ax in axes: 
    for i in range(nruns):
        run = int(runval[i])
        #find max B for first and second period
        ax.barh(yplot[i],var_results.loc[var_results['run']==run,'magoff_1']-var_results.loc[var_results['run']==run,'magon_1'],left=var_results.loc[var_results['run']==run,'magon_1'],color=barcol,alpha=maxB1_norm[i])
        ax.barh(yplot[i],var_results.loc[var_results['run']==run,'magoff_2']-var_results.loc[var_results['run']==run,'magon_2'],left=var_results.loc[var_results['run']==run,'magon_2'],color=barcol,alpha=maxB2_norm[i])
        ax.barh(yplot[i],var_results.loc[var_results['run']==run,'magoff_3']-var_results.loc[var_results['run']==run,'magon_3'],left=var_results.loc[var_results['run']==run,'magon_3'],color=barcol,alpha=maxB3_norm[i])
        ax.vlines(var_results.loc[var_results['run']==run,'tsolid_start'],yplot[i]-0.4,yplot[i]+0.4,color='black')
        if (ytick_lab[i]!='')&(ytick_lab[i-1]!='')&(i>0):
            ax.hlines(yplot[i]-1,0.8,210,color='grey',linestyle='dashed',alpha=0.5,linewidth=0.5) 
            
        #import run with standard variables for comparison
    if ref == True: #plot reference run
        #ref_results = pd.read_csv(f'../Results_combined/{folder}/Reference_run/run_results.csv',skiprows=[1])
        #ax.fill_betweenx(yplot[:-5],ref_results['magon_1'],ref_results['magoff_1'],alpha=0.1,color='gray')
        #ax.fill_betweenx(yplot[:-5],ref_results['magon_2'],ref_results['magoff_2'],alpha=0.1,color='gray')
        #values from Bryson 2019, single event model, 
        ax.fill_betweenx(yplot[:-5],5,9.7,alpha=0.1,color='gray')
        ax.fill_betweenx(yplot[:-5],100,400,alpha=0.1,color='gray') #estimate Rem>Remc in Fig 1 for lower limit from fig 3 for 300km for upper limit
    
    ax.set_xlabel('Time/Myr') 
    ytick_lab = ['0','$10^{-9}$','$10^{-8}$','$10^{-7}$','$6\\times 10^{-7}$']
    ax.set_yticks(yplot,ytick_lab)
    
    #ax.set_yticks(ytick_val2,ytick_lab2)
#change limit of second graph
axes[0].set_xlim([0,xlim])
axes[0].set_ylabel(f'{varlab}',fontsize=15)
#axes[0].set_ylabel('Variable')
axes[1].set_title('Full history - 300km body')
axes[0].set_title(f'First {xlim} Myr')

#first period colourbar
transparency_ticks = 50
color_out =[]
for trans in range(transparency_ticks):
    color_out.append(np.ndarray.tolist(mpl.colors.hsv_to_rgb((0.6667, trans/transparency_ticks, 0.77))))
cmap = mpl.colors.ListedColormap(color_out)
norm = mpl.colors.Normalize(vmin=0,vmax=maxB/1e-6)
if ron == True:
    cax = fig.add_axes([0.55, 0.17, 0.01, 0.3])
else:
    cax = fig.add_axes([1, 0.3, 0.01, 0.3])
mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='vertical',label='Max B/$\\mu T$')

# #second period colourbar
# transparency_ticks = 50
# color_out =[]
# for trans in range(transparency_ticks):
#     color_out.append(np.ndarray.tolist(mpl.colors.hsv_to_rgb((0.525, trans/transparency_ticks, 0.79))))
# cmap2 = mpl.colors.ListedColormap(color_out)
# norm2 = mpl.colors.Normalize(vmin=0,vmax=maxB/1e-6)
# if ron == True:
#     cax2 = fig.add_axes([0.63, 0.17, 0.01, 0.3])
# else:
#     cax2 = fig.add_axes([1, 0.5, 0.01, 0.3])
# mpl.colorbar.ColorbarBase(cax2, cmap=cmap2, norm=norm2, orientation='vertical')


if save == True:
    if ref == True:
        plt.savefig(f'../Plots/{folder}timing_bars_zoomref.png',dpi=450,bbox_inches='tight') 
    else:
        plt.savefig(f'../Plots/{folder}timing_bars_viscosity.png',dpi=450,bbox_inches='tight')