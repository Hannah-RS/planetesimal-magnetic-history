#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Heatmap of dynamo on times for Psyche
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

#%% Load data for chosen variable

folder = 'Psyche/'
savefolder='Psyche/'

ytick_lab = []

bcol = 'white'
ecol = 'gray'
barcol = '#0202c4'
barcol2= '#04acc9'
save = True

path = '../Results_combined/'+folder
var_data = pd.read_csv(path+'auto_params.csv',skiprows=[1])
var_results = pd.read_csv(f'../Results_combined/{folder}/run_results.csv',skiprows=[1])

maxrun =int(max(var_data['run']))

#%% Massive stacked barchart
yplot = np.arange(0,maxrun+1,1)
#yplot = np.arange(nruns+7)*2 #add extra blank space

fig, ax = plt.subplots(1,1,figsize=[10,5],tight_layout=True)

#for i in range(1,maxrun+1): #plot in order of run number
    #run = i
for i, run in enumerate(var_data['run']): #plot in order on spreadsheet
    ytick_lab.append(f" $\\frac{{r_c}}{{r}}$={var_data.loc[var_data['run']==run,'rcr'].values[0]},$\\eta_0$ ={var_data.loc[var_data['run']==run,'eta0'].values[0]}Pas,$X_{{S,0}}$ ={var_data.loc[var_data['run']==run,'Xs_0'].values[0]}wt %  ")
    ax.barh(yplot[i],var_results.loc[var_results['run']==run,'magoff_1']-var_results.loc[var_results['run']==run,'magon_1'],left=var_results.loc[var_results['run']==run,'magon_1'],color=barcol)
    ax.barh(yplot[i],var_results.loc[var_results['run']==run,'magoff_2']-var_results.loc[var_results['run']==run,'magon_2'],left=var_results.loc[var_results['run']==run,'magon_2'],color=barcol)
    ax.barh(yplot[i],var_results.loc[var_results['run']==run,'magoff_3']-var_results.loc[var_results['run']==run,'magon_3'],left=var_results.loc[var_results['run']==run,'magon_3'],color=barcol)
    ax.vlines(var_results.loc[var_results['run']==run,'tsolid_start'],yplot[i]-0.4,yplot[i]+0.4,color='cornflowerblue')
   
ax.set_xlabel('Time after CAI formation/Ma')     
ax.set_yticks(yplot[:-1],ytick_lab)
ax.set_title('Dynamo duration')
ax.annotate('Onset of core solidification',(9,1),(20,2),color='cornflowerblue',arrowprops=dict(arrowstyle='simple',facecolor='cornflowerblue',edgecolor='cornflowerblue'))
#ax.set(ylim=[-0.5,5.5])

if save == True:
    plt.savefig(f'../Plots/{savefolder}timing_bars.png',dpi=450,bbox_inches='tight')