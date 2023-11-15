#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rem as a function of time over all variables
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap, ListedColormap


#%% Load data for chosen variable
folder = 'Paper_run4/'
from plot_params import subfolders, labels, units, logs, variables, Myr
 

varlabels = []
paths = []
ytick_lab = []
ytick_lab2 = []
ytick_val2 = []
runval = np.array([])


save = True
xmax = 12 #how much to zoom by
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
    if logvar == True:
        ytick_lab.append(f'{varlab}={var1:.1e}{unit}')
        ytick_lab.extend(['']*(2*nrun-2))
        ytick_lab.append(f'{varlab}={var2:.1e}{unit}')
    else:
        ytick_lab.append(f'{varlab}={var1:.1g}{unit}')
        ytick_lab.extend(['']*(2*nrun-2))
        ytick_lab.append(f'{varlab}={var2:.1g}{unit}')   
    
    #for labelling midpoints
    ytick_lab2.append(varlab)   
    ytick_val2.append(len(ytick_lab)-nrun) #len(ytick_lab)-nrun/2 but yplot gets doubled later, so x2
        
#%% Massive stacked barchart
nruns = len(runval)
yplot = np.arange(2*nruns+1)
#convert yticks to array
ytick_arr = np.array([ytick_lab])

#make colormap for Rem data
#cmap = (mpl.colors.ListedColormap([ '#7D9CEB', '#DDAA33']).with_extremes(over='#BB5566', under='#FFFFFF'))
cmap = (mpl.colors.ListedColormap([ '#7D9CEB', '#A77B00']).with_extremes(over='#87263C', under='#EBEBEB'))

bounds = [10,40,100]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

#make colormap that is all white to indicate gaps in data
#cmap2 = (mpl.colors.ListedColormap([ '#004488', '#DDAA33']).with_extremes(over='#BB5566', under='#FFFFFF'))
cmap2 = (mpl.colors.ListedColormap([ '#004488', '#DDAA33']).with_extremes(over='#BB5566', under='#EBEBEB'))
bounds2 = [0.5,1]
norm2 = mpl.colors.BoundaryNorm(bounds2, cmap2.N)

fig, axes = plt.subplots(1,2,sharey='row',figsize=[7.5,5],gridspec_kw={'width_ratios': [3,1]},tight_layout=True)
for ax in axes: 
    ax.fill_betweenx([-1,yplot[-1]],0,500,color='#EBEBEB')
    for i in range(nruns):
        #load Rem data
        run = int(runval[i])
        #import data
        npzfile = np.load(f'{paths[i]}run_{run}.npz')
        Rem = npzfile['Rem']
        t = npzfile['t']/Myr
        
        #average Rem post onset of solidification to remove white bars
        tsolid_start = var_results.loc[var_results['run']==run,'tsolid_start'].values[0]
        toff2 = var_results.loc[var_results['run']==run,'magoff_2'].values[0]
        toff3 = var_results.loc[var_results['run']==run,'magoff_3'].values[0]
        if toff3 != 0: #3 dynamo periods
            Rem[(t>tsolid_start)&(t<toff3)&(Rem==0)]=np.average(Rem[(t>tsolid_start)&(t<toff3)])
        elif toff2 !=0: #2 dynamo periods
            Rem[(t>tsolid_start)&(t<toff2)&(Rem==0)]=np.average(Rem[(t>tsolid_start)&(t<toff2)])
        else: #1 dynamo period
            toff1 = var_results.loc[var_results['run']==run,'magoff_1'].values[0]
            Rem[(t>tsolid_start)&(t<toff1)&(Rem==0)]=np.average(Rem[(t>tsolid_start)&(t<toff1)])
        
        #make the plot
        ax.pcolormesh(t,yplot[2*i:2*i+2],[Rem,Rem],cmap=cmap,norm=norm)
        
        if (ytick_lab[2*i]!='')&(ytick_lab[2*i-1]!='')&(i>0):
            ax.hlines(yplot[2*i]-1,0.8,210,color='black',linestyle='dashed',linewidth=0.5)
        ax.pcolormesh(t,yplot[2*i+1:2*i+3],np.zeros([2,len(Rem)]),cmap=cmap2,norm=norm2) #add gaps between variables
        ax.vlines(var_results.loc[var_results['run']==run,'tsolid_start'],yplot[2*i]-0.5,yplot[2*i]+0.5,color='black')
         
    ax.set_xlabel('Time/Myr') 
    ax.set_ylim([yplot[0]-1,yplot[-1]+0.5])
axes[0].set_xlim([0.8,500])
axes[1].set_xlim([0.8,xmax])
axes[0].set_yticks(ytick_val2,ytick_lab2)
axes[0].set_title('Full history')
axes[1].set_title(f'First {xmax} Myr')
cax = fig.add_axes([0.6, 0.17, 0.01, 0.3])
fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm),cax=cax,extend='both',ticks=bounds,spacing='proportional',
             extendfrac=0.5, orientation='vertical', label='Rem')


if save == True:
    plt.savefig(f'../Plots/{folder}Rem_bars_ugly.png',dpi=450,bbox_inches='tight') 

#%% Continuous B version - this sucks so have commented it out
#if do want to use it need to make 0 values unfilled rather than lightest colour 
# need perceptially uniform colourmap
# maxB = max(var_results['max_B'])/1e-6
# minB = 0 #lowest B

# plt.figure()
# minrun = 1
# for i in range(nruns):
#     #load Rem data
#     run = int(runval[i])
#     #import data
#     npzfile = np.load(f'{paths[i]}run_{run}.npz')
#     B = npzfile['B']/1e-6 #B in muT
#     t = npzfile['t']/Myr
    
#     plt.pcolormesh(t,yplot[2*i:2*i+2],[B,B],vmin=minB,vmax=maxB,cmap='YlOrRd')
#     if i==(nruns-1):
#         plt.colorbar(label='B/ $\mu$T')
        
#     plt.pcolormesh(t,yplot[2*i+1:2*i+3],np.zeros([2,len(B)]),cmap=cmap2,norm=norm2) #add gaps between variables
#     plt.vlines(var_results.loc[var_results['run']==run,'tsolid_start'],yplot[2*i]-0.5,yplot[2*i]+0.5,color='black')
#     if (ytick_lab[2*i]!='')&(ytick_lab[2*i-1]!='')&(i>0):
#         plt.hlines(yplot[2*i]-1,0.8,500,color='grey',linestyle='dashed',alpha=0.5,linewidth=0.5) 
    

# plt.xlabel('Time/Myr') 
# #plt.xlim([0.8,20])
# plt.yticks(yplot[np.where(ytick_arr!='')[1]],ytick_arr[0,np.where(ytick_arr!='')[1]])

