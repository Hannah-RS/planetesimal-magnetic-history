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
folder = 'Paper_run300km/'
savefolder = 'EPSL_paper/'
from plot_params import subfolders, labels, units, logs, variables, Myr
 

varlabels = []
paths = []
ytick_lab = []
ytick_lab2 = []
ytick_val2 = []
runval = np.array([])


save = False
xmax = 12 #how much to zoom by
#%% Load a given variable

for i, var in enumerate(variables[:-1]):

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

#make colormap for B data
import copy
import matplotlib.colors as colors
cmap = copy.copy(plt.get_cmap('viridis'))
cmap.set_under('white')
Bmin = 3 #min B value in muT based on lower limit from data
#Bminfind=[] #use this line to find Bmin
#Bmaxfind = np.max([var_results['max_Btherm'].max(),var_results['max_Bcomp'].max()])/1e-6
Bmax = 40 #based on Bmaxfind
norm = colors.LogNorm(vmin=Bmin,vmax=Bmax)

#make colormap that is all white to indicate gaps in data
cmap2 = (mpl.colors.ListedColormap([ '#004488', '#DDAA33']).with_extremes(over='#BB5566', under='#FFFFFF'))
bounds2 = [0.5,1]
norm2 = mpl.colors.BoundaryNorm(bounds2, cmap2.N)

fig, axes = plt.subplots(1,2,sharey='row',figsize=[7.5,5],gridspec_kw={'width_ratios': [1,3]},tight_layout=True)
for ax in axes: 
    for i in range(nruns):
        #load Rem data
        run = int(runval[i])
        #import data
        npzfile = np.load(f'{paths[i]}run_{run}.npz')
        Rem = npzfile['Rem']
        B = npzfile['B']/1e-6
        t = npzfile['t']/Myr
        Rem_c = 10 #critical Rem
        #average B post onset of solidification to remove white bars
        tsolid_start = var_results.loc[var_results['run']==run,'tsolid_start'].values[0]
        toff2 = var_results.loc[var_results['run']==run,'magoff_2'].values[0]
        toff3 = var_results.loc[var_results['run']==run,'magoff_3'].values[0]
        if toff3 != 0: #3 dynamo periods
            B[(t>tsolid_start)&(t<toff3)]=np.average(B[(t>tsolid_start)&(t<toff3)])
        elif toff2 !=0: #2 dynamo periods
            B[(t>tsolid_start)&(t<toff2)]=np.average(B[(t>tsolid_start)&(t<toff2)])
        else: #1 dynamo period
            toff1 = var_results.loc[var_results['run']==run,'magoff_1'].values[0]
            B[(t>tsolid_start)&(t<toff1)]=np.average(B[(t>tsolid_start)&(t<toff1)])
        #set all B with subccritical Rem before solidification to 0
        B[(t<tsolid_start)&(Rem<Rem_c)]=0
        #Bminfind.append(np.min(B[B>0]))
        
        #make the plot
        ax.pcolormesh(t,yplot[2*i:2*i+2],[B,B],cmap=cmap,norm=norm)
        
        if (ytick_lab[2*i]!='')&(ytick_lab[2*i-1]!='')&(i>0):
            ax.hlines(yplot[2*i]-1,0.8,210,color='black',linestyle='dashed',linewidth=0.5)
        ax.pcolormesh(t,yplot[2*i+1:2*i+3],np.zeros([2,len(Rem)]),cmap=cmap2,norm=norm2) #add gaps between variables
        ax.vlines(var_results.loc[var_results['run']==run,'tsolid_start'],yplot[2*i]-0.5,yplot[2*i]+0.5,color='black')
         
    
    ax.set_ylim([yplot[0]-1,yplot[-1]+0.5])
axes[1].set_xlabel('Time after CAI formation /Ma') 
axes[0].set_xlim([0.8,xmax])
axes[1].set_xlim([xmax,300])
axes[0].set_yticks(ytick_val2,ytick_lab2)
axes[1].set_title(f'>{xmax} Ma')
axes[0].set_title(f'First {xmax} Ma')
cax = fig.add_axes([0.87, 0.17, 0.01, 0.3])
colorbar = fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm),cax=cax,
             orientation='vertical', label='magnetic field strength /$\\mu T$')
colorbar.set_ticks([3,10,20,30,40],labels=[3,10,20,30,40])


if save == True:
    plt.savefig(f'../Plots/{savefolder}Rem_bars_full.png',dpi=450,bbox_inches='tight') 

#%% Minimum value for colormap
#print(np.min(Bminfind)) #muT