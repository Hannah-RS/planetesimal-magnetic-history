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
    #get tick labels correct format
    labs = data[f'{var}']
    ytick_lab_new = []
    for lab in labs:
        if var == 'beta':
            ytick_lab_new.append(f'{lab:.5f} $K^{{-1}}$')
        elif var == 'eta0':
            ytick_lab_new.append(f'{lab:.0e} Pas')
        elif var=='Fe0':
            ytick_lab_new.append(f'{lab:.0e}')
        elif (var == 'alpha_n'):
            ytick_lab_new.append(f'{lab:.0f}')
        elif var == 'etal':
            ytick_lab_new.append(f'{lab:.0f} Pas')
        elif var == 'r':
            ytick_lab_new.append(f'{lab/1000:.0f} km')
        elif (var == 'Xs_0'):
            ytick_lab_new.append(f'{lab:.2f} wt %')
        else:
            ytick_lab_new.append(f'{lab:.2f}')
    ytick_lab = np.concatenate([ytick_lab,ytick_lab_new,['']]) #add a gap 
    
#%% Massive stacked barchart
nruns = len(runval)
yplot = np.arange(2*(nruns+7))

#make colormap for Rem data
cmap = (mpl.colors.ListedColormap([ '#7D9CEB', '#DDAA33']).with_extremes(over='#BB5566', under='#FFFFFF'))

bounds = [10,40,100]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

fig, axes = plt.subplots(1,2,sharey='row',figsize=[7.5,6],gridspec_kw={'width_ratios': [1,3]},tight_layout=True)
for ax in axes: 
    j = 0 
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
            Rem[(t>tsolid_start)&(t<toff3)]=np.average(Rem[(t>tsolid_start)&(t<toff3)])
        elif toff2 !=0: #2 dynamo periods
            Rem[(t>tsolid_start)&(t<toff2)]=np.average(Rem[(t>tsolid_start)&(t<toff2)])
        else: #1 dynamo period
            toff1 = var_results.loc[var_results['run']==run,'magoff_1'].values[0]
            Rem[(t>tsolid_start)&(t<toff1)]=np.average(Rem[(t>tsolid_start)&(t<toff1)])
        
        if (ytick_lab[j]==''): #don't plot in gaps
            ax.hlines(yplot[2*j],0.8,210,color='grey',linestyle='dashed',alpha=0.5,linewidth=0.5)
            j=j+1
            
        #make the plot
        ax.pcolormesh(t,yplot[2*j:2*j+2],[Rem,Rem],cmap=cmap,norm=norm)
        ax.barh(yplot[2*j+1],t[-1],color='white',height=1)
    
        ax.vlines(var_results.loc[var_results['run']==run,'tsolid_start'],yplot[2*j]-0.7,yplot[2*j]+0.7,color='black')
        j = j+1
        
    
    ax.set_ylim([yplot[0]-1,yplot[-1]+0.5]) 
axes[0].set_xlim([0.8,xmax])
axes[1].set_xlim([xmax,600])
axes[0].set_yticks(yplot[::2],ytick_lab[:-1])
axes[0].tick_params(axis='y', labelsize=7) #make y labels smaller
axes[1].set_title(f'>{xmax} Ma',fontsize=10)
axes[0].set_title(f'First {xmax} Ma',fontsize=10)
cax = fig.add_axes([0.87, 0.17, 0.01, 0.3])
fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm),cax=cax,extend='both',ticks=bounds,spacing='proportional',
             extendfrac=0.5, orientation='vertical', label='Re$_m$')
fig.suptitle('Time after CAI formation /Ma',y=0,fontsize=10)

if save == True:
    plt.savefig(f'../Plots/{savefolder}Rem_bars_full.png',dpi=500,bbox_inches='tight') 

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

