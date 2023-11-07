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
subfolders = {'rcmf':1,'eta0':2,'beta':3,'etal':4,'Xs_0':5,'Fe0':6,'alpha_n':7,'r':8}
labels = {'rcmf':'$\\phi_{{RCMF}}$','eta0':'$\\eta_0$','beta':'$\\beta$','etal':'$\\eta_l$ ','Xs_0':'$X_{{s,0}}$','Fe0':'$^{{60}}Fe/^{{56}}Fe$','alpha_n':'$\\alpha_n$','r':'radius'}
units = {'rcmf':'','eta0':'Pas','beta':'$K^{-1}$','etal':'Pas','Xs_0':'wt %','Fe0':'','alpha_n':'','r':'km'}
logs =[False,True,False,True,False,True,False,False]
variables = ['rcmf','eta0','beta','etal','Xs_0','Fe0','alpha_n','r']

varlabels = []
paths = []
ytick_lab = []
Myr = 365*24*3600*1e6 #number of s in Myr

save = False
zoom = True #do you want to zoom into early times
xmax = 10 #how much to zoom by
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
    paths.extend([path]*nrun) #save paths for opening files later
    
    if logvar == True:
        ytick_lab.append(f'{varlab}={var1:.1e}{unit}')
        ytick_lab.extend(['']*(2*nrun-2))
        ytick_lab.append(f'{varlab}={var2:.1e}{unit}')
    else:
        ytick_lab.append(f'{varlab}={var1:.1g}{unit}')
        ytick_lab.extend(['']*(2*nrun-2))
        ytick_lab.append(f'{varlab}={var2:.1g}{unit}')   

#%% Massive stacked barchart
nruns = len(var_results)
yplot = np.arange(2*nruns+1)
#convert yticks to array
ytick_arr = np.array([ytick_lab])

#make colormap for Rem data
cmap = (mpl.colors.ListedColormap([ '#004488', '#DDAA33']).with_extremes(over='#BB5566', under='#FFFFFF'))
bounds = [10,40,100]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

#make colormap that is all white to indicate gaps in data
cmap2 = (mpl.colors.ListedColormap([ '#004488', '#DDAA33']).with_extremes(over='#BB5566', under='#FFFFFF'))
bounds2 = [0.5,1]
norm2 = mpl.colors.BoundaryNorm(bounds2, cmap2.N)

fig, ax = plt.subplots(ncols=1)
minrun = 1
for i in range(nruns):
    #load Rem data
    run = int(minrun+i)
    #import data
    npzfile = np.load(f'{paths[i]}run_{run}.npz')
    #B = npzfile['B']
    Rem = npzfile['Rem']
    t = npzfile['t']/Myr
    #make the plot
    ax.pcolormesh(t,yplot[2*i:2*i+2],[Rem,Rem],cmap=cmap,norm=norm)
    if zoom == True:
        if (ytick_lab[2*i]!='')&(ytick_lab[2*i-1]!='')&(i>0):
            ax.hlines(yplot[2*i]-1,0.8,500,color='white',linewidth=1)
            
    else:
        if (ytick_lab[2*i]!='')&(ytick_lab[2*i-1]!='')&(i>0):
            ax.hlines(yplot[2*i]-1,0.8,500,color='grey',linestyle='dashed',alpha=0.5,linewidth=0.5)
            ax.pcolormesh(t,yplot[2*i+1:2*i+3],np.zeros([2,len(Rem)]),cmap=cmap2,norm=norm2) #add gaps between variables
    ax.vlines(var_results.loc[i,'tsolid_start'],yplot[2*i]-1,yplot[2*i]+1,color='black')
     


ax.set_xlabel('Time/Myr') 
if zoom ==True:
    ax.set_xlim([0.8,xmax])
ax.set_yticks(yplot[np.where(ytick_arr!='')[1]],ytick_arr[0,np.where(ytick_arr!='')[1]])
cax = fig.add_axes([1, 0.17, 0.01, 0.3])
fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm),cax=cax,extend='both',ticks=bounds,spacing='proportional',
             extendfrac=0.5, orientation='vertical', label='Rem')


if save == True:
    if zoom == True:
        plt.savefig(f'../Plots/{folder}Rem_bars_zoom_stripe.png',dpi=450,bbox_inches='tight') 
    else:
        plt.savefig(f'../Plots/{folder}Rem_bars_full.png',dpi=450,bbox_inches='tight') 

#%% Continuous B version - this sucks so have commented it out
maxB = max(var_results['max_B'])/1e-6
minB = 0 #lowest B

plt.figure()
minrun = 1
for i in range(nruns):
    #load Rem data
    run = int(minrun+i)
    #import data
    npzfile = np.load(f'{paths[i]}run_{run}.npz')
    B = npzfile['B']/1e-6 #B in muT
    t = npzfile['t']/Myr
    
    plt.pcolormesh(t,yplot[2*i:2*i+2],[B,B],vmin=minB,vmax=maxB,cmap='YlOrRd')
    if i==(nruns-1):
        plt.colorbar(label='B/ $\mu$T')
        
    plt.pcolormesh(t,yplot[2*i+1:2*i+3],np.zeros([2,len(B)]),cmap=cmap2,norm=norm2) #add gaps between variables
    plt.vlines(var_results.loc[i,'tsolid_start'],yplot[2*i]-0.5,yplot[2*i]+0.5,color='black')
    if (ytick_lab[2*i]!='')&(ytick_lab[2*i-1]!='')&(i>0):
        plt.hlines(yplot[2*i]-1,0.8,500,color='grey',linestyle='dashed',alpha=0.5,linewidth=0.5) 
    

plt.xlabel('Time/Myr') 
#plt.xlim([0.8,20])
plt.yticks(yplot[np.where(ytick_arr!='')[1]],ytick_arr[0,np.where(ytick_arr!='')[1]])

