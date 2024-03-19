#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Scatter plot of peak magnetic field strength for thermal vs compositional dynamos
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

#%% Load data for chosen variable
from plot_params import subfolders, variables, f0, Xs_eutectic

folder = 'Paper_run300km/'
savefolder='Icarus_paper/'
data  = pd.read_csv(f'../Results_combined/{folder}all_sucess_info.csv',skiprows=[1])
if folder == 'Paper_run300km/': #remove r variable runs
    data = data[data['run']<=39]
save = False

#filter out data for different radii
#data = data[data['r']==300e3]
#%% Make colorbar
cmap = plt.colormaps['viridis']
bounds = [-0.5e-9,0.5e-9,0.5e-8,1.5e-8,0.5e-7,1.5e-7,1e-6]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

#%% Data for average field strengths
#add blank columns to dataframe
data['Btherm_av']=""
data['Bcomp_av']=""
f0=0.999
for var in variables[:-1]: #includes r variation but won't be saved
    path = '../Results_combined/'+folder+f"params_{subfolders[var]}/"
    var_data = pd.read_csv(path+'auto_params.csv',skiprows=[1])

    for run in var_data['run']: 
    
        #import data
        npzfile = np.load(f'{path}run_{int(run)}.npz')
        B = npzfile['B']
        Rem = npzfile['Rem']
        f = npzfile['f']
        Xs = npzfile['Xs']
        Btherm = B[(f>=f0)&(Rem>=10)]
        Btherm_av = np.mean(Btherm)
        Bcomp = B[(f<f0)&(Xs<Xs_eutectic)]
        Bcomp_av = np.mean(Bcomp)
        #set to correct data frame column
        data.loc[data['run']==run,'Btherm_av']=Btherm_av
        data.loc[data['run']==run,'Bcomp_av']=Bcomp_av

#%% Make data for straight lines on first plot
x = np.linspace(0,max(data['Btherm_av'])+5e-6,50)/1e-6
#%% Compare average field strengths
#Make sure to adjust limits correctly
fig = plt.figure()
plt.scatter(data['Btherm_av']/1e-6,data['Bcomp_av']/1e-6,c=data['Fe0'],norm=norm,cmap=cmap)
plt.plot(x,0.5*x,linestyle='dashed',color='black',alpha=0.5)
plt.plot(x,x,linestyle='dashed',color='black',alpha=0.5)
plt.plot(x,2*x,linestyle='dashed',color='black',alpha=0.5)
#plt.ylim([5,15])
#plt.xlim([2.5,15])
plt.xlabel('Average thermal dynamo \n field strength /$\mu T$')
plt.ylabel('Average compositional dynamo \n field strength /$\mu T$')
# plt.text(6,13,"2B$_{therm}$",bbox=dict(facecolor='white', 
#                                                       edgecolor='black',alpha=0.7))
# plt.text(13,13,"B$_{therm}$",bbox=dict(facecolor='white', 
#                                                       edgecolor='black',alpha=0.7))
# plt.text(13,6.5,"0.5B$_{therm}$",bbox=dict(facecolor='white', 
#                                                       edgecolor='black',alpha=0.7))
cax = fig.add_axes([0.78, 0.3, 0.01, 0.3])
fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm),cax=cax,ticks=[0,1e-9,1e-8,1e-7,6e-7]
             ,spacing='uniform',orientation='vertical',label='$^{60}Fe/^{56}Fe$')

if save == True:
    plt.savefig(f'../Plots/{savefolder}Bav_comp.png',dpi=450,bbox_inches='tight') 

#%% Make data for straight lines on second plot
x2 = np.linspace(0,max(data['max_Btherm'])+5e-6,50)/1e-6    
#%% Compare peak field strengths
fig = plt.figure()
plt.scatter(data['max_Btherm']/1e-6,data['max_Bcomp']/1e-6,c=data['Fe0'],norm=norm,cmap=cmap)
plt.plot(x2,0.5*x2,linestyle='dashed',color='black',alpha=0.5)
plt.plot(x2,x2,linestyle='dashed',color='black',alpha=0.5)
plt.plot(x2,2*x2,linestyle='dashed',color='black',alpha=0.5)
plt.ylim([8,20])
plt.xlim([5,40])
plt.xlabel('Peak thermal dynamo \n field strength /$\mu T$')
plt.ylabel('Peak compositional dynamo \n field strength /$\mu T$')
plt.text(8,18.5,"2B$_{therm}$",bbox=dict(facecolor='white', 
                                                      edgecolor='black',alpha=0.7))
plt.text(17,18.5,"B$_{therm}$",bbox=dict(facecolor='white', 
                                                      edgecolor='black',alpha=0.7))
plt.text(33,18.5,"0.5B$_{therm}$",bbox=dict(facecolor='white', 
                                                      edgecolor='black',alpha=0.7))
cax = fig.add_axes([0.78, 0.3, 0.01, 0.3])
fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm),cax=cax,ticks=[0,1e-9,1e-8,1e-7,6e-7]
             ,spacing='uniform',orientation='vertical',label='$^{60}Fe/^{56}Fe$')

if save == True:
    plt.savefig(f'../Plots/{savefolder}Bpeak_comp.png',dpi=450,bbox_inches='tight') 
    
#%% Compositional dynamo duration as fraction of core solidification
data['frac_dur'] = (data['magoff_2']-data['tsolid_start'])/(data['tsolid']-data['tsolid_start'])
data2 = data[data['frac_dur']>0] #remove lowest eta0
data3 = data[data['frac_dur']<0]
#one solidification epoch
data.loc[data['frac_dur']<0,'frac_dur'] = (data['magoff_1']-data['tsolid_start'])/(data['tsolid']-data['tsolid_start'])
data3.loc[:,'frac_dur'] = (data3['magoff_1']-data3['tsolid_start'])/(data3['tsolid']-data3['tsolid_start'])
#three solidification epochs
data4 = data3[data3['frac_dur']<0]
data3 = data3[data3['frac_dur']>0]
data.loc[data['frac_dur']<0,'frac_dur'] = (data['magoff_3']-data['tsolid_start'])/(data['tsolid']-data['tsolid_start'])
data4.loc[:,'frac_dur'] = (data4['magoff_3']-data4['tsolid_start'])/(data4['tsolid']-data4['tsolid_start'])
#eutectic data
data4.loc[data4['Xs_0']==Xs_eutectic,'frac_dur']=np.nan
data.loc[data['Xs_0']==Xs_eutectic,'frac_dur']=np.nan

plt.figure()
plt.hist(data3['frac_dur'])
plt.hist(data2['frac_dur'])
plt.hist(data4['frac_dur'])
plt.xlabel('Fraction of core solidification for which a dynamo is generated')
plt.ylabel('Number of runs')

plt.figure()
plt.hist(data['frac_dur'])
plt.xlabel('Fraction of core solidification for which a dynamo is generated')
plt.ylabel('Number of runs')
#plt.savefig('../Plots/fraction_duration.png')