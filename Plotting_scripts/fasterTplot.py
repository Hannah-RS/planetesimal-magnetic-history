#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot things notebook is doing very slowly - need to tidy this up later and add more commenting
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

#find unstable indices for plotting core stratification
folder = 'Singlevar1/'
subfolders = {'rcmf':1,'eta0':2,'beta':3,'etal':4,'Xs_0':5,'Fe0':6,'alpha_n':7,'r':8}
labels = {'rcmf':'$\\phi_{{RCMF}}$','eta0':'$\\eta_0$','beta':'$\\beta$','etal':'$\\eta_l$ ','Xs_0':'$X_{{s,0}}$','Fe0':'$^{{60}}Fe/^{{56}}Fe$','alpha_n':'$\\alpha_n$','r':'radius'}
Myr = 365*24*3600*1e6 #number of s in Myr

var = 'r'
unit = 'km' #unit of variable
varlab = labels[var]
logvar = False
save = True
path = '../Results_combined/'+folder+f"params_{subfolders[var]}/"

#find run numbers
var_data = pd.read_csv(path+'auto_params.csv',skiprows=[1])
var_results = pd.read_csv(path+'run_results.csv',skiprows=[1])
data = pd.merge(var_data,var_results,on='run')
minrun = min(var_data['run'])
maxrun = max(var_data['run'])
nrun = len(var_data)
var2=var_data.loc[var_data['run']==maxrun,var].values[0]
#End member 2
run = int(maxrun)
npzfile = np.load(f'{path}run_{run}_diff.npz')
Tdiff = npzfile['Tdiff'] 
tdiff = npzfile['t_diff']/Myr
d0_diff = npzfile['d0']

npzfile = np.load(f'{path}run_{run}.npz')
T_profile = npzfile['T_profile']
t2 = npzfile['t']/Myr #time in Myr 
d02 = npzfile['d0'] 
min_unstable = npzfile['min_unstable'] 
Flux = npzfile['Flux']
Fs2 = Flux[0]
Fcmb2 = Flux[1]
Fad2 = Flux[2]
Frad2 = Flux[3]

#time for field to be on
on21=var_results.loc[var_results['run']==maxrun,'magon_1'].values[0]
off21=var_results.loc[var_results['run']==maxrun,'magoff_1'].values[0]
on22=var_results.loc[var_results['run']==maxrun,'magon_2'].values[0]
off22=var_results.loc[var_results['run']==maxrun,'magoff_2'].values[0]
#get time for switch to conduction
fcond_t2 = var_results.loc[var_results['run']==run,'fcond_t'].values[0]

#Concatenate
Tall_2 = np.hstack((Tdiff,np.transpose(T_profile)))
tall_2 = np.append(tdiff,t2)
d0_all_2 = np.append(d0_diff,d02)
dr = var_data.loc[var_data['run']==run,'dr'].values[0]
r2 = var_data.loc[var_data['run']==run,'r'].values[0]
rplot_2 = np.arange(0,r2+dr,dr)/1e3
r_unstable2=np.array([]) 
for ind in min_unstable:
    r_unstable2 = np.append(r_unstable2,rplot_2[int(ind)])
rc2 = r2/2
with sns.plotting_context('talk'):
    plt.figure(figsize=[20,7.5])
    plt.pcolormesh(tall_2,rplot_2,Tall_2,shading = 'nearest',vmin=200,vmax=1650,rasterized=True)
    plt.hlines(rc2/1e3,min(t2),max(tall_2),linestyle='--',color='black',label='CMB')
    plt.vlines(t2[0],0,r2/1e3,linestyle='-.',label='Differentiation')
    plt.fill_betweenx([0,rc2/5e3],on21,off21,alpha=0,hatch='/',label='dynamo on')
    if on22>0:
        plt.fill_betweenx([0,rc2/5e3],on22,off22,alpha=0,hatch='/')
    plt.plot(t2,r_unstable2,linestyle='dotted',label='Convecting core')
    if np.any(t2/Myr<fcond_t2):
        plt.plot(t2[t2<=fcond_t2],(r2-d02[t2<=fcond_t2])/1e3,linestyle='dashed',label='base of $\delta_0$',color='blue')
        plt.vlines(t2[t2<=fcond_t2][-1],r2/1e3,rc2/1e3,linestyle='dotted',label='conductive mantle',color='red')
    else: #i.e. never becomes conductive
        plt.plot(t2,(r2-d02)/1e3,linestyle='dashed',label='base of $\delta_0$',color='blue')
    #labels and limits
    plt.ylabel('r /km')
    plt.xlabel('t / Myr')
    plt.xlim([0.8,800])
    plt.colorbar(label='T/K')
    plt.title(f"{varlab}={var2} {unit}")
    plt.xscale('log')
    plt.legend(bbox_to_anchor=[0.65,0.1])
    if save == True:
        plt.savefig(f'../Plots/{folder}/Ttalk_{var}.pdf',bbox_inches='tight',format='pdf',dpi=300)