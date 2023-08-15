#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 15:06:22 2023

@author: exet5460
"""
import numpy as np
import pandas as pd

folder = 'Singlevar1/'
subfolders = {'rcmf':1,'eta0':2,'frht':3,'etal':4,'Xs_0':5,'Fe0':6,'alpha_n':7,'r':8}
labels = {'rcmf':'$\\phi_{{RCMF}}$','eta0':'$\\eta_0$','frht':'frht','etal':'$\\eta_l$','Xs_0':'$X_{{s,0}}$','Fe0':'$^{{60}}Fe/^{{56}}Fe$','alpha_n':'$\\alpha_n$','r':'radius'}
Myr = 365*24*3600*1e6 #number of s in Myr

for j in range(8):

    path = '../Results_combined/'+folder+f"params_{j+1}/"
    
    var_data = pd.read_csv(path+'auto_params.csv',skiprows=[1])
    var_results = pd.read_csv(path+'run_results.csv',skiprows=[1])
    
    minrun = min(var_data['run'])
    maxrun = max(var_data['run'])
    nrun = len(var_data)
    
    peakT = np.zeros([nrun])
    tmax = np.zeros([nrun])
    peak_coreT = np.zeros([nrun])
    tcoremax = np.zeros([nrun])
    
    for i in range(nrun):
        run = int(minrun + i)
        r = var_data.loc[var_data['run']==run,'r'].values[0]
        dr = 500
        nmantle = int((r/dr)/2)
        
        npzfile = np.load(f'{path}run_{run}.npz')
        T_profile = npzfile['T_profile']
        t = npzfile['t']/Myr #time in Myr
    
        peakT[i] = np.amax(T_profile[:,nmantle+1:])
        loc_max = np.where(T_profile[:,nmantle+1:]==peakT[i])[0][0] #take the set of time coordinates and first value (they should all be the same)
        tmax[i] = t[loc_max]
        peak_coreT[i] = np.amax(T_profile[:,:nmantle])
        loc_max = np.where(T_profile[:,:nmantle]==peak_coreT[i])[0][0] #take the set of time coordinates and first value (they should all be the same)
        tcoremax[i] = t[loc_max]
        
    var_results['peakT'] = peakT
    var_results['tmax'] = tmax
    var_results['peak_coreT'] =peak_coreT
    var_results['tcoremax'] = tcoremax
    
    var_results.to_csv(path+'run_results.csv',index=False)