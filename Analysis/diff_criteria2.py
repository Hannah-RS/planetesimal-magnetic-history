#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check if Ra becomes supercritical before d0<R in undifferentiated planetesimal
Save run parameters for which this isn't true
"""
#%% import modules
import pandas as pd
import numpy as np
import sys
sys.path.append('../Plotting_scripts/')
from plot_params import variables, subfolders
folders = ['Paper_run100km','Paper_run200km','Paper_run300km','Paper_run400km','Paper_run500km']

#%% Run loop
#bad_params = pd.DataFrame()
bad_runs = []
conv_runs = []
for j, folder in enumerate(folders):
    for var in variables[:-1]:
        path = '../Results_combined/'+folder+f"/params_{subfolders[var]}/"
        
        #find run numbers
        var_data = pd.read_csv(path+'auto_params.csv',skiprows=[1])
        minrun = min(var_data['run'])
        maxrun = max(var_data['run'])
        nrun = len(var_data)
    
        for i in range(nrun):
            run = int(minrun+i)
            #import data
            npzfile = np.load(f'{path}run_{run}_diff.npz')
            d0 = npzfile['d0']
            Ra = npzfile['Ra']
            Rac = npzfile['Ra_crit']
            t = npzfile['t_diff']
            
            conv_on = t[d0<0.99*(j+1)*100e3] #lid thickness less than radius
            Ra_super = t[Ra>Rac]
            if (len(Ra_super)>0)&(len(conv_on)>0):
                if Ra_super[0]>conv_on[0]:
                    bad_runs.append(f'{run} and size {j}')
            #how many runs have convection before differentiation
            if len(conv_on)>0:
                conv_runs.append(f'{run} and size {j}')       
#%% Print list length
print(f'Ra supercritical after lid thinner than domain for {len(bad_runs)} runs')      
print(f'{len(conv_runs)} runs ({len(conv_runs)*100/195:.0f} %) convect before differentiation')      