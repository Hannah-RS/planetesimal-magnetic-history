#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to find maximum field strength for each body size and max and min
durations for the first epoch of dynamo generation
"""
import numpy as np
import pandas as pd


#%% Load data for chosen variable
folders = ['Paper_run100km/','Paper_run200km/','Paper_run300km/','Paper_run400km/','Paper_run500km/']
peakBtherm =[]
peakBcomp =[]
max_dur = []
min_dur = []

fefilter = False #do you only want 1e-8 60Fe/56Fe

for i, folder in enumerate(folders):
    
    var_results = pd.read_csv(f'../Results_combined/{folder}/all_sucess_info.csv',skiprows=[1])
    if fefilter == True:
        var_results = var_results[var_results['Fe0']==1e-8]
    if folder == 'Paper_run300km/': #exclude radius runs
        peakBtherm.append(max(var_results.loc[var_results['r']==300000,'max_Btherm'])/1e-6)
        peakBcomp.append(max(var_results.loc[var_results['r']==300000,'max_Bcomp'])/1e-6)
        dur = var_results.loc[var_results['r']==300000,'magoff_1'] - var_results.loc[var_results['r']==300000,'magon_1']
        max_dur.append(max(dur))
        min_dur.append(min(dur[dur>0]))
    else:
        peakBtherm.append(max(var_results['max_Btherm'])/1e-6)
        peakBcomp.append(max(var_results['max_Bcomp'])/1e-6)
        dur = var_results['magoff_1'] - var_results['magon_1']
        max_dur.append(max(dur))
        min_dur.append(min(dur[dur>0]))
print('Peak thermal field strength is',np.round(peakBtherm,0))
print('Peak compositional field strength is',np.round(peakBcomp,0))
print('Minimum non-zero field duration is',np.round(min_dur,0))
print('Maximum field duration is',np.round(max_dur,0))