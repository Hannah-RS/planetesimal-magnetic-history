#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to find maximum field strength for each body size
"""
import numpy as np
import pandas as pd


#%% Load data for chosen variable
folders = ['Paper_run100km/','Paper_run200km/','Paper_run4/','Paper_run400km/','Paper_run500km/']
peakB =[]

for folder in folders:
    
    var_results = pd.read_csv(f'../Results_combined/{folder}/all_sucess_info.csv',skiprows=[1])
    if folder == 'Paper_run4/': #exclude radius runs
        peakB.append(max(var_results.loc[var_results['r']==300000,'max_B'])/1e-6)
    else:
        peakB.append(max(var_results['max_B'])/1e-6)
    
print(np.round(peakB,0))