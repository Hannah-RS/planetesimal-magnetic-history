#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stitch results of runs into one csv
"""
import pandas as pd

folder = 'Results_combined/Test/' #folder
files = 'params_' #subfolder
num_runs = 4 #number of subruns

for i in range(num_runs):
    #combine run info
    params = pd.read_csv(folder+files+f'{(i+1)}'+'/auto_params.csv',skiprows=[1])
    params.to_csv(f'{folder}auto_params.csv',index=False,mode='a',header=False)
    results = pd.read_csv(folder+files+f'{(i+1)}'+'/run_results.csv',skiprows=[1])
    results.to_csv(f'{folder}run_results.csv',index=False,mode='a',header=False)
