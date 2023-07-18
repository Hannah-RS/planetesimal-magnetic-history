#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stitch results of runs into one csv
"""
from load_info import combine_info
import pandas as pd

import sys
folder1 = sys.argv[1]
folder = 'Results_combined/'+folder1 #folder
subfolder = 'params_' #subfolder
num_runs = int(sys.argv[2]) #number of subruns
nfails = 0 #total number of fails
nsuc = 0 #total number of sucesses
for i in range(num_runs):
    #combine all run params and results
    # by default the merging here will remove unfinished runs which didn't save mag field parameters
    results = combine_info(folder+subfolder+f'{(i+1)}/','auto_params.csv','run_results.csv',['MAC_onoff.csv','CIA_onoff.csv','comp_onoff.csv','coreconv_onoff.csv'])   
    sucesses = results[results['status']==1]
    nsuc = nsuc + len(sucesses)
    sucesses.to_csv(f'{folder}all_sucess_info.csv',index=False,mode='a',header=False)
    #count number of failed runs and combine parameters for rerun
    run_info = pd.read_csv(folder+subfolder+f'{(i+1)}/auto_params.csv',delimiter=',',skiprows=[1]) 
    fails = run_info[run_info['status']!=1]
    nfails = nfails+len(fails)
    fails.to_csv(f'{folder}fail_params.csv',index=False,mode='a',header=False)

ntot = nfails+nsuc #total number of runs    
print(f'{nsuc} runs were sucessful - {nsuc/ntot*100:.1f}%')
print(f'{nfails} runs failed - {nfails/ntot*100:.1f}%')