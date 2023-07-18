#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stitch results of runs into one csv
"""
from load_info import combine_info

import sys
folder1 = sys.argv[1]
folder = 'Results_combined/'+folder1 #folder
subfolder = 'params_' #subfolder
num_runs = int(sys.argv[2]) #number of subruns
nfails = 0 #total number of fails
nsuc = 0 #total number of sucesses
for i in range(num_runs):
    #combine all run params and results
    results = combine_info(folder+subfolder+f'{(i+1)}/','auto_params.csv','run_results.csv',['MAC_onoff.csv','CIA_onoff.csv','comp_onoff.csv','coreconv_onoff.csv'])
    #count number of failed runs
    fails = results[results['status']!=1]
    nfails = nfails+len(fails)
    sucesses = results[results['status']==1]
    nsuc = nsuc + len(sucesses)
    sucesses.to_csv(f'{folder}all_sucess_info.csv',index=False,mode='a',header=False)
    fails.to_csv(f'{folder}all_fails_info.csv',index=False,mode='a',header=False)

ntot = nfails+nsuc #total number of runs    
print(f'{nsuc} runs were sucessful - {nsuc/ntot*100:.1f}%')
print(f'{nfails} runs failed - {nfails/ntot*100:.1f}%')