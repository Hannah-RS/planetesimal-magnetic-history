#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stitch results of runs into one csv - version with decorator to test memory usage
To see memory usage run python merge_results_profiled.py on the command line
"""
from load_info import combine_info
from memory_profiler import profile
import pandas as pd
folder = 'Results_combined/Fullrun1/' #folder
subfolder = 'params_' #subfolder
num_runs = 3 #number of subruns


@profile
def combine_csv(): 
    nfails = 0 #total number of fails
    nsuc = 0 #total number of sucesses
    ninval = 0 #number of invalid parameter combos
    for i in range(num_runs):
        #combine all run params and results
        # by default the merging here will remove unfinished runs which didn't save mag field parameters
        results = combine_info(folder+subfolder+f'{(i+1)}/','auto_params.csv','run_results.csv')   
        sucesses = results[results['status']==1]
        nsuc = nsuc + len(sucesses)
        sucesses.to_csv(f'{folder}all_sucess_info.csv',index=False,mode='a',header=False)
        #count number of failed runs and combine parameters for rerun
        run_info = pd.read_csv(folder+subfolder+f'{(i+1)}/auto_params.csv',delimiter=',',skiprows=[1]) 
        fails = run_info[run_info['status']==0]
        nfails = nfails+len(fails)
        fails.to_csv(f'{folder}fail_params.csv',index=False,mode='a',header=False)
        #count number of invalid viscosity profiles
        inval = run_info[run_info['status']==-1]
        ninval = ninval+len(inval)
        inval.to_csv(f'{folder}inval_params.csv',index=False,mode='a',header=False)
    ntot = nfails+nsuc+ninval #total number of runs    
    print(f'{nsuc} runs were sucessful - {nsuc/ntot*100:.1f}%')
    print(f'{nfails} runs failed - {nfails/ntot*100:.1f}%')
    print(f'{ninval} runs had invalid viscosity parameters - {ninval/ntot*100:.1f}%')

        
if __name__ == '__main__':
       combine_csv()
