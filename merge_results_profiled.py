#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stitch results of runs into one csv - version with decorator to test memory usage
To see memory usage run python merge_results_profiled.py on the command line
"""
import pandas as pd
from memory_profiler import profile

folder = 'Results_combined/Test/' #folder
files = 'params_' #subfolder
num_runs = 4 #number of subruns

@profile
def combine_csv():
    for i in range(num_runs):
        #combine run info
        params = pd.read_csv(folder+files+f'{(i+1)}'+'/auto_params.csv',skiprows=[1])
        params.to_csv(f'{folder}auto_params.csv',index=False,mode='a',header=False)
        results = pd.read_csv(folder+files+f'{(i+1)}'+'/run_results.csv',skiprows=[1])
        results.to_csv(f'{folder}run_results.csv',index=False,mode='a',header=False)
        
if __name__ == '__main__':
       combine_csv()
