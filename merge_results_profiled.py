#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stitch results of runs into one csv - version with decorator to test memory usage
To see memory usage run python merge_results_profiled.py on the command line
"""
from load_info import combine_info
from memory_profiler import profile

folder = 'Results_combined/Test/' #folder
subfolder = 'params_' #subfolder
num_runs = 2 #number of subruns

@profile
def combine_csv():
    
    for i in range(num_runs):
        #combine all run params and results
        results = combine_info(folder+subfolder+f'{(i+1)}/','auto_params.csv','run_results.csv',['MAC_onoff.csv','CIA_onoff.csv','comp_onoff.csv','coreconv_onoff.csv'])
        results.to_csv(f'{folder}all_info.csv',index=False,mode='a',header=False)

        
if __name__ == '__main__':
       combine_csv()
