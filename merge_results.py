#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stitch results of runs into one csv
"""
from load_info import combine_info

folder = 'Results_combined/Test/' #folder
subfolder = 'params_' #subfolder
num_runs = 2 #number of subruns

for i in range(num_runs):
    #combine all run params and results
    results = combine_info(folder+subfolder+f'{(i+1)}/','auto_params.csv','run_results.csv',['MAC_onoff.csv','CIA_onoff.csv','comp_onoff.csv','coreconv_onoff.csv'])
    results.to_csv(f'{folder}all_info.csv',index=False,mode='a',header=False)
