#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Predict whether Finmarken and Glorieta Mountain will record a dynamo based on sucessful model runs
"""
#%% Imports
import pandas as pd
import numpy as np
import sys
sys.path.append('../')
Myr = 1e6*365*24*3600 #seconds in a million years
from pallasite_dynamo_functions import dynamo_status
from pallasite_cooling_functions import find_temp_depth, find_abs_depth
from average_B import average_B_rem

# Load data
folder = sys.argv[1]
mdata = pd.read_csv(f'../Results/{folder}/pallasite_sucess_info.csv',skiprows=[1]) #sucessful model params
mdata = mdata[mdata['f3']==True] #only keep sucessful runs
edata = {'cr_yang_low': np.array([2.2, 17.5]),
         'cr_yang_up': np.array([2.8, 19.9]),} #[K/Myr] Yang et. al. 2010 cooling rates at 925K
Remc = 10 #critical magnetic Reynolds number
Xs_eutectic=33 #eutectic composition
#create dataframe to store output
pdata = pd.DataFrame(columns=['run','gm_dlow','gm_dup','f_dlow','f_dup','gm_tlow','gm_tup','f_tlow','f_tup','gm_mag','f_mag','gm_core_solid','f_core_solid'])

#loop over runs
i = 0 # index for saving to dataframe
for run in mdata['run']: 
    #navigate to correct subfolder
    if run%12 == 0:
        folder_num = int(run/12)
    else:
        folder_num = int(run/12) +1
    subfolder = f'params_{folder_num}'
    mout = pd.read_csv(f'../Results/{folder}/{subfolder}/run_results.csv',skiprows=[1]) #model outputs
    run = int(run)
    #load data
    npzfile = np.load(f'../Results/{folder}/{subfolder}/run_{run}_B.npz')
    B = npzfile['B']
    Rem = npzfile['Rem']
    Xs = npzfile['Xs']
    t = npzfile['t']/Myr
    temp = npzfile['T_profile']
    r = mdata.loc[mdata['run'] == run, 'r'].values[0]
    rcr = mdata.loc[mdata['run'] == run, 'rcr'].values[0]
    Xs0 = mdata.loc[mdata['run'] == run, 'Xs_0'].values[0]
    eta0 = mdata.loc[mdata['run'] == run, 'eta0'].values[0]
    dr = mdata.loc[mdata['run'] == run, 'dr'].values[0]
    dt = t[1] - t[0] #time step in data [Myr]
    tsolid_start = mout.loc[mout['run'] == run, 'tsolid_start'].values[0]
    #average B and Rem
    if Xs0!=Xs_eutectic: #if the core doesn't start at the eutectic composition
        B, Rem = average_B_rem(B, Rem, t, Xs, Xs_eutectic, tsolid_start)
    #calculate radial grid and truncate Tprofile
    n_cells = int(r/dr) +1 #number of cells needed to span the body including one at the centre
    i_core = round((n_cells-3)*rcr) +1 #index of CMB
    rprof = np.arange(i_core*dr,r+dr,dr)/1e3
    temp = temp[:,i_core:] #restrict to mantle cells
    #calculate dTdt
    tempdt = -np.gradient(temp,dt,axis=0) #cooling rate [K/Myr]
    
    #find contour for Yang et. al. 2010 data, find depths where cooling rates match
    rind = find_temp_depth(925,edata['cr_yang_low'],edata['cr_yang_up'], t, temp, tempdt) 
    #find depths of each sample and time of remanence acquisition
    depth, time = find_abs_depth(rind,rprof,t,temp,593)

    #find out if the dynamo is on/off at 593K for each sample
    Bmodel, c = dynamo_status(rind,t,temp,Rem,Remc,B,tsolid_start)
    

    #save to file
    pdata.loc[i] = [run, depth[0, 0], depth[0, 1], depth[1, 0], depth[1, 1], time[0,0], time[0,1], time[1,0], time[1,1], Bmodel[0], Bmodel[1],c[0],c[1]]
    i += 1

# merge with overall results
mdata['run'] = mdata['run'].astype(int) #convert to int
sucesses = pd.merge(mdata, pdata, on="run") #join together on matching runs
sucesses.to_csv(f'../Results/{folder}/gmf_sucess_info.csv',index=False,mode='a',header=False)