#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Find possible parameters for Eagle Station pallasites from a set of model runs
"""
#%% Imports
import pandas as pd
import numpy as np
import sys
sys.path.append('../')
Myr = 1e6*365*24*3600 #seconds in a million years
from pallasite_cooling_functions import find_temp_depth, find_abs_depth
from pallasite_dynamo_functions import dynamo_status
from average_B import average_B_rem

# Load data
folder = sys.argv[1]
subfolder = sys.argv[2]
mdata = pd.read_csv(f'../Results/{folder}/{subfolder}/auto_params.csv',skiprows=[1]) #model params
mout = pd.read_csv(f'../Results/{folder}/{subfolder}/run_results.csv',skiprows=[1]) #model outputs
edata = pd.read_csv('es_data.csv',skiprows=[1]) #experimental data
Remc = 10 #critical magnetic Reynolds number
Xs_eutectic=33 #eutectic composition
#create dataframe to store output
pdata = pd.DataFrame(columns=['run','d1_low','d1_up','d2_low','d2_up','d3_low',
                                      'd3_up','d4_low','d4_up','t1_low',
                                      't1_up','t2_low','t2_up','t3_low','t3_up','t4_low',
                                      't4_up','c1',
                                      'c2','c3','c4','B1','B2','B3','B4','f','subfolder'])
#%% Loop over runs
i = 0 # index for saving to dataframe
for run in mout['run']:
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
    n_cells = int(r/dr) +1 #number of cells needed to span the bodyenumerate(rind[:,0]) including one at the centre
    i_core = round((n_cells-3)*rcr) +1 #index of CMB
    rprof = np.arange(i_core*dr,r+dr,dr)/1e3
    temp = temp[:,i_core:] #restrict to mantle cells
    #calculate dTdt
    tempdt = -np.gradient(temp,dt,axis=0) #cooling rate [K/Myr]
    #find contour for Fish cooling rate data, find depths where cooling rates match
    rind = find_temp_depth(623,edata['cr_623_low'],edata['cr_623_up'], t, temp, tempdt)
    if np.any(rind==0): #if any meteorite doesn't reach the required cooling rate
        f = 0
        depth = np.zeros([4,2])
        time = np.zeros([4,2])
        Bmodel = np.zeros([4])
        c = [False,False,False,False]
    else:
        f = 1
        #find out if the dynamo is on/off at 593K for each sample
        Bmodel, c = dynamo_status(rind,t,temp,Rem,Remc,B,tsolid_start)
        #find depths of each sample
        depth, time = find_abs_depth(rind,rprof,t,temp,623)
    #save to file
    pdata.loc[i] = [run, depth[0, 0], depth[0, 1], depth[1, 0], depth[1, 1], depth[2, 0], depth[2, 1],
                depth[3,0], depth[3,1], time[0,0], time[0,1], time[1,0], time[1,1], time[2,0], 
                time[2,1], time[3,0], time[3,1], c[0], c[1], c[2], c[3],Bmodel[0],Bmodel[1],Bmodel[2],Bmodel[3],f,subfolder[7]]
    i += 1


# merge with overall results
mdata['run'] = mdata['run'].astype(int) #convert to int
sucesses = pd.merge(mdata, pdata, on="run") #join together on matching runs
sucesses.to_csv(f'../Results/{folder}/es_sucess_info.csv',index=False,mode='a',header=False)