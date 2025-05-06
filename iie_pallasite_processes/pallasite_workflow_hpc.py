#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Find possible parameters for pallasites from a set of model runs
"""
#%% Imports
import pandas as pd
import numpy as np
import sys
sys.path.append('../')
Myr = 1e6*365*24*3600 #seconds in a million years
from pallasite_dynamo_functions import dynamo_check, paleo_check
from pallasite_cooling_functions import find_temp_depth, check_623_cool
from average_B import average_B_rem

# Load data
folder = sys.argv[1]
subfolder = sys.argv[2]
mdata = pd.read_csv(f'../Results/{folder}/{subfolder}/auto_params.csv',skiprows=[1]) #model params
mout = pd.read_csv(f'../Results/{folder}/{subfolder}/run_results.csv',skiprows=[1]) #model outputs
edata = pd.read_csv('pallasite_data.csv',skiprows=[1]) #experimental data
Remc = 10 #critical magnetic Reynolds number
Xs_eutectic=33 #eutectic composition
#add relative paleointensity column
Bmax = edata['B_mid'].max()
Bmax_err = edata.loc[edata['B_mid'] == Bmax, 'B_err'].values[0]
edata['Brel'] = edata['B_mid']/edata['B_mid'].max()
edata['Brel_err'] = np.sqrt((edata['B_err']/Bmax)**2 + (Bmax_err**2*edata['B_mid']**2/Bmax**4))
#create dataframe to store output
pdata = pd.DataFrame(columns=['run','d1_low','d1_up','d2_low','d2_up','d3_low',
                                      'd3_up','d4_low','d4_up','d5_low','d5_up','t1_low',
                                      't1_up','t2_low','t2_up','t3_low','t3_up','t4_low',
                                      't4_up','t5_low','t5_up','f1','f2','f3','f4','c1',
                                      'c2','c3','c4','c5','subfolder'])
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
    n_cells = int(r/dr) +1 #number of cells needed to span the body including one at the centre
    i_core = round((n_cells-3)*rcr) +1 #index of CMB
    rprof = np.arange(i_core*dr,r+dr,dr)/1e3
    temp = temp[:,i_core:] #restrict to mantle cells
    #calculate dTdt
    tempdt = -np.gradient(temp,dt,axis=0) #cooling rate [K/Myr]
    # check if there is a gap in dynamo generation, f1 = False if not, 
    if mout.loc[mout['run'] == run, 'Bn1'].values[0] > 1:
        f1 = True
    else:
        f1 = False
    
    if f1 == True:
        #find contour for Yang et. al. 2010 data, find depths where cooling rates match
        rind = find_temp_depth(925,edata['cr_yang_low'],edata['cr_yang_up'], t, temp, tempdt)
        #check if Maurel cooling rates match at these positions,refine depths
        f2, rind = check_623_cool(rind,edata['cr_623_low'],edata['cr_623_up'],t, temp, tempdt)
    else: 
        f2 = False
    
    # check if dynamo is on/off at 593K at these positions for the correct samples
    if f2 == True:
        f3, depth, time, Bmodel, c = dynamo_check(edata['B_mid'],rind,t,temp,Rem,Remc,B,tsolid_start,rprof)
    else: 
        f3 = False
        depth = np.zeros([5,2])
        time = np.zeros([5,2])
        c = [False,False,False,False,False]

    #check relative paleointensities of experimental and model data, f4 = False if not
    if f3 == True:
        f4 = paleo_check(Bmodel, edata['Brel'], edata['Brel_err'])   
    else:
        f4 = False

    #save to file
    pdata.loc[i] = [run, depth[0, 0], depth[0, 1], depth[1, 0], depth[1, 1], depth[2, 0], depth[2, 1],
                depth[3,0], depth[3,1], depth[4,0], depth[4,1], time[0,0], time[0,1], time[1,0], time[1,1], time[2,0], 
                time[2,1], time[3,0], time[3,1], time[4,0], time[4,1], f1, f2, f3, f4, c[0], c[1], c[2], c[3], c[4],subfolder[7]]
    i += 1


# merge with overall results
mdata['run'] = mdata['run'].astype(int) #convert to int
sucesses = pd.merge(mdata, pdata, on="run") #join together on matching runs
sucesses.to_csv(f'../Results/{folder}/pallasite_sucess_info.csv',index=False,mode='a',header=False)