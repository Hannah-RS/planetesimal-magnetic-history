#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate the field strength at the time of remanace acquisition for the main group pallasites.
Doesn't currently include the runs that had to be redone.
"""
#%% Imports
import pandas as pd
import numpy as np
import sys
sys.path.append('../')
Myr = 1e6*365*24*3600 #seconds in a million years
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
pdata = pd.DataFrame(columns=['run','Bav3','Bav4','Bav5'])

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
    r = mdata.loc[mdata['run'] == run, 'r'].values[0]
    rcr = mdata.loc[mdata['run'] == run, 'rcr'].values[0]
    Xs0 = mdata.loc[mdata['run'] == run, 'Xs_0'].values[0]
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
    
    Bout = np.zeros([3])

    for j in range(3,6): #loop over the three pallasite depths that record a remnanance
        tlow = mdata.loc[mdata['run'] == run, f't{j}_low'].values[0]
        tup = mdata.loc[mdata['run'] == run, f't{j}_up'].values[0]
        depth = 1e3*(mdata.loc[mdata['run'] == run, f'd{j}_low'].values[0] + mdata.loc[mdata['run'] == run, f'd{j}_up'].values[0])/2 #average depth below surface kilometres.
        #Find Bfield in that time.
        Bsurf = B[(t>=tlow) & (t<=tup)].mean()
        #Scale by distance from the surface
        Bout[j] = Bsurf*(r/(r-depth))**3
    #save to file
    pdata.loc[i] = [run, Bout[0], Bout[1], Bout[2]]
    i += 1

# merge with overall results
mdata['run'] = mdata['run'].astype(int) #convert to int
sucesses = pd.merge(mdata, pdata, on="run") #join together on matching runs
sucesses.to_csv(f'../Results/{folder}/bstrength_sucess_info.csv',index=False,mode='a',header=False)