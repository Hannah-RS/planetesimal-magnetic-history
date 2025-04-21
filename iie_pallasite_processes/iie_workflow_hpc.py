#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Find possible parameters for IIE irons from a set of model runs on hpc
"""
#%% Imports
import pandas as pd
import numpy as np
from csv import writer
from iie_dynamo_functions import dynamo_check, paleo_check
from iie_cooling_functions import find_593_depth, find_623_cool
import sys
sys.path.append('../')
from average_B import average_B_rem
Myr = 1e6*365*24*3600 #seconds in a million years


# Load data
folder = sys.argv[1]
subfolder = sys.argv[2]
mdata = pd.read_csv(f'../Results/{folder}/{subfolder}/auto_params.csv',skiprows=[1]) #model params
mout = pd.read_csv(f'../Results/{folder}/{subfolder}/run_results.csv',skiprows=[1]) #model outputs
edata = pd.read_csv('iie_data.csv',skiprows=[1]) #experimental data
Remc = 10 #critical magnetic Reynolds number
Xs_eutectic=33 #eutectic composition
#add relative paleointensity column
Bmax = edata['B_mid'].max()
Bmax_err = edata.loc[edata['B_mid'] == Bmax, 'B_err'].values[0]
edata['Brel'] = edata['B_mid']/edata['B_mid'].max()
edata['Brel_err'] = np.sqrt((edata['B_err']/Bmax)**2 + (Bmax_err**2*edata['B_mid']**2/Bmax**4))
#%% Loop over runs
for run in mout['run']: #only do for sucessful runs
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
    nmcells = n_cells - round((n_cells-3)*rcr) -1 #number of mantle cells
    rprof = np.linspace(rcr*r,int(r),nmcells)/1e3
    temp = temp[:,(nmcells-1):] #restrict to mantle cells
    #calculate dTdt
    tempdt = -np.gradient(temp,dt,axis=0) #cooling rate [K/Myr]
    # check if dynamo is on at the radiogenic times, f1 = False if not, 
    # return c1-3 for whether core is solidifying at each time
    f1, c = dynamo_check(edata['Ar_age_low'],edata['Ar_age_up'], t, Rem, Remc, tsolid_start)

    #check relative paleointensities of experimental and model data, f3 = False if not
    if f1 == True:
        f3 = paleo_check(edata['Ar_age_low'],edata['Ar_age_up'], t, B, edata['Brel'], edata['Brel_err'])
        # find 593K contour and depths for each radiometric date
        rind = find_593_depth(edata['Ar_age_low'],edata['Ar_age_up'], t, temp)
        #find 623K cooling rate and check if it matches data, f2 = False if not
        f2, depth = find_623_cool(rind, temp, tempdt, edata['cr_623_low'], edata['cr_623_up'], rprof)    
    else:
        f2 = False
        f3 = False
        depth = np.zeros([3,2])

    #save to file
    var_list = [run, r, rcr, Xs0, eta0, depth[0, 0], depth[0, 1], depth[1, 0], depth[1, 1], depth[2, 0], depth[2, 1], f1, f2, f3, c[0], c[1], c[2],subfolder[7]]

    with open(f'../Results/{folder}/{subfolder}/iie_results.csv','a') as f_object:
        writer_object = writer(f_object) #pass file object to csv.writer
        writer_object.writerow(var_list) # pass list as argument into write row
        f_object.close() #close file

# merge with overall results
iiedata = pd.read_csv(f'../Results/{folder}/{subfolder}/iie_results.csv',skiprows=[1]) #model outputs
mdata['run'] = mdata['run'].astype(int) #convert to int
sucesses = pd.merge(mdata, iiedata, on="run") #join together on matching runs
#drop duplicate columns
sucesses = sucesses.drop(columns = ['r_y','rcr_y', 'Xs0', 'eta0_y'])
sucesses = sucesses.rename(columns={'r_x':'r','rcr_x':'rcr','eta0_x':'eta0'}) #rename columns
sucesses.to_csv(f'../Results/{folder}/iie_sucess_info.csv',index=False,mode='a',header=False)
#count number of failed runs and combine parameters for rerun
fails = mdata[mdata['status']==0]
fails['subfolder'] = subfolder[7] #just save folder number
fails.to_csv(f'../Results/{folder}/iiepallasite_fail_params.csv',index=False,mode='a',header=False)
