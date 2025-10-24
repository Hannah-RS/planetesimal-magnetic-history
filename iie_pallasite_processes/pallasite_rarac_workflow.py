#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Find possible parameters for pallasites from a set of model runs but applying a critical Ra/Rac rather than Rem
"""
#%% Imports
import pandas as pd
import numpy as np
import sys
sys.path.append('../')
Myr = 1e6*365*24*3600 #seconds in a million years
from pallasite_cooling_functions import find_temp_depth, check_623_cool
from crit_Ra_funcs import Ra_calc, Rac_calc, rarac_check, omega, G
from fe_fes_liquidus import fe_fes_density

# Load data
folder = sys.argv[1]
subfolder = sys.argv[2]
mdata = pd.read_csv(f'../Results/{folder}/{subfolder}/auto_params.csv',skiprows=[1]) #model params
mout = pd.read_csv(f'../Results/{folder}/{subfolder}/run_results.csv',skiprows=[1]) #model outputs
edata = pd.read_csv('pallasite_data.csv',skiprows=[1]) #experimental data
Remc = 10 #critical magnetic Reynolds number
Xs_eutectic=33 #eutectic composition
#create dataframe to store output
pdata = pd.DataFrame(columns=['run','d1_low','d1_up','d2_low','d2_up','d3_low',
                                      'd3_up','d4_low','d4_up','d5_low','d5_up','t1_low',
                                      't1_up','t2_low','t2_up','t3_low','t3_up','t4_low',
                                      't4_up','t5_low','t5_up','f1','f2','f3','c1',
                                      'c2','c3','c4','c5','subfolder'])
#%% Loop over runs
i = 0 # index for saving to dataframe
for run in mout['run']:
    run = int(run)
    #load data
    npzfile = np.load(f'../Results/{folder}/{subfolder}/run_{run}_B.npz')
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
    ri = rcr*r*(Xs0/Xs)**(1/3) #liquid inner core radius
    rhoc = fe_fes_density(Xs0) 
    gc = 4/3*np.pi*rcr*r*rhoc*G #gravitational field strength at CMB [m/s^2]
    #calculate radial grid and truncate Tprofile
    n_cells = int(r/dr) +1 #number of cells needed to span the body including one at the centre
    i_core = round((n_cells-3)*rcr) +1 #index of CMB
    rprof = np.arange(i_core*dr,r+dr,dr)/1e3
    temp = temp[:,i_core:] #restrict to mantle cells
    #calculate dTdt
    tempdt = -np.gradient(temp,dt,axis=0) #cooling rate [K/Myr]
    #calculate CMB dTdr
    tempdr = np.gradient(temp,dr,axis=1) #temperature gradient [K/m]
    dTdr_cmb = tempdr[:,0] #CMB temperature gradient [K/m]
    #calculate core Ra and Rac
    Rac = Rac_calc(omega, ri)
    Ra = Ra_calc(dTdr_cmb, ri, gc)
    racrit = 1e4 #critical Ra/Rac value
    # check if there is a gap in dynamo generation, f1 = False if not, 
    if (np.any(Ra/Rac < racrit))&(np.any(Ra/Rac > racrit)): #check there are values above and below critical
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
        f3, depth, time, c = rarac_check(edata['B_mid'],rind,t,temp,Ra/Rac,racrit,tsolid_start,rprof)
    else: 
        f3 = False
        depth = np.zeros([5,2])
        time = np.zeros([5,2])
        c = [False,False,False,False,False]

    #save to file
    pdata.loc[i] = [run, depth[0, 0], depth[0, 1], depth[1, 0], depth[1, 1], depth[2, 0], depth[2, 1],
                depth[3,0], depth[3,1], depth[4,0], depth[4,1], time[0,0], time[0,1], time[1,0], time[1,1], time[2,0], 
                time[2,1], time[3,0], time[3,1], time[4,0], time[4,1], f1, f2, f3, c[0], c[1], c[2], c[3], c[4],subfolder[7]]
    i += 1


# merge with overall results
mdata['run'] = mdata['run'].astype(int) #convert to int
sucesses = pd.merge(mdata, pdata, on="run") #join together on matching runs
sucesses.to_csv(f'../Results/{folder}/rarac_pallasite_sucess_info.csv',index=False,mode='a',header=False)