#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for post processing magnetic field duration data
"""
import numpy as np
import pandas as pd
import sys
sys.path.append('../')
from duration_calc import on_off_test
folder = 'Singlevar1/'
subfolders = {'rcmf':1,'eta0':2,'beta':3,'etal':4,'Xs_0':5,'Fe0':6,'alpha_n':7,'r':8}
labels = {'rcmf':'$\\phi_{{RCMF}}$','eta0':'$\\eta_0$','beta':'beta','etal':'$\\eta_l$','Xs_0':'$X_{{s,0}}$','Fe0':'$^{{60}}Fe/^{{56}}Fe$','alpha_n':'$\\alpha_n$','r':'radius'}

from parameters import save_interval_t, Myr, f0

for j in range(8):

    path = '../Results_combined/'+folder+f"params_{j+1}/"
    
    var_data = pd.read_csv(path+'auto_params.csv',skiprows=[1])
    var_results = pd.read_csv(path+'run_results.csv',skiprows=[1])
    # #drop uneeded columns
    var_results.drop(columns=list(var_results.columns.values[14:]),inplace=True)
    minrun = min(var_data['run'])
    maxrun = max(var_data['run'])
    nrun = len(var_data)
    
    Bn1 = np.zeros([nrun])
    Bn2 = np.zeros([nrun])
    Bn3 = np.zeros([nrun])
    max_B = np.zeros([nrun])
    max_Bt = np.zeros([nrun])
    max_R = np.zeros([nrun])
    max_Rt = np.zeros([nrun])
    magon_1 = np.zeros([nrun])
    magoff_1 = np.zeros([nrun])
    magon_2 = np.zeros([nrun])
    magoff_2 = np.zeros([nrun])
    magon_3 = np.zeros([nrun])
    magoff_3 = np.zeros([nrun])
    magon_4 = np.zeros([nrun])
    magoff_4 = np.zeros([nrun])
    magon_5 = np.zeros([nrun])
    magoff_5 = np.zeros([nrun])
    magon_6 = np.zeros([nrun])
    magoff_6 = np.zeros([nrun])
    tsolid_start = np.zeros([nrun])
    
    for i in range(nrun):
        run = int(minrun + i)
        rc = var_data.loc[var_data['run']==run,'r'].values[0]/2
        npzfile = np.load(f'{path}run_{run}.npz')
        f = npzfile['f']
        t = npzfile['t']/Myr
        tsolid_start[i] = t[f<f0][0]
        #calculate dynamo generation times
        therm = npzfile['therm']
        comp = npzfile['comp']
        Rem = np.concatenate([therm[0,:],comp[0,:]])
        B = np.concatenate([therm[1,:],comp[1,:]])
        t = npzfile['t']/Myr
            
        ################# Threshold for magnetic field being on ########################
        threshold1 = 10 
        threshold2 = 40
        threshold3 = 100
    
        ########### Dynamo maxima - must be greater than threshold1 to be on ##########
        #thermal dynamo
        if np.any(Rem>threshold1):
            max_B[i] = max(B[Rem>threshold1])
            max_Bt[i] = t[B==max_B[i]][0]
        else:
            max_B[i] = 0
            max_Bt[i] = np.nan
            
        max_R[i] = max(Rem)
        max_Rt[i] = t[Rem==max_R[i]][0]
    
        ########################## on and off times - calculate and save ####################
    
        #Rem > 10  
        on, off, dur = on_off_test(t,Rem,threshold1,100*save_interval_t/Myr) #use 10 Myr interval to split up dynamo generation periods
        Bn1[i] = len(on) #number of on periods
        if len(np.where(on>0)) > 0:
            magon_1[i] = on[0]
            magoff_1[i] = off[0]
        else:
            magon_1[i] = 0
            magoff_1[i] = 0    
        if len(on) > 1:
            magon_2[i] = on[1]
            magoff_2[i] = off[1]
        else:
            magon_2[i] = 0
            magoff_2[i] = 0
    
        #Rem > 40
        on, off, dur = on_off_test(t,Rem,threshold2,100*save_interval_t/Myr) #use 10 Myr interval to split up dynamo generation periods
        Bn2[i] = len(on) #number of on periods
        if len(np.where(on>0)) > 0: #i.e. there is one on value > nan
            magon_3[i] = on[0]
            magoff_3[i] = off[0]
        else:
            magon_3[i] = 0
            magoff_3[i] = 0
            
        if len(on) > 1:
            magon_4[i] = on[1]
            magoff_4[i] = off[1]
        else:
            magon_4[i] = 0
            magoff_4[i] = 0
            
        # Rem > 100
        on, off, dur = on_off_test(t,Rem,threshold3,100*save_interval_t/Myr) #use 10 Myr interval to split up dynamo generation periods
        Bn3[i] = len(on) #number of on periods
        if len(np.where(on>0)) > 0:
            magon_5[i] = on[0]
            magoff_5[i] = off[0]
        else:
            magon_5[i] = 0
            magoff_5[i] = 0
            
        if len(on) > 1:
            magon_6[i] = on[1]
            magoff_6[i] = off[1]
        else:
            magon_6[i] = 0
            magoff_6[i] = 0

    var_results['max_R'] = max_R
    var_results['max_Rt'] = max_Rt
    var_results['max_B'] = max_B
    var_results['max_Bt'] = max_Bt
    var_results['Bn1'] = Bn1
    var_results['magon_1'] = magon_1
    var_results['magoff_1'] = magoff_1
    var_results['magon_2'] =magon_2
    var_results['magoff_2'] = magoff_2
    var_results['Bn2'] = Bn2
    var_results['magon_3'] = magon_3
    var_results['magoff_3'] = magoff_3
    var_results['magon_4'] =magon_4
    var_results['magoff_4'] = magoff_4
    var_results['Bn3'] = Bn3
    var_results['magon_5'] = magon_5
    var_results['magoff_5'] = magoff_5
    var_results['magon_6'] =magon_6
    var_results['magoff_6'] = magoff_6

    var_results.insert(14, "tsolid_start", tsolid_start)
    #add a unit row
    unit_row = ['','Myr','s','Myr','K','K','Myr','K','Myr','Myr','Myr','Myr','Myr','K','Myr','','Myr','T','Myr','','Myr','Myr','Myr','Myr','','Myr','Myr','Myr','Myr','','Myr','Myr','Myr','Myr'] #first row is units
    var_results.loc[-1] = unit_row
    var_results.sort_index(inplace=True)
    var_results.to_csv(path+'run_results.csv',index=False)