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
from fe_fes_liquidus import fe_fes_density
folder = 'Singlevar2/'
subfolders = {'rcmf':1,'eta0':2,'beta':3,'etal':4,'Xs_0':5,'Fe0':6,'alpha_n':7,'r':8}
labels = {'rcmf':'$\\phi_{{RCMF}}$','eta0':'$\\eta_0$','beta':'beta','etal':'$\\eta_l$','Xs_0':'$X_{{s,0}}$','Fe0':'$^{{60}}Fe/^{{56}}Fe$','alpha_n':'$\\alpha_n$','r':'radius'}
threshold = 10
from parameters import save_interval_t, Myr, alpha_c, cpc, rhofe_s, rho_exp, f0

for j in range(8):

    path = '../Results_combined/'+folder+f"params_{j+1}/"
    
    var_data = pd.read_csv(path+'auto_params.csv',skiprows=[1])
    var_results = pd.read_csv(path+'run_results.csv',skiprows=[1])
    # #drop bad columns
    # var_results.drop(columns=list(var_results.columns.values[22:]),inplace=True)
    minrun = min(var_data['run'])
    maxrun = max(var_data['run'])
    nrun = len(var_data)
    
    # magon_1 = np.zeros([nrun])
    # magoff_1 = np.zeros([nrun])
    # magbuoy_1 = np.zeros([nrun])
    # magon_2 = np.zeros([nrun])
    # magoff_2 = np.zeros([nrun])
    # magbuoy_2 = np.zeros([nrun])
    tsolid_start = np.zeros([nrun])
    
    for i in range(nrun):
        run = int(minrun + i)
        rc = var_data.loc[var_data['run']==run,'r'].values[0]/2
        npzfile = np.load(f'{path}run_{run}.npz')
        f = npzfile['f']
        t = npzfile['t']/Myr
        tsolid_start[i] = t[f<f0][0]
        #calculate dynamo generation times
    #     therm = npzfile['therm']
    #     comp = npzfile['comp']
    #     Rem = np.concatenate([therm[0,:],comp[0,:]])
    #     t = npzfile['t']/Myr
        
    #     on, off, dur = on_off_test(t,Rem,threshold,100*save_interval_t/Myr) #use 10 Myr interval to split up dynamo generation periods
    #     if len(on) > 0:
    #         magon_1[i] = on[0]
    #         magoff_1[i] = off[0]
            
    #     if len(on) > 1:
    #         magon_2[i] = on[1]
    #         magoff_2[i] = off[1]
        
    #     if len(on) > 2:
    #         print(f'Run {run} has more than two dynamo generation periods')
    
    #     #calculate average buoyancy flux ratios
    #     Flux = npzfile['Flux']
    #     f = npzfile['f']
    #     Xs = npzfile['Xs']
        
    #     therm = alpha_c*cpc*(Flux[1,:]-Flux[2,:]) #Fcmb - Fad
    #     dfdt = np.ediff1d(f)/save_interval_t
    #     #add one element in to make dfdt same length
    #     dfdt = np.insert(dfdt,0,0)
    #     rhol = fe_fes_density(Xs)*rho_exp
    #     drho = rhofe_s - rhol
    #     comp = -drho*rc*dfdt
    #     buoy = comp/therm
        
    #     #slice buoyancy along where dynamo is on
    #     on_ind = np.zeros([len(on)])
    #     off_ind = np.zeros([len(on)])
    #     buoy_av = np.zeros([len(on)])
    #     for k in range(len(on)):
    #         on_ind[k] = np.where(t==on[k])[0]
    #         off_ind[k] = np.where(t==off[k])[0]
    #         buoy_av[k] = np.average(buoy[int(on_ind[k]):int(off_ind[k])])
    #         if k == 0:
    #             magbuoy_1 = buoy_av[0]
    #         if k == 1:
    #             magbuoy_2 = buoy_av[1]
    #         if k == 2:
    #             magbuoy_2 = buoy_av[2]
                
    #     np.savez_compressed(path+f'buoy_ratio_{run}', buoy = buoy)
        
    # var_results['magon_1'] = magon_1
    # var_results['magoff_1'] = magoff_1
    # var_results['magbuoy_1'] = magbuoy_1
    # var_results['magon_2'] =magon_2
    # var_results['magoff_2'] = magoff_2
    # var_results['magbuoy_2'] = magbuoy_2
    var_results.insert(15, "tsolid_start", tsolid_start)
    #add a unit row
    unit_row = ['','Myr','s','Myr','K','K','Myr','K','Myr','Myr','Myr','Myr','Myr','K','Myr','','Myr','T','Myr','','Myr','Myr','Myr','Myr','','Myr','Myr','Myr','Myr','','Myr','Myr','Myr','Myr'] #first row is units
    var_results.loc[-1] = unit_row
    var_results.sort_index(inplace=True)
    var_results.to_csv(path+'run_results.csv',index=False)
