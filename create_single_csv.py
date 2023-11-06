#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to create auto_params.csv for chosen parameters ranges
"""
import pandas as pd
import numpy as np
 
folder = 'Run_params/Paper_run100km/' #folder to save parameters in   
#for parameter list
#each variable

#%%
#rcmf
nrcmf = 6
minrcmf = 0.2
maxrcmf = 0.5
rcmf = np.linspace(minrcmf,maxrcmf,nrcmf)

#%%
#reference viscosity
neta0 = 11
mineta0 = 14 #10^{mineta0}
maxeta0 = 24 
eta0 = np.logspace(mineta0,maxeta0,neta0)

#%%
#beta
nbeta = 5
minbeta = 0.01 #from Dodds 2021
maxbeta = 0.035 #for E = 570 kj/mol 
beta = np.linspace(minbeta,maxbeta,nbeta)

#%%
#liquid viscosity
netal = 3
minetal = 0 #log10(etal)
maxetal = 2 #log10(etal)
etal = np.logspace(minetal,maxetal,3)

#%%
#sulfur content
nxs = 3
minxs = 27.1 #minimum for a 300km radius, rcmf = 0.3 
maxxs = 32
Xs_0 = np.linspace(minxs,maxxs,nxs)

#%%
#iron content - specified manually for now
nfe = 5
minfe = 0
maxfe = 6e-7
Fe0 = np.array([0,1e-9,1e-8,1e-7,6e-7])

#%%
#diffusion vs dislocation creep
nalpha = 5
minalpha = 25
maxalpha = 45
alpha_n = np.linspace(minalpha,maxalpha,nalpha) #diffusion vs (diffusion+dislocation) creep (S&C)

#%%
# #radius [m]
nr = 5 #number of points
minr = 100e3 #min val
maxr = 500e3 #max val
r = np.array([100e3,200e3,300e3,400e3,500e3])

#%%
# parameters that are always fixed
default = 'vary' #viscosity profile will change
t_acc_m = 0.8 #accretion time [Myr]
t_end_m = 1500 #maximum end time of simulation [Myr]
w = 5 #width of linear region [K]
dr = 500 # grid space [m]
dt = 0.075 #timestep [fraction of core conductive timestep]

#%%
csv_nums = {'rcmf':1,'eta0':2,'beta':3,'etal':4,'Xs_0':5,'Fe0':6,'alpha_n':7,'r':8}
nruns = [nrcmf,neta0,nbeta,netal,nxs,nfe,nalpha,nr] 
variables = ['rcmf','eta0','beta','etal','Xs_0','Fe0','alpha_n','r']
varlab = {'rcmf':rcmf,'eta0':eta0,'beta':beta,'etal':etal,'Xs_0':Xs_0,'Fe0':Fe0,'alpha_n':alpha_n,'r':r}

run_nums = np.zeros([len(nruns)])
for i in range(len(nruns)):
    if i==0:
        run_nums[i] = 1
    else:
        run_nums[i] = run_nums[i-1]+nruns[i-1]
        

#%%
# create dictionaries and write to csv - change these three lines

for i, runs in enumerate(run_nums):
    run = runs
    csv_num = csv_nums[variables[i]]
    var = varlab[variables[i]] #choose your variable
    
    #set constant parameters - middle value of each range
    rcmfval= 0.3
    eta0val= 1e19
    betaval = 0.0225
    etalval = 100 #liquid viscosity [Pas]
    xsval = 29.5
    feval = 1e-8
    alpha_nval = 30 #melt weakening, diffusion creep
    rval=100e3

    run_info = pd.DataFrame(columns=['run','r','default','rcmf','eta0','beta','w','etal','alpha_n','Xs_0','Fe0','t_acc_m','t_end_m','dr','dt','status']) #create columns of dataframe
    unit_row = ['','m','','','Pas','K^-1','K','Pas','','wt %','60Fe/56Fe','Myr','Myr','m','t_cond_core',''] #first row is units
    run_info.loc[len(run_info)] = unit_row
    
    
    for val in var:
        if val == var[0]: #create csv headers
            run_info = pd.DataFrame(columns=['run','r','default','rcmf','eta0','beta','w','etal','alpha_n','Xs_0','Fe0','t_acc_m','t_end_m','dr','dt','status']) #create columns of dataframe
            unit_row = ['','m','','','Pas','K^-1','K','Pas','','wt %','60Fe/56Fe','Myr','Myr','m','t_cond_core',''] #first row is units
            run_info.loc[len(run_info)] = unit_row
        #make data frame of correct length will all same values
        run_info = pd.concat([run_info, pd.DataFrame({"run":[run],"r":[rval],"default":[default],"rcmf":[rcmfval],"eta0":[eta0val],"beta":[betaval],"w":[w],"etal":[etalval],"alpha_n":[alpha_nval],"Xs_0":[xsval], "Fe0":[feval], "t_acc_m":[t_acc_m], "t_end_m":[t_end_m], "dr":[dr],"dt":[dt],"status":""})],ignore_index=True)
        run = run+1
    #replace variable values
    run_info.loc[1:,variables[i]]=var
    #create csv                      
    run_info.to_csv(f'{folder}auto_params_{csv_num}.csv',index=False)                        
