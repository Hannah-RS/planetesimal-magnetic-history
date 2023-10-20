#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to create auto_params.csv for chosen parameters ranges
"""
import pandas as pd
import numpy as np
 
folder = 'Run_params/Paper_run1/' #folder to save parameters in   
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
minxs = 27 #minimum for a 300km radius, rcmf = 0.3 
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

#parameters fixed in this combination - take middle value of each range
rcmfval= 0.3
eta0val= 1e19
betaval = 0.0225
etalval = 100 #liquid viscosity [Pas]
xsval = 29.5
feval = 1e-8
alpha_nval = 30 #melt weakening, diffusion creep
#rval=300e3

#%%
# create dictionaries and write to csv - change these three lines
run = 41
csv_num = 8
var = r #choose your variable

run_info = pd.DataFrame(columns=['run','r','default','rcmf','eta0','beta','w','etal','alpha_n','Xs_0','Fe0','t_acc_m','t_end_m','dr','dt','status']) #create columns of dataframe
unit_row = ['','m','','','Pas','K^-1','K','Pas','','wt %','60Fe/56Fe','Myr','Myr','m','t_cond_core',''] #first row is units
run_info.loc[len(run_info)] = unit_row


for val in var:
    if val == var[0]: #create csv headers
        run_info = pd.DataFrame(columns=['run','r','default','rcmf','eta0','beta','w','etal','alpha_n','Xs_0','Fe0','t_acc_m','t_end_m','dr','dt','status']) #create columns of dataframe
        unit_row = ['','m','','','Pas','K^-1','K','Pas','','wt %','60Fe/56Fe','Myr','Myr','m','t_cond_core',''] #first row is units
        run_info.loc[len(run_info)] = unit_row
    #change this line
    rval = val #assign variable to val for addition to data frame
    run_info = pd.concat([run_info, pd.DataFrame({"run":[run],"r":[rval],"default":[default],"rcmf":[rcmfval],"eta0":[eta0val],"beta":[betaval],"w":[w],"etal":[etalval],"alpha_n":[alpha_nval],"Xs_0":[xsval], "Fe0":[feval], "t_acc_m":[t_acc_m], "t_end_m":[t_end_m], "dr":[dr],"dt":[dt],"status":""})],ignore_index=True)
    run = run+1

#create csv                      
run_info.to_csv(f'{folder}auto_params_{csv_num}.csv',index=False)                        
