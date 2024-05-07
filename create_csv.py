#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to create auto_params.csv for all possible combinations of parameters
"""
import pandas as pd
import numpy as np
 
folder = 'Run_params/Fullrun2/' #folder to save parameters in   
#for parameter list
#each variable
#radius [m]
nr = 4 #number of points
minr = 100e3 #min val
maxr = 400e3 #max val
r = np.linspace(minr,maxr,nr) 

#rcmf
nrcmf = 4
minrcmf = 0.2
maxrcmf = 0.5
rcmf = np.linspace(minrcmf,maxrcmf,nrcmf)

#reference viscosity
neta0 = 7
mineta0 = 14 #10^{mineta0}
maxeta0 = 21 
eta0 = np.logspace(mineta0,maxeta0,neta0)

#beta
nbeta = 5
minbeta = 0.005
maxbeta = 0.08
beta = np.linspace(minbeta,maxbeta,nbeta)

#sulfur content
nxs = 3
minxs = 28.5 #will vary by body size later
maxxs = 32
Xs_0 = np.linspace(minxs,maxxs,nxs)

#iron content - specified manually for now
nfe = 3
minfe = 0
maxfe = 1e-7
Fe0 = np.array([0,1e-8,1e-7])

#liquid viscosity
netal = 3
minetal = 1 #log10(etal)
maxetal = 3 #log10(etal)
etal = np.logspace(minetal,maxetal,3)

#diffusion vs dislocation creep
#alpha_n = np.array([25,30]) #diffusion vs (diffusion+dislocation) creep (S&C)


#parameters fixed in this combination
default = 'vary' #viscosity profile will change
t_acc_m = 0.8 #accretion time [Myr]
t_end_m = 1000 #maximum end time of simulation [Myr]
w = 5 #width of linear region [K]
#etal = 100 #liquid viscosity [Pas]
dr = 500 # grid space [m]
dt = 0.075 #timestep [fraction of core conductive timestep]
alpha_n = 25 #melt weakening, diffusion creep

# create dictionaries and write to csv
run = 1
csv_num = 1
run_info = pd.DataFrame(columns=['run','r','default','rcmf','eta0','beta','w','etal','alpha_n','Xs_0','Fe0','t_acc_m','t_end_m','dr','dt','status']) #create columns of dataframe
unit_row = ['','m','','','Pas','K^-1','K','Pas','','wt %','60Fe/56Fe','Myr','Myr','m','t_cond_core',''] #first row is units
run_info.loc[len(run_info)] = unit_row
for n, feval in enumerate(Fe0):
    for j, rcmfval in enumerate(rcmf):
        for k, eta0val in enumerate(eta0):
            for l, betaval in enumerate(beta):
                for p, etalval in enumerate(etal):
                    if run > 1:
                        #skip this step on the first step down
                        run_info.to_csv(f'{folder}auto_params_{csv_num}.csv',index=False)
                        csv_num = csv_num + 1 #start making a new csv
                        run_info = pd.DataFrame(columns=['run','r','default','rcmf','eta0','beta','w','etal','alpha_n','Xs_0','Fe0','t_acc_m','t_end_m','dr','dt','status']) #create columns of dataframe
                        unit_row = ['','m','','','Pas','K^-1','K','Pas','','wt %','60Fe/56Fe','Myr','Myr','m','t_cond_core',''] #first row is units
                        run_info.loc[len(run_info)] = unit_row
                    for m, xsval in enumerate(Xs_0):
                        for i, rval in enumerate(r): #run r as last so all sets of files take similar time
                            run_info = pd.concat([run_info, pd.DataFrame({"run":[run],"r":[rval],"default":[default],"rcmf":[rcmfval],"eta0":[eta0val],"beta":[betaval],"w":[w],"etal":[etalval],"alpha_n":[alpha_n],"Xs_0":[xsval], "Fe0":[feval], "t_acc_m":[t_acc_m], "t_end_m":[t_end_m], "dr":[dr],"dt":[dt],"status":""})],ignore_index=True)
                            run = run+1

#call once more at the end for the final csv                        
run_info.to_csv(f'{folder}auto_params_{csv_num}.csv',index=False)                        
