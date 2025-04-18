#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to create auto_params.csv for all possible combinations of parameters
"""
import pandas as pd
import numpy as np
 
folder = 'Run_params/Fullrun/' #folder to save parameters in   
#for parameter list
#each variable
#radius [m]
nr = 4 #number of points
minr = 200e3 #min val
maxr = 500e3 #max val
r = np.linspace(minr,maxr,nr) 

#rcmf
nrcmf = 3
minrcmf = 0.3
maxrcmf = 0.5
rcmf = np.linspace(minrcmf,maxrcmf,nrcmf)

#reference viscosity
neta0 = 5
mineta0 = 15 #10^{mineta0}
maxeta0 = 23 
eta0 = np.logspace(mineta0,maxeta0,neta0)

#beta
nbeta = 5
minbeta = 0.01
maxbeta = 0.05
beta = np.linspace(minbeta,maxbeta,nbeta)

#sulfur content
nxs = 3
minxs = 28.7 #will vary by body size later
maxxs = 32
Xs_0 = np.linspace(minxs,maxxs,nxs)

#melt weakening exponent
nalphan = 5
minalphan = 25
maxalphan = 45
alpha_n = np.linspace(minalphan,maxalphan,nalphan) 

#core radius fraction
nrc = 5
minrc = 0.1
maxrc = 0.9
rc = np.linspace(minrc,maxrc,nrc)

#water content
nxw = 4
xw = np.array([0,0.02,0.05,0.07]) #wt %

#parameters fixed in this combination
accrete = False #start at differentiation
default = 'vary' #viscosity profile will change
t_start_m = 1 #accretion time [Myr]
t_end_m = 1000 #maximum end time of simulation [Myr]
w = 5 #width of linear region [K]
etal = 100 #liquid viscosity [Pas]
Fe0 = 1e-8 #initial 60Fe/56Fe
dr = 500 # grid space [m]
dt = 0.075 #timestep [fraction of core conductive timestep]
icfrac = 1 #core solidification endmember

# create dictionaries and write to csv
run = 1
csv_num = 1
run_info = pd.DataFrame(columns=['run','r','rcr','default','rcmf','eta0','beta','w','etal','alpha_n','Xs_0','Fe0','t_start_m','t_end_m','dr','dt','icfrac','xwater','accrete','status']) #create columns of dataframe
unit_row = ['','m','','','','Pas','K^-1','K','Pas','','wt %','60Fe/56Fe','Myr','Myr','m','t_cond_core','','wt %','',''] #first row is units
run_info.loc[len(run_info)] = unit_row
for n, xwval in enumerate(xw):
    for j, rcmfval in enumerate(rcmf):
        for k, eta0val in enumerate(eta0):
            for l, betaval in enumerate(beta):
                for p, rcval in enumerate(rc):
                    for q, alphanval in enumerate(alpha_n):
                        if run > 1:
                            #skip this step on the first step down
                            run_info.to_csv(f'{folder}auto_params_{csv_num}.csv',index=False)
                            csv_num = csv_num + 1 #start making a new csv
                            run_info = pd.DataFrame(columns=['run','r','rcr','default','rcmf','eta0','beta','w','etal','alpha_n','Xs_0','Fe0','t_start_m',
                                                             't_end_m','dr','dt','icfrac','xwater','accrete','status']) #create columns of dataframe
                            unit_row = ['','m','','','','Pas','K^-1','K','Pas','','wt %','60Fe/56Fe','Myr','Myr','m','t_cond_core','','wt %','',''] #first row is units
                            run_info.loc[len(run_info)] = unit_row
                        for m, xsval in enumerate(Xs_0):
                            for i, rval in enumerate(r): #run r as last so all sets of files take similar time
                                run_info = pd.concat([run_info, pd.DataFrame({"run":[run],"r":[rval],"rcr":[rcval],"default":[default],
                                                                              "rcmf":[rcmfval],"eta0":[eta0val],"beta":[betaval],"w":[w],"etal":[etal],
                                                                              "alpha_n":[alphanval],"Xs_0":[xsval], "Fe0":[Fe0], "t_start_m":[t_start_m],
                                                                                "t_end_m":[t_end_m], "dr":[dr],"dt":[dt],"icfrac":[icfrac],"xwater":[xwval],
                                                                                "accrete":[accrete],"status":""})],ignore_index=True)
                                run = run+1

#call once more at the end for the final csv                        
run_info.to_csv(f'{folder}auto_params_{csv_num}.csv',index=False)                        
