#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create auto_params.csv for paramter space exploration in rcr and Xw.  
"""
import pandas as pd
import numpy as np
 
folder = 'Run_params/nc_cc_natast/' #folder to save parameters in  
 
#for parameter list
#each variable

#%%
#rcr
nrcr = 9
minrcr = 0.1
maxrcr = 0.9
rcr = np.linspace(minrcr,maxrcr,nrcr)
#Xw
xwater = np.array([0,0.02,0.03,0.04,0.05,0.06, 0.07])  #water content of mantle [wt %]
nxwater = len(xwater)
#radius [m]
r = np.array([100e3,300e3,500e3])
nr = len(r)
#%% constant parameters middle values for Sanderson et. al. 2024 unless specified
default = 'vary' 
rcmf = 0.5 #critical melt fraction - max range of S
eta0 = 1e19 #reference viscosity [Pas]
w = 5 #width of linear region [K]
beta = 0.0225 #Arrhenius slope
etal = 10 #liquid viscosity [Pas]
alpha_n = 30 #melt weakening
Xs0 = 23 #lowest value valid for whole parameter space
Fe0 = 1e-8 #initial 60Fe/56Fe ratio in accreting material
t_start_m = 2 #differentiation time [Myr]
t_end_m = 1500 #maximum end time of simulation [Myr]
icfrac = 1 #fraction of solidified mass that ends up in inner core
dr = 500 # grid space [m]
dt = 0.075 #timestep [fraction of core conductive timestep]
accrete = False #start from differentiation
#%%
csv_num = 1
for i, rval  in enumerate(r):
    for j, rcrval in enumerate(rcr):
        run_info = pd.DataFrame(columns=['run','r','rcr','default','rcmf','eta0','beta','w','etal','alpha_n','Xs_0','Fe0','t_start_m','t_end_m','dr','dt','icfrac','xwater','accrete','status']) #create columns of dataframe
        unit_row = ['','m','','','','Pas','K^-1','K','Pas','','wt %','60Fe/56Fe','Myr','Myr','m','t_cond_core','','wt %','',''] #first row is units
        run_info.loc[len(run_info)] = unit_row
        for k, xwval in enumerate(xwater):
            run = i*nr*nxwater + j*nxwater + k +1 #calculate run number
            run_info = pd.concat([run_info, pd.DataFrame({"run":[run],"r":[rval],"rcr":[rcrval],"default":[default],"rcmf":[rcmf],"eta0":[eta0],
                                                          "beta":[beta],"w":[w],"etal":[etal],"alpha_n":[alpha_n],"Xs_0":[""], "Fe0":[Fe0],
                                                            "t_start_m":[t_start_m], "t_end_m":[t_end_m], "dr":[dr],"dt":[dt],"icfrac":[icfrac],
                                                            "xwater":[xwval],"accrete":accrete,"status":""})],ignore_index=True)
        #create csv                     
        run_info.to_csv(f'{folder}auto_params_{csv_num}.csv',index=False)
        csv_num = csv_num + 1 #start making a new csv
