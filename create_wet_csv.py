#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create auto_params.csv for variations in size and differentiation time for wet and dry bodies. Creates the same number 
of csv files as accretion times.  
"""
import pandas as pd
import numpy as np
from solidus_calc import solidus
from parameters import rhom, R

folder = 'Run_params/nc_cc/' #folder to save parameters in 

#%% wet and dry params - first value is dry, second is wet
xwater = [0,0.05]  #water content of mantle [wt %]
rcr = [0.5,0.4] #core radius as a fraction of asteroid radius
eta0 = [1e19, 1e18] #reference viscosity at Tms [Pas]
eac =[440e3,510e3] #creep activation energy [J/mol]
etal = [10,1] #liquid viscsoity [Pas]
Xs_wet = np.array([19.8,19.4,17.8,14.9,12.5])
Xs_0 = np.vstack((Xs_wet+9,Xs_wet)) # initial wt % sulfur in core - first row is dry, second is wet, first column is 100km, last is 500km

#%%
# #radius [m]
nr = 5 #number of points
minr = 100e3 #min val
maxr = 500e3 #max val
r = np.array([100e3,200e3,300e3,400e3,500e3])

#%%  differentiation time [Myr]
tdiff = [2,3] #differentiation time midpoints when varying other params [Myr] Hellmann et. al. 2024
nt = 7 #number of points
mint = 1 #min val
maxt = 4 #max val
t_start_m = np.linspace(mint,maxt,nt)

#%%
#reference viscosity
neta0 = 3
mineta0 = [17,16] #10^{mineta0}
maxeta0 = [21,20]

#%%
# parameters that are always fixed
#numerical parameters
default = 'vary' #viscosity profile will change
t_end_m = 1500 #maximum end time of simulation [Myr]
w = 5 #width of linear region [K]
dr = 500 # grid space [m]
dt = 0.075 #timestep [fraction of core conductive timestep]
icfrac = 1 #fraction of solidified mass that ends up in inner core

#parameters which could be varied
rcmf = 0.5 #rheologically critical melt fraction - melting required for differentiation
alpha_n = 30 #melt weakening (median value from Sanderson et. al. 2024)
Fe0 = 1e-8 # 60Fe/56FE ratio in accreting material 
accrete = False #start from differentiation

#%%
csv_nums = {'r':1,'t_start_m':2,'eta0':3}
nruns = [nr,nt,neta0] 
variables = ['r','t_start_m','eta0']
varlab = {'r':r,'t_start_m':t_start_m,'eta0':eta0}
mod = len(variables) #how many things vary within a wet or dry body

run_nums = np.zeros([len(nruns)*2])
for i in range(len(nruns)*2): #create continuous run numbers so no duplicates across all bodies
    if i==0:
        run_nums[i] = 1
    else:
        run_nums[i] = run_nums[i-1]+nruns[(i%3)-1]
            

#%%
# create dictionaries and write to csv - change these three lines
for i, runs in enumerate(run_nums):
    if i<3:
        m = 0 #dry first
    else:
        m = 1 #wet second
    #calculated wet-dry parameters 
    rav = (minr+maxr)/2
    Tms = solidus(rav,rav/2,xwater[m],Xs_0[m,2],rhom) # mantle solidus [K] Katz et al. 2003
    beta = eac[m]/(R*Tms**2) #E/RTref^2 - Arrhenius slope
    eta0 = np.logspace(mineta0[m],maxeta0[m],neta0)
    varlab = {'r':r,'t_start_m':t_start_m,'eta0':eta0}

    run = runs
    csv_num = csv_nums[variables[i%mod]]
    var = varlab[variables[i%mod]] #choose your variable
    
    #set constant parameters - middle value of each range for wet or dry
    t_start_val = tdiff[m]
    eta0val= 10**((maxeta0[m]+mineta0[m])/2)
    rval=300e3
    Xsval = Xs_0[m,2] #midpoint sulfur content
    run_info = pd.DataFrame(columns=['run','r','rcr','default','rcmf','eta0','beta','w','etal','alpha_n','Xs_0','Fe0','t_start_m','t_end_m','dr','dt','icfrac','xwater','accrete','status']) #create columns of dataframe
    unit_row = ['','m','','','','Pas','K^-1','K','Pas','','wt %','60Fe/56Fe','Myr','Myr','m','t_cond_core','','wt %','',''] #first row is units
    run_info.loc[len(run_info)] = unit_row
    
    for val in var:
        if val == var[0]: #create csv headers
            run_info = pd.DataFrame(columns=['run','r','rcr','default','rcmf','eta0','beta','w','etal','alpha_n','Xs_0','Fe0','t_start_m','t_end_m','dr','dt','icfrac','xwater','accrete','status']) #create columns of dataframe
            unit_row = ['','m','','','','Pas','K^-1','K','Pas','','wt %','60Fe/56Fe','Myr','Myr','m','t_cond_core','','wt %','',''] #first row is units
            run_info.loc[len(run_info)] = unit_row
        #make data frame of correct length will all same values
        run_info = pd.concat([run_info, pd.DataFrame({"run":[run],"r":[rval],"rcr":[rcr[m]],"default":[default],"rcmf":[rcmf],"eta0":[eta0val],"beta":[beta],"w":[w],"etal":[etal[m]],"alpha_n":[alpha_n],"Xs_0":[Xsval], "Fe0":[Fe0], "t_start_m":[t_start_val], "t_end_m":[t_end_m], "dr":[dr],"dt":[dt],"icfrac":[icfrac],"xwater":[xwater[m]],"accrete":[False],"status":""})],ignore_index=True)
        run = run+1
    #replace variable values
    run_info.loc[1:,variables[i%mod]]=var
    if variables[i%mod] == 'r':
        run_info.loc[1:,'Xs_0']=Xs_0[m,:] #replace min sulfur content when r varies
    #create csv                      
    run_info.to_csv(f'{folder}auto_params_{i+1}.csv',index=False)                        
