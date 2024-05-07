#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test conductive stencil
"""
from temp_cond import Tm_cond_calc
from stencil import cond_stencil_mantle
import scipy.sparse as sp
import numpy as np
from parameters import r, dr, Myr, step_m, km, rhom, rc
dT_mat_m = cond_stencil_mantle(r,rc,dr,km/rhom)  
sparse_mat_m = sp.dia_matrix(dT_mat_m)


n_cells = round((int(r/dr)+1)/2) #number of cells needed to span body - add one to include the centre
T_old_mantle = np.ones([n_cells+1])*1600
T_old_mantle[-1] = 200
nsave = 138
Tout = np.zeros([nsave,n_cells+1])
tout = np.zeros([nsave])
tsolve_new = 10*Myr
i = 0
dt = step_m

while tsolve_new < 110*Myr:
    i = i + 1
    tsolve_new = tsolve_new + dt
    T_new_mantle = Tm_cond_calc(tsolve_new,dt,T_old_mantle,sparse_mat_m)
    T_old_mantle = T_new_mantle
    if i%10000 ==0:
        Tout[int(i/10000)-1,:] = T_old_mantle
        tout[int(i/10000)-1] = tsolve_new/Myr
#%%
## make plot
import matplotlib.pyplot as plt
rplot = np.arange(rc,int(r)+dr,int(dr))/1e3
plt.figure(rasterized=True,figsize=[15,10])
plt.pcolormesh(tout,rplot,np.transpose(Tout),shading = 'gouraud')
plt.ylabel('r /km')
plt.xlabel('t / Myr')
plt.xscale('log')
#plt.ylim([0,r/1e3])
plt.colorbar(label='T/K')


#%%
plt.figure()
plt.plot(tout,Tout[:,0])
plt.plot(tout,Tout[:,1])
plt.ylabel('CMB temperature')
plt.xlabel('Time /Myr')