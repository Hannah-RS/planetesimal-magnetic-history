#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for solving the thermal evolution of an asteroid. Doesn't include differentiation (i.e commented out)'
"""
# import modules
import numpy as np
import matplotlib.pyplot as plt
import time #use this to time the integration

#import time constants and initial conditions
from parameters import  Myr, Ts, f0, r, rc, dr, kappa_c, save_interval_d, save_interval_t, km, Vm, As, cpm_p, rhom, default, kappa, rcmf, Xs_0, Fe0


#calculate the stencil for the conductive profile, save so can be reloaded in later steps
from stencil import cond_stencil_core, cond_stencil_mantle
from diff_function import differentiation
from full_euler import thermal_evolution

import scipy.sparse as sp
kappa_m = km/(rhom*cpm_p) #use modified specific heat capacity to account for mantle melting
dT_mat_c = cond_stencil_core(r,rc,dr,kappa_c) 
dT_mat_m = cond_stencil_mantle(r,rc,dr,kappa_m)  
sparse_mat_m = sp.dia_matrix(dT_mat_m)
sparse_mat_c = sp.dia_matrix(dT_mat_c)


# define the run number, start and end times
run =12

t_acc=0.8*Myr  #Accretion time
t_end_m=30#end time in Myr
t_end=t_end_m*Myr
t_cond_core = dr**2/kappa_c #conductive timestep for core
t_cond_mantle = dr**2/kappa #conductive timestep for mantle
step_m=0.1*t_cond_core  #max timestep must be smaller than conductive timestep
n_save_d = int(save_interval_d/step_m)
n_save_t = int(save_interval_t/step_m)

# set initial temperature profile
n_cells = int(r/dr) #number of cells needed to span body
Tint = np.ones([n_cells])*Ts#first element in the array is at r=0, accrete cold at surface temp 
Tint[-1]=Ts

print('Initial conditions set')

###############  sintering and compaction code will go here  ##################

########################### Differentiation ###################################
tic = time.perf_counter()
Tdiff, Xfe, Xsi, cp, Ra, Ra_crit, convect, d0, t_diff, H  = differentiation(Tint,t_acc,r, dr, step_m)
toc = time.perf_counter()
diff_time = toc - tic  

# update user on progress and plot differentiated temperature profile 
rplot= np.arange(0,r,dr)/1e3

plt.figure()
plt.scatter(rplot,Tdiff[:,-1])
plt.xlabel('r/km')
plt.ylabel('Temperature/K')
plt.title('Temperature profile post differentiation')

#rescale data and save here in case thermal evolution crashes
Tdiff = Tdiff[:,0::n_save_d]
Xfe = Xfe[:,0::n_save_d]
Xsi = Xsi[:,0::n_save_d]
cp = cp[:,0::n_save_d]
Ra = Ra[0::n_save_d]
Ra_crit = Ra_crit[0::n_save_d]
convect = convect[0::n_save_d]
d0 = d0[0::n_save_d]
t_diff = t_diff[0::n_save_d]
H = H[0::n_save_d]

np.savez(f'Results_combined/run_{run}_diff', Tdiff = Tdiff, Xfe = Xfe, Xsi = Xsi, cp = cp, Ra = Ra, Ra_crit = Ra_crit, convect = convect, d0=d0, t_diff = t_diff, H=H)

print('Differentiation complete. It took', time.strftime("%Hh%Mm%Ss", time.gmtime(diff_time)))
######################## Thermal evolution ####################################

#integrate
tic = time.perf_counter()
Tc, Tc_conv, Tcmb, Tm_mid, Tm_conv, Tm_surf, Tprofile, f, Xs, dl, dc, d0, min_unstable, Ra, RaH, RanoH, Racrit, Fs, Flid, Fad, Fcmb, Rem_c, t, cond_t = thermal_evolution(t_diff[-1],t_end,step_m,Tdiff[:,-1],f0,sparse_mat_c,sparse_mat_m) 
toc = time.perf_counter()
int_time = toc - tic    

#update on progress
plt.figure()
plt.scatter(rplot,Tprofile[-1,:])
plt.xlabel('r/km')
plt.ylabel('Temperature/K')
plt.title('Temperature profile post thermal evolution')

print('Thermal evolution complete', time.strftime("%Hh%Mm%Ss", time.gmtime(int_time)))

############################# Process data ####################################
########## Current comparitive parameters ####################

from parameters import convect_ratio
nmantle = int((r/dr)/2)
diff_time = t_diff[-1]/Myr
diff_T = Tdiff[int(nmantle),-1]
peakT = np.amax(Tprofile[:,nmantle:])
loc_max = np.where(Tprofile[:,nmantle:]==peakT)[1][0] #take the set of time coordinates and first value (they should all be the same)
tmax = t[loc_max]/Myr
if np.all(Tprofile[:,int(nmantle)-2]<Tcmb):
    tstrat_remove = np.inf
else:
    tstrat_remove = t[Tcmb < Tprofile[:,int(nmantle)-2]][0]/Myr
    
if np.any(min_unstable==0):
     strat_end = t[np.where(min_unstable==0)[0]][0]/Myr
else:
    strat_end = np.inf
    
if np.all(Fcmb < Fad):
    super_ad_start = np.inf
    super_ad_end = np.inf
else:
    super_ad_start = t[np.where(Fcmb>Fad)[0]][0]/Myr
    super_ad_end = t[np.where(Fcmb>Fad)[0]][-1]/Myr


# Frad - radiogenic heat flux, normalised to surface of body
from heating import Al_heating
h = Al_heating(t) 
Frad = h*rhom*Vm/As #radiogenic heatflux

#combine these in a single array
Flux = [Fs, Fcmb, Fad, Frad]

# calculate thermal magnetic reynolds number
from Rem_calc import Rem_therm
from parameters import Xs_eutectic
# need Fdrive for Rem_therm
#only calculate this for Fdrive >0
Fdrive = Fcmb - Fad
Fdrive_nn = Fdrive.copy()
Fdrive_nn[Fdrive<0]=0
Fdrive_nn[Xs>=Xs_eutectic]=0 #no dynamo in eutectic solidification
Rem_t = Rem_therm(Fdrive_nn,f,min_unstable) # magnetic Reynolds number for thermal convection - tuple of MAC and CIA balance

print('Fluxes and magnetic Reynolds number calculated.')

############################ Save results #####################################
# save variables to file
np.savez('Results_combined/run_{}'.format(run), Tc = Tc, Tc_conv = Tc_conv, Tcmb = Tcmb,  Tm_mid = Tm_mid, Tm_conv = Tm_conv, Tm_surf = Tm_surf, T_profile = Tprofile, Flid = Flid, f=f, Xs = Xs, dl = dl, dc=dc, d0 = d0, min_unstable=min_unstable, Ra = Ra, RaH= RaH, RanoH = RanoH, Racrit = Racrit, t=t, Rem_t = Rem_t, Rem_c = Rem_c, Flux = Flux) 

#write parameters to the run file
from csv import writer

var_list = [run, r, dr, t_acc/Myr, t_end_m, step_m/Myr, max(t)/Myr, cond_t, int_time,  save_interval_d/Myr, save_interval_t/Myr, default, rcmf, Xs_0, Fe0]

    
with open('run_info.csv','a') as f_object:         
    writer_object = writer(f_object) #pass file object to csv.writer
    writer_object.writerow(var_list) # pass list as argument into write row
    f_object.close() #close file
  
# #comparative parameters
# var_list2 = [convect_ratio, diff_time, peakT, tmax, tstrat_remove, strat_end, super_ad_start, super_ad_end, diff_T, run, r, rcmf]

# with open('lid_test.csv','a') as f_object:
#     writer_object = writer(f_object) #pass file object to csv.writer
#     writer_object.writerow(var_list2) # pass list as argument into write row
#     f_object.close() #close file

# print('Results and run parameters saved')
