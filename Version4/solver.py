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
from parameters import  Myr, Tm0, Tc0, Ts, f0, r, rc, dr, kappa_c, out_interval, km, cpm_p, rhom, save_interval, default, kappa


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
run =66

t_acc=1*Myr  #Accretion time
t_end_m=200#end time in Myr
t_end=t_end_m*Myr
t_cond_core = dr**2/kappa_c #conductive timestep for core
t_cond_mantle = dr**2/kappa #conductive timestep for mantle
step_m=0.1*t_cond_core  #max timestep must be smaller than conductive timestep
n_save = int(save_interval/step_m)

# set initial temperature profile
n_cells = int(r/dr) #number of cells needed to span body
Tint = np.ones([n_cells])*1600 #first element in the array is at r=0, accrete cold at surface temp 
Tint[-1]=Ts
print('Initial conditions set')

###############  sintering and compaction code will go here  ##################

########################### Differentiation ###################################
# tic = time.perf_counter()
# Tdiff, Tdiff_profile, k_profile, Xfe, rho_profile, Ra, Ra_crit, convect, t_diff = differentiation(Tint,t_acc,r, dr, step_m)
# toc = time.perf_counter()
# diff_time = toc - tic  

# # update user on progress and plot differentiated temperature profile 
# rplot= np.arange(0,r,dr)/1e3

# plt.figure()
# plt.plot(rplot,Tdiff_profile)
# plt.xlabel('r/km')
# plt.ylabel('Temperature/K')
# plt.title('Temperature profile post differentiation')

# print('Differentiation complete. It took', time.strftime("%Hh%Mm%Ss", time.gmtime(diff_time)))
# np.savez(f'Results/diff_run_{run}', Tdiff = Tdiff, Tdiff_profile = Tdiff_profile, k_profile = k_profile, Xfe = Xfe, rho_profile = rho_profile, Ra = Ra, Ra_crit = Ra_crit, convect = convect, t_diff = t_diff)
# raise ValueError('Differentiation passed')
######################## Thermal evolution ####################################

#integrate
tic = time.perf_counter()
Tc, Tc_conv, Tcmb, Tm_mid, Tm_conv, Tm_surf, Tprofile, f, Xs, dl, dc, d0, Ra, Fs, Fad, Fcmb, t, cond_i = thermal_evolution(t_acc,t_end,step_m,Tint,f0,sparse_mat_c,sparse_mat_m) 
toc = time.perf_counter()
int_time = toc - tic    
print('Thermal evolution complete', time.strftime("%Hh%Mm%Ss", time.gmtime(int_time)))

############################# Process data ####################################

#Reduce data points - as model saves more often than needed
# take every nth point at an interval specified by save_interval in parameters.py
Tc= Tc[0::n_save]
Tc_conv = Tc_conv[0::n_save]
Tcmb = Tcmb[0::n_save]
Tm_mid = Tm_mid[0::n_save]
Tm_conv = Tm_conv[0::n_save]
Tm_surf = Tm_surf[0::n_save] 
f = f[0::n_save]
Xs = Xs[0::n_save]
dl = dl[0::n_save]
dc = dc[0::n_save]
d0 = d0[0::n_save]
Ra = Ra[0::n_save]
Fs = Fs[0::n_save]
Fad = Fad[0::n_save]
Fcmb = Fcmb[0::n_save]
t = t[0::n_save] 
if cond_i != 'nan':
    cond_i = int(cond_i/n_save) #scale cond_i too

# calculate Frad, Rem
from parameters import km, kc, G, rhoc, alpha_c, cpc, gamma


# Frad - radiogenic heat flux, normalised to surface of body
from parameters import h0, Al0, XAl, thalf_al, rhom, Vm, As
h = h0*Al0*XAl*np.exp(-np.log(2)*t/thalf_al) 
Frad = h*rhom*Vm/As #radiogenic heatflux

#combine these in a single array
Flux = [Fs, Fcmb, Fad, Frad]

# calculate both magnetic reynolds numbers and merge only keeping the larger value (the larger velocity will dominate)
from Rem_calc import Rem_comp, Rem_therm, ucomp_nimmo, ucomp_aubert, p_nichols

# need Fdrive for Rem_therm
#only calculate this for Fdrive >0
Fdrive = Fcmb - Fad
Fdrive_nn = Fdrive.copy()
Fdrive_nn[Fdrive<0]=0
Rem_t = Rem_therm(Fdrive_nn) # magnetic Reynolds number for thermal convection

#calculate compositional convection two ways
unimmo = ucomp_nimmo(t,f)
power = p_nichols(t,f) #convective power density
uaubert = ucomp_aubert(power,f)
Rem_cn = Rem_comp(unimmo,f) # magnetic Reynolds number for compositional convection
Rem_ca = Rem_comp(uaubert,f) # magnetic Reynolds number for compositional convection
Rem1 = Rem_cn
Rem2 = Rem_ca
Rem1[Rem_cn<Rem_t] = Rem_t[Rem_cn<Rem_t] #replace values where Rem_t < Rem_c
Rem2[Rem_ca<Rem_t] = Rem_t[Rem_ca<Rem_t] #replace values where Rem_t < Rem_c

print('Fluxes and magnetic Reynolds number calculated.')

############################ Save results #####################################
# save variables to file
np.savez('Results/run_{}'.format(run), Tc = Tc, Tc_conv = Tc_conv, Tcmb = Tcmb,  Tm_mid = Tm_mid, Tm_conv = Tm_conv, Tm_surf = Tm_surf, T_profile = Tprofile, f=f, Xs = Xs, dl = dl, dc=dc, d0 = d0, Ra = Ra, t=t, Rem1 = Rem1, Rem2 = Rem2, Flux = Flux) 

#write parameters to the run file
from csv import writer
from parameters import r, Tm0, default

var_list = [run, r, Tm0, t_acc/Myr, t_end_m, step_m/Myr, max(t)/Myr, cond_i, int_time, dr, out_interval/Myr, default]

    
with open('run_info4.csv','a') as f_object:
    writer_object = writer(f_object) #pass file object to csv.writer
    writer_object.writerow(var_list) # pass list as argument into write row
    f_object.close() #close file

print('Results and run parameters saved')
