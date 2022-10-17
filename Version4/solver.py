#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for solving the thermal evolution of an asteroid
"""
# import modules
import numpy as np
import matplotlib.pyplot as plt
import time #use this to time the integration

#import time constants and initial conditions
from parameters import  Myr, Tm0, Tc0, Ts, f0, r, rc, dr, kappa, kappa_c, out_interval


#calculate the stencil for the conductive profile, save so can be reloaded in later steps
from stencil import cond_stencil_core, cond_stencil_mantle
import scipy.sparse as sp
dT_mat_c = cond_stencil_core(r,rc,dr,kappa_c) 
dT_mat_m = cond_stencil_mantle(r,rc,dr,kappa)  
sparse_mat_m = sp.dia_matrix(dT_mat_m)
sparse_mat_c = sp.dia_matrix(dT_mat_c)


# define the run number, start and end times
run = 40

t_start=1*Myr #start after the end of stage 2
t_end_m=1000 #end time in Myr
t_end=t_end_m*Myr
t_cond = dr**2/kappa #conductive timestep
step_m=0.1*t_cond  #max timestep must be smaller than conductive timestep

# set initial temperature profile
n_cells = int(r/dr) #number of cells needed to span body
Tint = np.ones([n_cells]) #first element in the array is at r=0 

#set core temperature
n_core = int(rc/dr)
Tint[:n_core] = Tc0 #core initially isothermal

#set mantle temperature - assume isothermal below stagnant lid
Tint[n_core:] =Tm0

#find initial stagnant lid
from Rayleigh_def import Rayleigh_calc
Ram, d0 = Rayleigh_calc(Tm0)
nlid_cells = int(d0/dr)
for i in range(nlid_cells):
    Tint[-i-1]=Ts+(Tm0-Ts)*i*dr/d0 #surface at Ts, initial linear profile in small lid

# update user on progress and plot initial temperature profile so can check
rplot= np.arange(0,r,dr)/1e3

plt.figure()
plt.plot(rplot,Tint)
plt.xlabel('r/km')
plt.ylabel('Temperature/K')
plt.title('Initial temperature profile')

print('Initial conditions set')

# set solver running  
from full_euler import thermal_evolution

#integrate
tic = time.perf_counter()
Tc, Tcmb, Tm_mid, Tm_conv, Tm_surf, Tprofile, f, Xs, bl, d0, Ra, Fs, Fad, Fcmb, t, cond_i = thermal_evolution(t_start,t_end,step_m,Tint,f0,sparse_mat_c,sparse_mat_m) 
toc = time.perf_counter()
int_time = toc - tic    
print('Integration finished')

# calculate Fs, Fcmb, Fad, Frad, Rem


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

#save variables to file
np.savez('Results/run_{}'.format(run), Tc = Tc, Tcmb = Tcmb, Tm_mid = Tm_mid, Tm_conv = Tm_conv, Tm_surf = Tm_surf, T_profile = Tprofile, f=f, Xs = Xs, bl = bl, d0 = d0, Ra = Ra, t=t, Rem1 = Rem1, Rem2 = Rem2, Flux = Flux) 

#write parameters to the run file
from csv import writer
from parameters import r, Tsolidus, Tm0

var_list = [run, r, Tsolidus, Tm0, t_start/Myr, t_end_m, step_m/Myr, max(t)/Myr, cond_i, int_time, dr, out_interval]

    
with open('run_info3.csv','a') as f_object:
    writer_object = writer(f_object) #pass file object to csv.writer
    writer_object.writerow(var_list) # pass list as argument into write row
    f_object.close() #close file

print('Results and run parameters saved')
