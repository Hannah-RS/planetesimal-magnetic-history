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
from stencil import cond_stencil
dT_mat = cond_stencil(r,rc,dr,kappa,kappa_c)  
np.savez('Stencil',dT_mat=dT_mat)

# define the run number, start and end times
run = 22

t_start=1*Myr #start after the end of stage 2
t_end_m=400 #end time in Myr
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
Ra, d0, Tprofile, Tc, Tm_base, Tm_surf, f, t, cond_i = thermal_evolution(t_start,t_end,step_m,Tint,f0) #Tm_base is temp at base of mantle, Tmsurf is temp one node below surface
toc = time.perf_counter()
int_time = toc - tic    
print('Integration finished')

# calculate Fs, Fcmb, Fad, Frad, Rem


from parameters import km, kc, G, rhoc, alpha_c, cpc, gamma

#calculate Fs
Fs = km*(Tm_base - Ts)/d0

#Fcmb
Fcmb = km*(Tc - Tm_base)/dr #expression for CMB flux in conductive regime
Fcmb[:cond_i] = Fs[:cond_i]*(Tc[:cond_i] - Tm_base[:cond_i])*gamma #expression for CMB flux in convective regime

# Fad
gc = G*4/3*np.pi*rc*rhoc #g at CMB
Fad = kc*Tc*alpha_c*gc/cpc

# Frad - radiogenic heat flux, normalised to surface of body
from parameters import h0, Al0, XAl, thalf_al, rhom, Vm, As
h = h0*Al0*XAl*np.exp(-np.log(2)*t/thalf_al) 
Frad = h*rhom*Vm/As #radiogenic heatflux

#combine these in a single array
Flux = [Fs, Fcmb, Fad, Frad]

# calculate both magnetic reynolds numbers and merge only keeping the larger value (the larger velocity will dominate)
from Rem_calc import Rem_comp, Rem_therm

# need Fdrive for Rem_therm
#only calculate this for Fdrive >0
Fdrive = Fcmb - Fad
Fdrive_nn = Fdrive.copy()
Fdrive_nn[Fdrive<0]=0
Rem_t = Rem_therm(Fdrive_nn) # magnetic Reynolds number for thermal convection
Rem_c = Rem_comp(t,f) # magnetic Reynolds number for compositional convection
Rem = Rem_c
Rem[Rem_c<Rem_t] = Rem_t[Rem_c<Rem_t] #replace values where Rem_t < Rem_c

print('Fluxes and magnetic Reynolds number calculated.')

#save variables to file
np.savez('Results/run_{}'.format(run), Tm_base = Tm_base, Tm_surf = Tm_surf, Tc = Tc, f = f, T_profile = Tprofile, t=t, Rem = Rem, Flux = Flux, Ra = Ra, d0 = d0 ) 

#write parameters to the run file
from csv import writer
from parameters import r, Tsolidus, Tm0

var_list = [run, r, Tsolidus, Tm0, t_start/Myr, t_end_m, step_m/Myr, max(t)/Myr, cond_i, int_time, dr, out_interval]

    
with open('run_info3.csv','a') as f_object:
    writer_object = writer(f_object) #pass file object to csv.writer
    writer_object.writerow(var_list) # pass list as argument into write row
    f_object.close() #close file

print('Results and run parameters saved')
