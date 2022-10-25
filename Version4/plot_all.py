#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for plotting all the results of the integration to check what happened
"""
#import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#choose your run
run=46
conduction = True # does the mantle switch to conduction?
save = False # do you want to save your figures?
#scale time to Myr
from parameters import Myr, r, Tm0, Tsolidus

#import data from npz file
npzfile = np.load('Results/run_{}.npz'.format(run))
Tc= npzfile['Tc'] 
Tc_conv = npzfile['Tc_conv']
Tcmb = npzfile['Tcmb']
Tm_mid = npzfile['Tm_mid']
Tm_conv = npzfile['Tm_conv']
Tm_surf = npzfile['Tm_surf'] 
T_profile = npzfile['T_profile']
f = npzfile['f'] 

t = npzfile['t'] #time in s
Rem1 = npzfile['Rem1'] # magnetic Reynolds number from compositional (Nimmo) and thermal convection (whatever is larger at each time step)
Rem2 = npzfile['Rem2'] # magnetic Reynolds number from compositional (Nichols) and thermal convection (whatever is larger at each time step)
Flux = npzfile['Flux']
Ra = npzfile['Ra'] 
d0 = npzfile['d0'] 
Xs = npzfile['Xs']
bl = npzfile['bl']
Fs = Flux[0]
Fcmb = Flux[1]
Fad = Flux[2]
Frad = Flux[3]

t_plot = t/Myr


#import label info - read in from correct row in csv
run_info = pd.read_csv('run_info3.csv',delimiter=',')

row = run_info[run_info['run']==run]
r = row.iloc[0,1] #radius [m]
Tsolidus = row.iloc[0,2] #temperature for onset of solidification [K]
Tm0 = row.iloc[0,3] # initial mantle and core temp [K]
tstart = row.iloc[0,4]
tend = row.iloc[0,5] # max possible time of simulation [Myr]
tstep = row.iloc[0,6] #max timestep [Myr]
tsolid = row.iloc[0,7] #time at which solidifcation finishes [Myr] if tsolid == tend then core may not have finished solidifying
cond_i= row.iloc[0,8] #index in array at which mantle started conducting

if conduction == True:
    cond_i = int(cond_i)
    cond_t = t[cond_i]/Myr #time at which mantle switched to conduction
    
dr=row.iloc[0,10] #cell spacing
dt =row.iloc[0,11] #T_profile output frequency

##################### Temperature and Flux Plot ###############################
# make  log-log plot - similar to Bryson 2019

plt.figure(tight_layout=True,figsize=[10,5])
plt.suptitle('Thermal evolution of a {:.0f}km asteroid \n Tm0 = {}K, Tsolidus ={}K \n run {} new viscosity'.format(r/1e3, Tm0, Tsolidus,run))

xmin=tstart

#temperatures as function of time
plt.subplot(2,1,1)
plt.semilogx(t_plot,Tc,label='T$_c$')
plt.semilogx(t_plot[Tc_conv!=0],Tc_conv[Tc_conv!=0],label='convective T$_c$')
plt.semilogx(t_plot,Tcmb,label='T$_{cmb}$')
plt.semilogx(t_plot,Tm_mid,label='T$_m$')
plt.semilogx(t_plot[Tm_conv!=0],Tm_conv[Tm_conv!=0],label='convective T$_m$')
plt.semilogx(t_plot,Tm_surf,label='T$_m$ - surface')
if conduction == True:
    plt.vlines(cond_t,ymin=min(Tm_surf),ymax=1600,color='black',linestyle='--',label='conduction')
plt.xlim([xmin,max(t_plot)])
#plt.ylim([1400,1650]) #use these limits when comparing runs
plt.ylabel('T/K')
plt.legend(loc='lower left', ncol= 2)

#fluxes as function of time
plt.subplot(2,1,2)
plt.loglog(t_plot,Fs,label='$F_s$')
plt.loglog(t_plot,Fcmb,label='$F_{CMB}$')
plt.loglog(t_plot,Fad,label='$F_{ad}$')
plt.loglog(t_plot,Frad,label='$F_{rad}$')
plt.xlabel('Time/ Myr')

plt.ylim([1e-3,1e2])   #use these limits when comparing runs
plt.xlim([xmin,max(t_plot)])
plt.ylabel('Flux/ W$m^{-2}$')
plt.legend(loc='upper right',ncol=2)


if save == True:
    plt.savefig('Plots/run_{}_Tflux.png'.format(run),dpi=450)

################### Temperature gradient - plot in progress ###################
Tgrad = np.gradient(Tc_conv)
plt.figure()
#plt.semilogy(t_plot[Tgrad>0],Tgrad[Tgrad>0])
plt.semilogy(t_plot,Tgrad)
#plt.xlim([200,300])
#plt.ylim([-0.1,0.1])

################### Magnetic Reynolds number and core size plots ##############
plt.figure(tight_layout=True)
plt.suptitle('Thermal evolution of a {:.0f}km asteroid \n Tm0 = {}K, Tsolidus ={}K \n run {} new viscosity'.format(r/1e3, Tm0, Tsolidus,run))
plt.subplot(2,1,1)
plt.loglog(t_plot,Rem1,label='Nimmo')
plt.loglog(t_plot,Rem2,label='Nichols')
plt.xlim([200,300])
plt.hlines(10,xmin=0,xmax=t_plot[len(Rem1)-1],color='k',linestyle='--')
#plt.xlabel('Time/Myr')
plt.ylabel('Rem')
plt.legend(loc='upper left')
plt.ylim([1,100])

plt.subplot(2,1,2)
plt.semilogx(t_plot,f,label='f')
plt.xlim([200,300])
plt.xlabel('Time/ Myr')
plt.ylabel('f')

if save == True:
    plt.savefig('Plots/run_{}_Remf.png'.format(run),dpi=450)

########################### Rayleigh number ###################################

plt.figure()
if conduction == True:
    plt.loglog(t_plot[:cond_i],Ra[:cond_i])
    plt.hlines(1000,min(t_plot),max(t_plot[:cond_i]),color='k',linestyle='--',label='Critical Rayleigh number')
else:
    plt.loglog(t_plot,Ra)
    plt.hlines(1000,min(t_plot),max(t_plot),color='k',linestyle='--',label='Critical Rayleigh number')
plt.xlabel('Time/Myr')
plt.ylabel('Rayleigh number')
plt.legend()

if save == True:
    plt.savefig('Plots/run_{}_Rayleigh.png'.format(run))


######################## Stagnant lid thickness ###############################

plt.figure()

if conduction == True:
    plt.plot(t_plot[:cond_i],2*d0[:cond_i]/r)
    plt.hlines(1,min(t_plot),max(t_plot[:cond_i]),color='k',linestyle='--')
    plt.xlim([xmin,t_plot[cond_i]])
else:
    plt.plot(t_plot,2*d0/r)

plt.xlabel('Time/Myr')
plt.ylabel('Stagnant Lid thickness/ mantle thickness')

if save == True:
    plt.savefig('Plots/run_{}_lid_thickness.png'.format(run))

######################## CMB mantle boundary layer thickness ##################
plt.figure()

if conduction == True:
    plt.plot(t_plot[:cond_i],2*bl[:cond_i]/r)
    #plt.hlines(1,min(t_plot),max(t_plot[:cond_i]),color='k',linestyle='--')
    plt.xlim([xmin,t_plot[cond_i]])
else:
    plt.plot(t_plot,2*bl/r)
plt.xlabel('Time/Myr')
plt.ylabel('CMB boundary layer thickness/ mantle thickness')
if save == True:
    plt.savefig('Plots/run_{}_cmb_thickness.png'.format(run))

################## Sulfur fraction in liquid part of the core #################
plt.figure()
#plt.loglog(t_plot,Xs)
plt.scatter(t_plot,Xs)
#plt.xlim([297.5,299.5])
#plt.xlim([xmin,max(t_plot)])
plt.xlabel('Time/Myr')
plt.ylabel('Liquid core sulfur content/ wt%')

if save == True:
    plt.savefig('Plots/run_{}_sulfur.png'.format(run))

###################  Temperature profiles first 100 Myr #######################
plt.figure()
rplot= np.arange(0,r,dr)
for i in range(10):
    plt.plot(rplot/1e3, T_profile[i,:],label='{:.1f} Myr'.format(10*i)) #approximate temp profile times
plt.xlabel('Distance from centre of asteroid /km')
plt.ylabel('Temperature / K')
plt.title('{:.0f}km asteroid \n Tm0 = {}K, Tsolidus ={}K, run {}'.format(r/1e3, Tm0, Tsolidus, run))
plt.legend()
if save == True:
    plt.savefig('Plots/run_{}_initial_temp.png'.format(run))

################# Temperature profiles across whole simulation ################

plt.figure()
rplot= np.arange(0,r,dr)
n = np.shape(T_profile)[0]
l = len(t_plot)
n_plot = 5 #how many plots do you want
for i in range(n_plot):
    plt.plot(rplot/1e3, T_profile[i*int(n/n_plot),:],label='{:.0f} Myr'.format(t_plot[i*int(l/n_plot)])) #approximate temp profile times
plt.xlabel('Distance from centre of asteroid /km')
plt.ylabel('Temperature / K')
plt.title('{:.0f}km asteroid \n Tm0 = {}K, Tsolidus ={}K, run {}'.format(r/1e3, Tm0, Tsolidus, run))
plt.legend()
if save == True:
    plt.savefig('Plots/run{}_Tprofile.png'.format(run))
