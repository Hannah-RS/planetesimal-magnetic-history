#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for plotting the results of the integration to produce a figure similar to Fig 1 in Bryson 2019
"""
#import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#choose your run
run=22

#scale time to Myr
from parameters import Myr

#import data from npz file
npzfile = np.load('Results/run_{}.npz'.format(run))
Tm = npzfile['Tm_base'] 
Tm_surf = npzfile['Tm_surf'] 
Tc= npzfile['Tc'] 
f = npzfile['f'] 
T_profile = npzfile['T_profile']
t = npzfile['t'] #time in s
Rem = npzfile['Rem'] # magnetic Reynolds number from compositional and thermal convection (whatever is larger at each time step)
Flux = npzfile['Flux']
Ra = npzfile['Ra'] 
d0 = npzfile['d0'] 

Fs = Flux[0]
Fcmb = Flux[1]
Fad = Flux[2]
Frad = Flux[3]

t_plot = t/Myr

# #optional script to concatentate 2 runs
# run=23
# npzfile = np.load('Results/run_{}.npz'.format(run))
# Tm2 = npzfile['Tm_base'] 
# Tm_surf2 = npzfile['Tm_surf'] 
# Tc2= npzfile['Tc'] 
# f2 = npzfile['f'] 
# T_profile2 = npzfile['T_profile']
# t2 = npzfile['t'] #time in s
# Rem2 = npzfile['Rem'] # magnetic Reynolds number from compositional and thermal convection (whatever is larger at each time step)
# Flux2 = npzfile['Flux']
# Ra2 = npzfile['Ra'] 
# d02 = npzfile['d0'] 

# Fs2 = Flux[0]
# Fcmb2 = Flux[1]
# Fad2 = Flux[2]
# Frad2 = Flux[3]

# t_plot2 = t2/Myr

# #concatenate
# Tm = np.hstack([Tm,Tm2])
# Tc = np.hstack([Tc,Tc2])
# Flux = np.hstack([Flux,Flux2])
# Fs = Flux[0]
# Fcmb = Flux[1]
# Fad = Flux[2]
# Frad = Flux[3]
# Rem = np.hstack([Rem,Rem2])
# f = np.hstack([f,f2])
# t_plot_comb = np.hstack([t,t2])/Myr
# t_plot = t_plot_comb
# #import label info - read in from correct row in csv
# run_info = pd.read_csv('run_info3.csv',delimiter=',')

# row = run_info[run_info['run']==run]
# r = row.iloc[0,1] #radius [m]
# Tsolidus = row.iloc[0,2] #temperature for onset of solidification [K]
# Tm0 = row.iloc[0,3] # initial mantle and core temp [K]
# tstart = row.iloc[0,4]
# tend = row.iloc[0,5] # max possible time of simulation [Myr]
# tstep = row.iloc[0,6] #max timestep [Myr]
# tsolid = row.iloc[0,7] #time at which solidifcation finishes [Myr] if tsolid == tend then core may not have finished solidifying
# cond_i= int(row.iloc[0,8]) #index in array at which mantle started conducting
# cond_t = t[cond_i]/Myr #time at which mantle switched to conduction
# dr=row.iloc[0,10] #cell spacing
# dt =row.iloc[0,11] #T_profile output frequency

################### Main Plot #########################################
# make  log-log plot - similar to Bryson 2019

plt.figure(tight_layout=True)
plt.suptitle('Thermal evolution of a {:.0f}km asteroid \n Tm0 = {}K, Tsolidus ={}K \n run {}'.format(r/1e3, Tm0, Tsolidus,run))

xmin=tstart

#temperatures as function of time
plt.subplot(2,1,1)
plt.semilogx(t_plot,Tm,label='T$_m$ - base')
#plt.semilogx(t_plot,Tm_surf,label='T$_m$ - surface')
plt.semilogx(t_plot,Tc,label='T$_c$')
#plt.vlines(cond_t,ymin=min(Tm),ymax=1600,color='black',linestyle='--',label='conduction')
#plt.xlim([xmin,max(t_plot)])
#plt.xlabel('Time/ Myr')
plt.ylabel('T/K')
plt.legend(loc='lower right')

#fluxes as function of time
plt.subplot(2,1,2)
plt.loglog(t_plot,Fs,label='$F_s$')
plt.loglog(t_plot,Fcmb,label='$F_{CMB}$')
plt.loglog(t_plot,Fad,label='$F_{ad}$')
plt.loglog(t_plot,Frad,label='$F_{rad}$')
#plt.xlabel('Time/ Myr')
#plt.xlim([xmin,max(t_plot)])
plt.ylim([1e-3,1e2])
plt.ylabel('Flux/ W$m^{-2}$')
plt.legend()
plt.savefig('Plots/Tflux_run{}.png'.format(run),dpi=450)

plt.figure(tight_layout=True)
plt.suptitle('Thermal evolution of a {:.0f}km asteroid \n Tm0 = {}K, Tsolidus ={}K \n run {}'.format(r/1e3, Tm0, Tsolidus,run))
plt.subplot(2,1,1)
plt.loglog(t_plot,Rem)
#plt.xlim([xmin,max(t_plot)])
plt.hlines(10,xmin=0,xmax=t_plot[len(Rem)-1],color='k',linestyle='--')
#plt.xlabel('Time/Myr')
plt.ylabel('Rem')
plt.ylim([1,100])

plt.subplot(2,1,2)
plt.semilogx(t_plot,f,label='f')
#plt.xlim([xmin,max(t_plot)])
plt.xlabel('Time/ Myr')
plt.ylabel('f')

plt.savefig('Plots/Remf_run{}.png'.format(run),dpi=450)

####################### Exploratory Plots #################################
#make a plot to check Rayleigh number

plt.figure()
plt.loglog(t_plot[:cond_i],Ra[:cond_i])
plt.hlines(1000,min(t_plot),max(t_plot[:cond_i]),color='k',linestyle='--')
plt.xlim([xmin,max(t_plot)])
plt.xlabel('Time/Myr')
plt.ylabel('Rayleigh number')
#plt.savefig('Plots/Rayleigh_run_{}.png'.format(run))


#lid thickness
plt.figure()
plt.plot(t_plot[:cond_i],d0[:cond_i]/1e3)
#plt.hlines(1,min(t_plot),max(t_plot[:cond_i]),color='k',linestyle='--')
plt.xlim([xmin,t_plot[cond_i]])
plt.xlabel('Time/Myr')
plt.ylabel('Lid thickness/ km')
#plt.savefig('Plots/lid_thickness_run_{}.png'.format(run))

plt.figure()
plt.plot(t_plot[:cond_i],2*d0[:cond_i]/r)
#plt.hlines(1,min(t_plot),max(t_plot[:cond_i]),color='k',linestyle='--')
plt.xlim([xmin,t_plot[cond_i]])
plt.xlabel('Time/Myr')
plt.ylabel('Lid thickness/ mantle thickness')
#plt.savefig('Plots/lid_thickness_run_{}.png'.format(run))

#plot temperature profiles 
plt.figure()
rplot= np.arange(0,r,dr)
n = np.shape(T_profile)[0]
l = len(t_plot)
n_plot = 4 #how many plots do you want
for i in range(n_plot):
    plt.plot(rplot/1e3, T_profile[i*int(n/n_plot),:],label='{:.0f} Myr'.format(t_plot[i*int(l/n_plot)])) #approximate temp profile times
plt.xlabel('Distance from centre of asteroid /km')
plt.ylabel('Temperature / K')
plt.xlim([80,100])
plt.title('{:.0f}km asteroid \n Tm0 = {}K, Tsolidus ={}K, run {}'.format(r/1e3, Tm0, Tsolidus, run))
plt.legend()
#plt.savefig('Plots/Tprofile_run{}.png'.format(run))

#initial temperature profiles
plt.figure()
rplot= np.arange(0,r,dr)
for i in range(10):
    plt.plot(rplot/1e3, T_profile[i,:],label='{:.1f} Myr'.format(10*i)) #approximate temp profile times
plt.xlabel('Distance from centre of asteroid /km')
plt.ylabel('Temperature / K')
plt.title('{:.0f}km asteroid \n Tm0 = {}K, Tsolidus ={}K, run {}'.format(r/1e3, Tm0, Tsolidus, run))
plt.legend()

#look at temp evolution at different depths
# is it zig zaggy at all depths?
t_plot2 = np.arange(tstart+10,tsolid,dt/Myr)
plt.figure()
for i in range(0,5):
    plt.plot(t_plot2,T_profile[:,2*(-i-1)],label='{} km'.format(i))
plt.xlabel('time /Myr')
plt.ylabel('T/K')
#plt.xlim([2.5,10])
plt.legend(loc='lower right')
#plt.savefig('Plots/profile_dt_comp.png')

#when is T[-2] zig-zaggy?
plt.figure(tight_layout=True,figsize=[6,6])
plt.title('Thermal evolution of a {:.0f}km asteroid \n Tm0 = {}K, Tsolidus ={}K \n run {}'.format(r/1e3, Tm0, Tsolidus,run))
xmin=tstart
#plt.plot(t_plot,Tm,label='T$_m$ - base')
plt.scatter(t_plot,Tm_surf,label='T$_m$ - surface')
#plt.plot(t_plot,Tc,label='T$_c$')
plt.vlines(cond_t,ymin=min(Tm),ymax=1600,color='black',linestyle='--',label='conduction')
plt.xlim([0,5])
#plt.xlabel('Time/ Myr')
plt.ylabel('T/K')
plt.legend(loc='lower right')
#plt.savefig('Plots/scatter_16.png')