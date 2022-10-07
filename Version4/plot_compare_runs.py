#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare outputs of two models by plotting their fluxes and temperatures with time overlaid,
as well as magnetic Reynolds number and solid core fraction
"""
#import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#choose your runs
run1 = 38
run2 = 39

#scale time to Myr
from parameters import Myr, r, Tm0, Tsolidus

#import data from npz file - run1
npzfile = np.load('Results/run_{}.npz'.format(run1))
Tm1 = npzfile['Tm_base'] 
Tm_surf1 = npzfile['Tm_surf'] 
Tc1= npzfile['Tc'] 
t1 = npzfile['t'] #time in s
Flux1 = npzfile['Flux']
f1 = npzfile['f']
Rem1 = npzfile['Rem1'] # magnetic Reynolds number from compositional (Nimmo) and thermal convection (whatever is larger at each time step)

Fs1 = Flux1[0]
Fcmb1 = Flux1[1]
Fad1 = Flux1[2]
Frad1 = Flux1[3]

t_plot1 = t1/Myr

#import data from npz file - run2
npzfile = np.load('Results/run_{}.npz'.format(run2))
Tm2 = npzfile['Tm_base'] 
Tm_surf2 = npzfile['Tm_surf'] 
Tc2= npzfile['Tc'] 
t2 = npzfile['t'] #time in s
Flux2 = npzfile['Flux']
f2 = npzfile['f']
Rem2 = npzfile['Rem1'] # magnetic Reynolds number from compositional (Nimmo) and thermal convection (whatever is larger at each time step)

Fs2 = Flux2[0]
Fcmb2 = Flux2[1]
Fad2 = Flux2[2]
Frad2 = Flux2[3]

t_plot2 = t2/Myr

# #import label info - read in from correct row in csv
run_info = pd.read_csv('run_info3.csv',delimiter=',')

row = run_info[run_info['run']==run1]
r = row.iloc[0,1] #radius [m]
Tsolidus = row.iloc[0,2] #temperature for onset of solidification [K]
Tm0 = row.iloc[0,3] # initial mantle and core temp [K]
tstart = row.iloc[0,4]
tend = row.iloc[0,5] # max possible time of simulation [Myr]
tstep = row.iloc[0,6] #max timestep [Myr]
tsolid = row.iloc[0,7] #time at which solidifcation finishes [Myr] if tsolid == tend then core may not have finished solidifying
cond_i= int(row.iloc[0,8]) #index in array at which mantle started conducting
cond_t = t1[cond_i]/Myr #time at which mantle switched to conduction
dr=row.iloc[0,10] #cell spacing
dt =row.iloc[0,11] #T_profile output frequency

################### Main Plot #########################################
#choose model labels
model1='FK approx.'
model2='Arrhenius'
# make  log-log plot - similar to Bryson 2019

plt.figure(tight_layout=True)
plt.suptitle('Thermal evolution of a {:.0f}km asteroid \n Tm0 = {}K, Tsolidus ={}K '.format(r/1e3, Tm0, Tsolidus))

xmin=tstart

#temperatures as function of time
plt.subplot(2,1,1)
plt.semilogx(t_plot1,Tm1,label='T$_m$ - base {}'.format(model1),color='royalblue')
plt.semilogx(t_plot1,Tc1,label='T$_c$',color='firebrick')
plt.semilogx(t_plot2,Tm2,color='royalblue',linestyle='--',label=model2)
plt.semilogx(t_plot2,Tc2,color='firebrick',linestyle='--')
plt.vlines(cond_t,ymin=min(Tm1),ymax=1600,color='black',linestyle='dotted',label='conduction')
#plt.xlim([xmin,max(t_plot)])
#plt.xlabel('Time/ Myr')
#plt.xlim([xmin,500])  #use these limits when comparing runs
#plt.ylim([1400,1650]) #use these limits when comparing runs
plt.ylabel('T/K')
plt.legend(loc='upper right',ncol=2,fontsize='small')

#fluxes as function of time
plt.subplot(2,1,2)
plt.loglog(t_plot1,Fs1,label='$F_s$',color='royalblue')
plt.loglog(t_plot1,Fcmb1,label='$F_{CMB}$',color='firebrick')
plt.loglog(t_plot1,Fad1,label='$F_{ad}$',color='forestgreen')
plt.loglog(t_plot1,Frad1,label='$F_{rad}$',color='orange')

plt.loglog(t_plot2,Fs2,color='royalblue',linestyle='--')
plt.loglog(t_plot2,Fcmb2,color='firebrick',linestyle='--')
plt.loglog(t_plot2,Fad2,color='forestgreen',linestyle='--')
plt.loglog(t_plot2,Frad2,color='orange',linestyle='--')
plt.vlines(cond_t,ymin=min(Fad1),ymax=max(Frad1),color='black',linestyle='dotted',label='conduction')
plt.xlabel('Time/ Myr')
plt.ylim([1e-3,1e2])   #use these limits when comparing runs
# plt.xlim([xmin,500])   #use these limits when comparing runs
plt.ylabel('Flux/ W$m^{-2}$')
plt.legend(loc='upper right',ncol=2,fontsize='small')
#plt.savefig('Plots/Tflux_comp.png',dpi=450)

### Do same thing for Rem and f 
plt.figure(tight_layout=True)
plt.suptitle('Thermal evolution of a {:.0f}km asteroid \n Tm0 = {}K, Tsolidus ={}K'.format(r/1e3, Tm0, Tsolidus))
plt.subplot(2,1,1)
plt.loglog(t_plot1,Rem1,label=model1,color='royalblue')
plt.loglog(t_plot2,Rem2,label=model2,color='orange',linestyle='--')
#plt.xlim([xmin,max(t_plot)])
plt.hlines(10,xmin=0,xmax=t_plot2[len(Rem1)-1],color='k',linestyle='--')
#plt.xlabel('Time/Myr')
plt.ylabel('Rem')
plt.legend(loc='upper left')
plt.ylim([1,100])

plt.subplot(2,1,2)
plt.semilogx(t_plot1,f1,label=model1,color='royalblue')
plt.semilogx(t_plot2,f2,label=model2,color='orange',linestyle='--')
#plt.xlim([xmin,max(t_plot)])
plt.xlabel('Time/ Myr')
plt.ylabel('f')
plt.legend()

#plt.savefig('Plots/Remf_comp.png',dpi=450)