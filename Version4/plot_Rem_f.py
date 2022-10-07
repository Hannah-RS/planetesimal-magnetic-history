#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot of magnetic reynolds number and core size with time
The second plot from make_plots.py
"""

#import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#choose your run
run=37

#scale time to Myr
from parameters import Myr, r, Tm0, Tsolidus

#import data from npz file
npzfile = np.load('Results/run_{}.npz'.format(run))
f = npzfile['f']
t = npzfile['t'] #time in s
Rem1 = npzfile['Rem1'] # magnetic Reynolds number from compositional (Nimmo) and thermal convection (whatever is larger at each time step)
Rem2 = npzfile['Rem2'] # magnetic Reynolds number from compositional (Nichols) and thermal convection (whatever is larger at each time step)
Ra = npzfile['Ra'] 
d0 = npzfile['d0'] 

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
cond_i= int(row.iloc[0,8]) #index in array at which mantle started conducting
cond_t = t[cond_i]/Myr #time at which mantle switched to conduction
dr=row.iloc[0,10] #cell spacing
dt =row.iloc[0,11] #T_profile output frequency

################### Main Plot #########################################
# make  log-log plot - similar to Bryson 2019

plt.figure(tight_layout=True)
plt.suptitle('Thermal evolution of a {:.0f}km asteroid \n Tm0 = {}K, Tsolidus ={}K \n run {} new viscosity'.format(r/1e3, Tm0, Tsolidus,run))
plt.subplot(2,1,1)
plt.loglog(t_plot,Rem1,label='Nimmo')
plt.loglog(t_plot,Rem2,label='Nichols')
#plt.xlim([xmin,max(t_plot)])
plt.hlines(10,xmin=0,xmax=t_plot[len(Rem1)-1],color='k',linestyle='--')
#plt.xlabel('Time/Myr')
plt.ylabel('Rem')
plt.legend(loc='upper left')
plt.ylim([1,100])

plt.subplot(2,1,2)
plt.semilogx(t_plot,f,label='f')
#plt.xlim([xmin,max(t_plot)])
plt.xlabel('Time/ Myr')
plt.ylabel('f')

#plt.savefig('Plots/Remf_run{}.png'.format(run),dpi=450)
