#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plots of lid thickness and Rayleigh number - first set of exploratory plots in make_plots.py
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
t = npzfile['t'] #time in s
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

xmin=tstart
####################### make a plot to check Rayleigh number

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