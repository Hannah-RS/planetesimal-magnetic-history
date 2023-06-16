#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot of temperature profiles with depth final set of plots in make_plots.py
"""
#import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#choose your run
run=38

#scale time to Myr
import sys
# setting path
sys.path.append('../')
from parameters import Myr, r, Tm0, Tsolidus

#import data from npz file
npzfile = np.load('../Results_combined/run_{}.npz'.format(run))

T_profile = npzfile['T_profile']
t = npzfile['t'] #time in s


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
dr=row.iloc[0,10] #cell spacing
dt =row.iloc[0,11] #T_profile output frequency



#################  plot temperature profiles 
plt.figure()
rplot= np.arange(0,r,dr)
n = np.shape(T_profile)[0]
l = len(t_plot)
n_plot = 4 #how many plots do you want
for i in range(n_plot):
    plt.plot(rplot/1e3, T_profile[i*int(n/n_plot),:],label='{:.0f} Myr'.format(t_plot[i*int(l/n_plot)])) #approximate temp profile times
plt.xlabel('Distance from centre of asteroid /km')
plt.ylabel('Temperature / K')
#plt.xlim([80,100])
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