#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot central temperature as a function of time for 4 values of 60Fe
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
# setting path
sys.path.append('../')
from parameters import Myr

#load run 
runs = [8,7,6,5]
fe60 = [9e-7,6e-7,4e-7,1e-7]
lines = ['dotted','dashed','-.','-']

plt.figure()
for i, run in enumerate(runs):
    #import data from npz file
    npzfile = np.load('../Results_combined/run_{}.npz'.format(run))
    Tdiff = npzfile['Tdiff']
    tdiff = npzfile['t_diff']
    Tc1 = Tdiff[0,:]
    
    npzfile = np.load('Results_combined/run_{}.npz'.format(run))
    Tc2= npzfile['Tc'] 
    t = npzfile['t'] #time in s
    
    Tcall = np.append(Tc1,Tc2)
    tall = np.append(tdiff,t)/Myr
    
    plt.plot(tall,Tcall,linestyle=f'{lines[i]}',label=f'{fe60[i]}',color='black')
    
plt.xlabel('t/Myr')
plt.ylabel('core central temperature/K')
plt.legend(title='$^{60}Fe/^{56}Fe$')
plt.xlim([0,25])
plt.savefig('Plots/fe60.png')