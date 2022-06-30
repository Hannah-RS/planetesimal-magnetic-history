#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for investigating effect of Qcmb on solidification timescale
Qcmb is 10 values log spaced between 20 GW (min value for core cooling) and 1000GW

Code runs until f=0.5 (need to alter stop_def.py to change this)
Max timestep is fixed at 0.002 Myr so get same resolution even if model runs for longer (change in solver_func)

"""
import numpy as np
import matplotlib.pyplot as plt

from parameters import Myr

#set values of Qcmb - two logspace distributions with a break
n=10 #number of runs
Qmin=np.log10(20e9) # W
Qbreak=np.log10(30e9) #W
Qmax=np.log10(1000e9) #W
Qcmb=np.concatenate((np.logspace(Qmin,Qbreak,n),np.logspace(Qbreak,Qmax,n)))

#create run array
run=np.linspace(36,36+2*n-1,2*n,dtype=int)

t_end=200 #end time in Myr
#import solver func
from solver_func import solver_func

#run solver for all values
for i in range(2*n):
    solver_func(run[i],t_end,Qcmb[i])
    
#reimport total times
t_final=np.zeros([2*n])
f_final=np.zeros([2*n]) # check stopping condition worked
dTdt_initial=np.zeros([2*n])

#import data from npz file

for i in range(2*n):
    npzfile = np.load('Results/run_{}.npz'.format(run[i]))

    t = npzfile['t'] #time in s
    f = npzfile['f'] #fractional inner core radius
    dTdt=npzfile['dTcdt']
    t_final[i]=max(t)
    f_final[i]=max(f)
    dTdt_initial[i]=dTdt[0]

#make plots
plt.figure()
plt.scatter(Qcmb,t_final/Myr,marker='x')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Qcmb /W')
plt.ylabel('Time to reach f=0.5 /Myr')
plt.savefig('Solidification_all_newQr.png')

plt.figure()
plt.scatter(Qcmb,-dTdt_initial,marker='x')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Qcmb /W')
plt.ylabel('dTdt initial')
#plt.savefig('Solidification_all_dt0p002.png')

#calculate gradients
m1=(np.log10(t_final[2*n-1])-np.log10(t_final[n-1]))/(np.log10(Qcmb[2*n-1])-np.log10(Qcmb[n-1])) # 2nd straight line section
m2=(np.log10(t_final[2*n-1])-np.log10(t_final[0]))/(np.log10(Qcmb[2*n-1])-np.log10(Qcmb[0])) # all points
m3=(np.log10(t_final[n-1])-np.log10(t_final[0]))/(np.log10(Qcmb[n-1])-np.log10(Qcmb[0])) #1st straight line section