#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code for investigating the two methods of calculating heat fluxes
Uses results from run 16, which output a temp profile very 0.1Myr
"""
import numpy as np
import matplotlib.pyplot as plt

import sys
# setting path
sys.path.append('../')
from parameters import Myr, km, dr, Ts,r 

run = 16

npzfile = np.load('../Results/run_{}.npz'.format(run))
Tm = npzfile['Tm_base'] 
Tm_surf = npzfile['Tm_surf'] 
T_profile = npzfile['T_profile']
t = npzfile['t'] #time in s
Flux = npzfile['Flux']
Ra = npzfile['Ra'] 
d0 = npzfile['d0'] 

Fs = Flux[0]
Fcmb = Flux[1]
Fad = Flux[2]
Frad = Flux[3]

tplot = t/Myr
cond_i = 8359 #index when switch from convection to conduction

tplot2 = np.arange(1,np.round(tplot[cond_i],1),0.01) #for shorter arrays
#tplot2 = np.arange(1,5,0.1) #for shorter arrays

#calculate d0 just at required times
from Rayleigh_def import Rayleigh_calc
Ra_short, d0_short = Rayleigh_calc(T_profile[:len(tplot2),int(r/(2*dr))+1])

#check they agree
plt.figure()
plt.plot(tplot[:cond_i],2*d0[:cond_i]/r,label='full')
plt.plot(tplot2,2*d0_short/r,label='short')
plt.xlabel('Time/Myr')
plt.ylabel('Lid thickness/ mantle thickness')
plt.legend()

#calculate the heatflux across the stagnant lid
# flux = km(Tm-Ts)/d0
T_iso = np.zeros(len(tplot2))
Fapprox= np.zeros(len(tplot2))
for i in range(len(tplot2)): #10Myr in 0.1 Myr steps
    nlid = int(d0_short[i]/dr) #1. Calculate number of cells in lid, take
    #2. Calculate index of Tprofile to sample
    T_iso[i] = T_profile[i,-nlid-2] #make sure definitely below lid in case of rounding error
    #3. Calculate flux
    Fapprox[i] = km*(T_iso[i]-Ts)/d0_short[i]

#How does it differ depending on number of cells approximated over?
n=5
Fvary=np.zeros([len(tplot2),n])
for i in range(n):
    Fvary[:,i] = km*(T_profile[:len(tplot2),-2-i]-Ts)/((1+i)*dr)

#compare with model Fs, F_last, F_approx
plt.figure()
plt.plot(tplot2,Fapprox,label='approx')
for i in range(n):
    plt.plot(tplot2,Fvary[:,i],label='{} cells below'.format(1+i))
#plt.plot(tplot,Fs,label='model')
plt.legend()
plt.xlabel('Time/Myr')
plt.ylabel('Surface flux')
plt.title('dr = {} m'.format(dr))
plt.xlim([0,5])
#plt.savefig('Plots/flux_approximations2.png')