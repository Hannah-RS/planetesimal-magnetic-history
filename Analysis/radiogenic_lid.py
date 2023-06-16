#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for comparing the two radiogenic lid thicknesses
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
# setting path
sys.path.append('../')
from Rayleigh_def import Rayleigh_H, Rayleigh_noH
from parameters import Myr, Ts, gamma, r, rc, Rac

#Choose your run to compare
run=71

#load your data
npzfile = np.load('Results/run_{}.npz'.format(run))
Tcmb = npzfile['Tcmb']
Tm_mid = npzfile['Tm_mid']
Tm_conv = npzfile['Tm_conv']


t = npzfile['t'] #time in s
RaH = npzfile['RaH'] 
RanoH = npzfile['RanoH'] 


f = (gamma*(Tm_conv-Ts))

#analytical lid thickness
dan = f*(r-rc)*(Rac/RaH)**(1/6)

#empirical lid thickness from Deschamps & Villela (2021) - note this uses non radio Ra
dem = 0.65*(r-rc)*(f**1.21)/(RanoH**0.27) #will get an error here as haven't masked for Tm_conv=0 yet

#Filter so only get values when mantle convecting and scale by mantle thickness
dan = dan[Tm_conv!=0]/(r-rc)
dem = dem[Tm_conv!=0]/(r-rc)
tplot = t[Tm_conv!=0]/Myr

#plot and compare
plt.figure()
plt.plot(tplot,dan,color='black',label='Analytical lid thickness')
plt.plot(tplot,dem,color='blue',label='Empirical lid thickness',linestyle='dashed')
plt.xlabel('t/Myr')
plt.ylabel('$\delta_0$/mantle thickness')
plt.legend()
plt.ylim([0,1])
#plt.savefig('Plots/radio_d0.png')

plt.figure()
plt.semilogy(tplot,RaH[Tm_conv!=0]**(-1/6),color='black',label='RaH n=-1/6')
plt.semilogy(tplot,RanoH[Tm_conv!=0]**(-0.27),color='blue',label='RanoH n=-0.27')
plt.xlabel('t/Myr')
plt.ylabel('$Ra^n$')
plt.legend()
#plt.savefig('Plots/radio_Ra.png')
