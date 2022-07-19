#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for investigating the effect of stagnant lid thickness on Rayleigh number
and temperature on both
Also compares Dodds (2021) and Bryson (2019) model
"""
import numpy as np
import matplotlib.pyplot as plt

n=20
m=20
Tm =np.linspace(1400,1600,m)
Ra_plot = np.zeros([m,n])
d0_plot = np.zeros([m,n])
iter_num = np.linspace(0,n-1,n)

from parameters import eta0, gamma, dt_gamma, alpha_n, Tms, Tml, rhom, alpha_m, g, r, rc, kappa, Rac, Ts

from viscosity_def import viscosity
eta1 = viscosity(Tm)

eta2 = eta0*np.exp(-gamma*dt_gamma)*np.exp(-alpha_n*(Tm-Tms)/(Tml-Tms))

#look at change in vicosity with temp
plt.figure()
plt.semilogy(Tm,eta1,label='Dodds 2021')
plt.semilogy(Tm,eta2,label='Bryson 2019')
plt.legend()

#look at stagnant lid thickness and Rayleigh number
d0_plot2 = (gamma/8)**(4/3)*(Tm-Ts)*((Rac*kappa*eta1)/(rhom*g*alpha_m))**(1/3)
drh_plot2 = d0_plot2/(gamma*(Tm-Ts)) #lower layer thickness
Ra_plot2 = rhom*g*alpha_m*(Tm-Ts)*(r-rc-(d0_plot2-drh_plot2))**3/(kappa*eta1) 
Ra_plot3 = Ra_plot2[Ra_plot2>1000] #look at values above critical Rayleigh number
d0_plot3 = d0_plot2[Ra_plot2>1000]
Tm_plot = Tm[Ra_plot2>1000]

plt.figure(tight_layout=True)
plt.plot(Tm_plot,d0_plot3*2/r,linestyle='--',color='black')  
plt.xlabel('Mantle temp')
plt.ylabel('stagnant lid thickness /mantle thickness')
#plt.savefig('Plots/itercomp_Ra.png')

plt.figure(tight_layout=True)
plt.semilogy(Tm_plot,Ra_plot3,label='single',linestyle='--',color='black') 
plt.legend()
plt.xlabel('Mantle temp')
plt.ylabel('Ra')

#plt.savefig('Plots/itercomp_Ra.png')