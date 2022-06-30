#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code to recreate the inner core lines in Figure 1 in Nimmo 2009 
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt

# array of parameters
Tc=2000 #[K]
rho=7019 # [kg m^-3]
phi=np.array([3e-12,3e-13,3e-14]) #dissipation [W kg^-1]
phiv=phi*rho/Tc #volumetric Ohmic dissipation [W m^-3 K^-1]
n=100 # number of radii to sample
rc=np.linspace(1e3,300e3,n) #core radius [m]
f0=0.5

year=60*60*24*365 #number of seconds in a year
Myr=1e6*year

from dtcdt_fig1 import dTcdt
"""
dTcdt(phiv,f0,rc)
Parameters
----------
phiv: float,
     Volumetric Ohmic dissipation
f0 : float
    intial inner core fraction.
rc: float
    core radius, [m]

Returns
-------
float
   dTc/dt

"""

dTc=np.zeros([3,n])
for i in range(3):
    dTc[i,:]=dTcdt(phiv[i],f0,rc)

dT_plot=np.log10(dTc*Myr) #convert to K/Myr rather than K/s
#dT_plot=dTc

#create figure
plt.figure()
for i in range(3):
    plt.plot(rc/1e3,dT_plot[i,:],label='$\Phi$={}'.format(phi[i]))
plt.xlabel('Core radius/ km')
plt.ylabel('Log$_{10}$(Cooling rate K/Myr)')
plt.legend()
    
