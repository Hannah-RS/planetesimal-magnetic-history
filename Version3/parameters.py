#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script containing all the constants and parameters required to run the model

"""
#Constants
G=6.67e-11 # gravitational constant [N kg^-2 m^2]
year=60*60*24*365 #number of seconds in a year
Myr=1e6*year

#Parameters (defaults taken from Nimmo, F. (2009). Energetics of asteroid dynamos and the role of compositional convection.)
#defaults are shown after the units in square brackets
rho=7019 #[kg m^-3] 7019
cp=835 #[J kg^-1 K^-1] 835
k=30 #[Wm^-1 K^-1] 30
alpha=9.2e-5 #[K^-1] 9.2e-5
Lh = 750000 # latent heat [J kg^-1] 750,000

#radioactivity
#assume core is pure iron with ppm amounts of K that has a negligible effect on mass
# energy generation rate per unit mass iron = de_fe * f60 * e^(-t/thalf_fe) 
# energy generation rate per unit mass potassium = de_k * ppm * e^(-t/thalf_k) 
thalf_fe=2.62*Myr #Henke 2013
thalf_k = 1.25e9*year #Wikipedia
f60 = 1e-6 # ratio of Fe60 to Fe56 Tachibana 2006 (Nimmo 2009)
ppm_k=400 # Nimmo 2004
de_fe = 4.9e9
de_k = 0.382

# Initial conditions
Tc0=2000 # intial CMB temp [K]
f0=0.01 #initial fractional inner core radius [N09 use fixed value of 0.5]

#constants for magnetic Reynolds number
eta = 2 # magnetic diffusivity [m^2 s^{-1}]
Omega = 3e-4 # rotational frequency (6 hr period) [s^{-1}]
beta=0.85

#Parameters (my choices)
#change to arrays later
phiv=1e-12# volumetric ohmic dissipation [WK^-1m^-3] [0.2-2 x1e-12 from N09]
drho=0.05 # \delta rho/rho 0.05
rc=1e5# core radius [m] 1e5
Delta=1.2 #dTm/dP/dT/dP 1.2

import numpy as np

#Derived constants
D= np.sqrt(3*cp/(2*np.pi*alpha*rho*G)) # scale height [m]

# df/dt = D^2/(2*Tc*rc^2*f*(Delta-1)dTc/dt=B dTc/dt) Equation 7 in N09
B=D**2/(2*rc**2*(Delta-1))
