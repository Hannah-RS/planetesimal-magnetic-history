#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Expression for dTc/dt from Nimmo 2009, used for recreating Figure 1
Not for use with an ODE solver as the arguments are written differently
Negative sign in final expression has been omitted as assume it is a COOLING rate and we are just interested in the magnitude

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
def dTcdt(phiv,f0,rc):
    import numpy as np
    
    #Parameters (taken from Nimmo 2009 - N09)
    rho=7019 #[kg m^-3]
    cp=835 #[J kg^-1 K^-1]
    k=30 #[Wm^-1 K^-1]
    alpha=9.2e-5 #[K^-1]
    drho=0.05 # \delta rho/rho 
    rc=1e5# core radius [m] 
    Delta=1.2 #dTm/dP/dT/dP 
    Tc0=2000 #[K]
    Tc=2000 #[K] no time variation
    G=6.67e-11 # gravitational constant [N kg^-2 m^2]
    
    D= np.sqrt(3*cp/(2*np.pi*alpha*rho*G)) # scale height [m]
    B=D**2/(2*rc**2*(Delta-1)) # constant to simplify expression for f in terms of Tc

    from F_def import F_calc
    from f_def import f_calc
    
    f=f_calc(Tc,Tc0,f0,B)
    F=F_calc(f)
    return Tc**2/F*(Delta-1)*phiv/(3*np.pi*rho**2*drho*G*D**2)