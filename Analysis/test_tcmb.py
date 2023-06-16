#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test for Tcmb code
"""
import numpy as np
import scipy.optimize as sco
import sys
# setting path
sys.path.append('../')
from parameters import km, kc, c1, gamma, kappa_c, kappa, eta_c, alpha_m, alpha_c, rhom, rhoc, g, gc, Rac, Ts
from viscosity_def import viscosity

Tm = 1450
Tc = 1570

def Tcmb_func(Tcmb,Tm,Tc):
    """
    Equation for Tcmb, = 0 when Tcmb has been found

    Parameters
    ----------
    Tm : float
        well-mixed mantle temperature
    Tc : float
        well mixed core temperature

    Returns
    -------
    f1-f2  where f1 is eqn. 29 and f2 is eqn 26 in Dodds (2021)

    """
    #Eqn 26 - flux out of core to CMB
    delta_c = ((kappa_c*eta_c)/(rhoc*alpha_c*gc*(Tc-Tcmb)))**(1/3)
    f2 = -kc*(Tcmb-Tc)/delta_c
    
    #Eqn 28 - flux from CMB to mantle
    eta_m = viscosity(Tm)
    delta_l = (gamma/c1)**(4/3)*(Tcmb-Tm)**(4/3)*(Tm-Ts)**(-1/3)*((kappa*eta_m*Rac)/(alpha_m*rhom*g))**(1/3)
    f1 = -km*(Tm-Tcmb)/delta_l
    
    return f1 - f2
    
Tcmb_fit = sco.root_scalar(Tcmb_func,(Tm,Tc),x0=(1+1e2)*Tm,x1=(1-1e2)*Tc)
Tcmb = Tcmb_fit.root
Tcmb=np.real(Tcmb) # not sure why this is giving me a complex root, but complex part is tiny

Tcmb_func(Tcmb,Tm,Tc)


