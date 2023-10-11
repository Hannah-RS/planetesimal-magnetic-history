#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for calculating CMB boundary layer thicknesses
"""
from parameters import kappa_c, eta_c, rhoc, alpha_c, gc, r, rc
from parameters import Ts, alpha_m, kappa, rhom, g
from viscosity_def import viscosity

def delta_l(Tm,Tcmb,Ur):
    """
    Eqn 21 in Dodds 2020 simplified to cancel out length-scale of conduction
    Parameters
    ----------
    Tm : float
        mantle temperature [K]
    Tcmb : float
        CMB temperature [K]
    Ur : float
        Urey ratio

    Returns
    -------
    mantle bottom boundary layer thickness  - Thiriet et. al. 2019 eqn. 13, 14, 16, 3
    """
    eta1 = viscosity(Tm)
    eta2 = viscosity((Tm+Tcmb)/2)
    delta_l = 0.65*abs(Tcmb-Tm)**(-1/3)*(Tm-Ts)**0.07*(r-rc)**0.21*(kappa/(alpha_m*rhom*g))**0.26*eta1**(-0.07)*eta2**(1/3) #Thiriet 2019 I think

    return delta_l

def delta_c(Tc,Tcmb):
    """
    Eqn on pg 14 in Dodds 2020 
    Parameters
    ----------
    Tc : float
        core temperature [K]
    Tcmb : float
        CMB temperature [K]
    Returns
    -------
    core cmb boundary layer thickness
    """
    return ((kappa_c*eta_c)/(rhoc*alpha_c*gc*abs(Tc-Tcmb)))**(1/3)