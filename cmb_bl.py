#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for calculating CMB boundary layer thicknesses
"""
from parameters import gamma, c1, kappa_c, eta_c, rhoc, alpha_c, gc, r, rc, default
from Rayleigh_def import Rayleigh_noH

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
    mantle bottom boundary layer thickness
    """

    RanoH, d0_noH = Rayleigh_noH(Tm,default)
    if Ur > 1:
        delta_l = 0.667*(r-rc)*(gamma*abs(Tcmb-Tm)/c1)**(1.21)*RanoH**(-0.27) #eqn 26 Deschamps & Villela (2021) 
    else:
        delta_l = 0.633*(r-rc)*(gamma*abs(Tcmb-Tm)/c1)**(1.21)*RanoH**(-0.27)

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