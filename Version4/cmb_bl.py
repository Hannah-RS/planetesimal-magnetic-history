#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for calculating CMB boundary layer thicknesses
"""
from parameters import gamma, c1, Rac, kappa_c, eta_c, rhoc, alpha_c, gc, r, rc
from Rayleigh_def import Rayleigh_calc

def delta_l(t,Tm,Tcmb):
    """
    Eqn 21 in Dodds 2020 simplified to cancel out length-scale of conduction
    Parameters
    ----------
    t : float
        time [s]
    Tm : float
        mantle temperature [K]
    Tcmb : float
        CMB temperature [K]
    Returns
    -------
    mantle bottom boundary layer thickness
    """

    Ra, d0, RaH, RanoH = Rayleigh_calc(t,Tm)
    if RaH > RanoH:
        delta_l = 0.65*(r-rc)*(gamma*abs(Tcmb-Tm)/c1)**(1.21)*RanoH**(-0.27) #eqn 26 Deschamps & Villela (2021) using average for alid
    else:
        delta_l = (r-rc)*(gamma*abs(Tcmb-Tm)/c1)**(4/3)*(Ra/Rac)**(-1/3) #eqn 21 from Dodds (2021)

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