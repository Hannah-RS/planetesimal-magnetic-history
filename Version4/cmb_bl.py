#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for calculating CMB boundary layer thicknesses
"""
import numpy as np
from parameters import gamma, c1, kappa, Rac, g, alpha_m, rhom, Ts, kappa_c, eta_c, rhoc, alpha_c, gc, r, rc
from viscosity_def import viscosity
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
    eta = viscosity(Tm)
    Ra, d0, RaH, RanoH = Rayleigh_calc(t,Tm)
    delta_l = (r-rc)*(gamma*abs(Tcmb-Tm)/c1)**(4/3)*(Ra/Rac)**(-1/3)

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