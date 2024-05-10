#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for calculating CMB boundary layer thicknesses
"""
from parameters import kappa_c, eta_c, rhoc, alpha_c, gc, r, rc
from parameters import Ts, alpha_m, kappa, rhom, g, Rac
from viscosity_def import viscosity

def delta_l(Tm,Tcmb):
    """
    Mantle bottom CMB boundary layer thickness.
    Eqn. 11 in Sanderson et. al. (2024).
    Parameters
    ----------
    Tm : float
        mantle temperature [K]
    Tcmb : float
        CMB temperature [K]

    Returns
    -------
    mantle bottom boundary layer thickness [m] 
    """
    eta1 = viscosity(Tm)
    eta2 = viscosity((Tm+Tcmb)/2)
    delta_l = 0.65*abs(Tcmb-Tm)**(-1/3)*(Tm-Ts)**0.07*(r-rc)**0.21\
        *(kappa/(alpha_m*rhom))**0.26*(eta1/g)**(-0.07)*(eta2/gc)**(1/3) 

    return delta_l

def delta_c(Tc,Tcmb):
    """
    Core CMB boundary layer thickness.
    Eqn. 21 in Sanderson et. al. (2024).
    Parameters
    ----------
    Tc : float
        core temperature [K]
    Tcmb : float
        CMB temperature [K]
    Returns
    -------
    core cmb boundary layer thickness [m]
    """
    return ((kappa_c*eta_c*Rac)/(rhoc*alpha_c*gc*abs(Tc-Tcmb)))**(1/3)