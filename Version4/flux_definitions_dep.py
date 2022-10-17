#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
List of different functions for calculating the CMB flux. All written to have same arguments so they
can all be put into scipy.optimize.root_scalar so in some cases arguments are unused
"""
import numpy as np
from parameters import km, kc, dr, kappa_c, eta_c,rhoc, alpha_c, gc, gamma, c1, Rac, r
from viscosity_def import viscosity
from Rayleigh_def import Rayleigh_calc

def f1_conductive(Tm,Tc,Tcmb):
    """
    Conductive flux from CMB into mantle

    Parameters
    ----------
    Tm : float
        mantle temperature one cell above CMB
    Tc : float
        core temperature one cell below CMB
    Tcmb : float
        CMB temperature

    Returns
    -------
    f1 - eqn 23 of Dodds et al. (2020)

    """
    return -km*(Tm-Tcmb)/dr

def f1_convective(Tm,Tc,Tcmb):
    """
    Convective flux from CMB into mantle

    Parameters
    ----------
    Tm : float
        mantle temperature one cell above CMB
    Tc : float
        core temperature one cell below CMB
    Tcmb : float
        CMB temperature

    Returns
    -------
    f1 - eqn 29 of Dodds et al. (2020)

    """
    eta_m = viscosity(Tm)
    Ra, d_lid = Rayleigh_calc(Tm)
    d = r/2-d_lid #length scale for convection = length of mantle - lid thickness
    delta_l = d*(gamma*abs(Tcmb-Tm)/c1)**(4/3)*(Ra/Rac)**(-1/3) #eqn 21 #use absolute value as just interested in thickness
    
    return -km*(Tm-Tcmb)/delta_l

def f2_conductive(Tm,Tc,Tcmb):
    """
    Conductive flux from core into CMB

    Parameters
    ----------
    Tm : float
        mantle temperature one cell above CMB
    Tc : float
        core temperature one cell below CMB
    Tcmb : float
        CMB temperature

    Returns
    -------
    f2 - eqn 24 of Dodds et al. (2020)

    """
    return -kc*(Tcmb-Tc)/dr

def f2_convective(Tm,Tc,Tcmb):
    """
    Conductive flux from core into CMB

    Parameters
    ----------
    Tm : float
        mantle temperature one cell above CMB
    Tc : float
        core temperature one cell below CMB
    Tcmb : float
        CMB temperature

    Returns
    -------
    f2 - eqn 26 of Dodds et al. (2020)

    """
    delta_c = ((kappa_c*eta_c)/(rhoc*alpha_c*gc*abs(Tc-Tcmb)))**(1/3) #use absolute value as just interested in thickness
    
    return -kc*(Tcmb-Tc)/delta_c

def flux_balance(Tcmb,Tc, Tm, f1,f2):
    return f1(Tm,Tc,Tcmb)-f2(Tm,Tc,Tcmb) #returns 0 when fluxes balance