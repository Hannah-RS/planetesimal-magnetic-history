#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for calculating the Rayleigh number, expression for d0 comes from substituting in the expression for Rayleigh number
so that length of the domain cancels
Also script for Rayleigh number and critical Rayleigh number for differentiation section
"""
import numpy as np
from viscosity_def import viscosity #import viscosity model

from parameters import gamma, rhom, alpha_m, g, r, rc, kappa, km, Rac, Ts, default, c1, G, E, R, Tref, h0, Al0, XAl, thalf_al, t_transition

def Rayleigh_crit(Tb):
    """
    Critical Rayleigh number from eqn 26 of Solomatov 1995

    Parameters
    ----------
    Tb : float
        temperature at base of convecting region

    Returns
    -------
    Ra_crit : float
        critical Rayleigh number

    """
    Ra_crit = 20.9*(gamma*(Tb-Ts))**4 
    return Ra_crit
    
def Rayleigh_calc(t,Tb,model=default):
    """
    

    Parameters
    ----------
    t : float
        time [s]
    Tb : float
        mantle base temperature

    Returns
    -------
    Rayleigh number, stagnant lid thickness

    """
    if t < t_transition:
        # use radiogenic Ra
        RamH = Rayleigh_H(t,Tb,model=model)
        eta = viscosity(Tb,model)
        Ram= rhom*g*alpha_m*(Tb-Ts)*(r-rc)**3/(kappa*eta)
        d0 = 0.65*(r-rc)*(gamma*(Tb-Ts))**(1.21)*Ram**(-0.27) #eqn 26 Deschamps & Villela (2021) using average for alid
    else:
        RamH, d0 = Rayleigh_noH(Tb,model)
    
    return RamH, d0

def Rayleigh_noH(Tb,model=default): 
    """
    

    Parameters
    ----------
    t : float
        time [s]
    Tb : float
        mantle base temperature

    Returns
    -------
    Rayleigh number, stagnant lid thickness

    """
    eta = viscosity(Tb,model)
    d0 = (gamma)**(4/3)*(Tb-Ts)*((Rac*kappa*eta)/(rhom*g*alpha_m))**(1/3) #upper bl
    
    #Ram= rhom*g*alpha_m*(Tb-Ts)*(r-rc-d0)**3/(kappa*eta)
    Ram= rhom*g*alpha_m*(Tb-Ts)*(r-rc-d0)**3/(kappa*eta)
    return Ram, d0
    
def Rayleigh_H(t,Tb,rcore = rc, model=default):
    """
    Rayleigh number for radiogenic heating

    Parameters
    ----------
    t : float
        time [s]
    Tb : float
        temperature at the base of the convecting region
    rcore: float
        core radius, defaults to value in parameters [m]   
    model : str, optional
        viscosity model The default is default (set in parameters.py).

    Returns
    -------
    Ra : float
        Rayleigh number

    """
    eta = viscosity(Tb,model)
    g = 4*np.pi*r*rhom*G/3
    h = h0*Al0*XAl*np.exp(-np.log(2)*t/thalf_al)
    Ra = rhom**3*alpha_m*h*G*(r-rcore)**6/(km*kappa*eta) #Internally heated sphere (Schubert 2001)
    
    return Ra


def Rayleigh_differentiate(t,Tb,model=default):
    """
    Check for onset of differentiation

    Parameters
    ----------
    t : float
        time [s]
    Tb : float
        temperature at the base of the convecting region
    model : str, optional
        viscosity model The default is default (set in parameters.py).

    Returns
    -------
    Ra : float
        Rayleigh number
    Ra_crit : float
        critical Rayleigh number
    convect: bool
        True if cell convects, false if not

    """
   
    Ra = Rayleigh_H(t,Tb,0,model)
    Ra_crit = Rayleigh_crit(Tb)
    convect = Ra>Ra_crit
    
    return Ra, Ra_crit, convect