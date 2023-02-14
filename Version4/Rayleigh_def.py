#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for calculating the Rayleigh number, expression for d0 comes from substituting in the expression for Rayleigh number
so that length of the domain cancels
Also script for Rayleigh number and critical Rayleigh number for differentiation section
"""
import numpy as np
from viscosity_def import viscosity #import viscosity model

from parameters import gamma, rhom, alpha_m, g, r, rc, kappa, km, Rac, Ts, default, c1, G, E, R, Tref, h0, Al0, XAl, thalf_al

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
    
def Rayleigh_calc(Tm,model=default):
    """
    

    Parameters
    ----------
    Tm : float
        mantle temperature

    Returns
    -------
    Rayleigh number, stagnant lid thickness

    """
         
    eta = viscosity(Tm,model)
    d0 = (gamma/c1)**(4/3)*(Tm-Ts)*((Rac*kappa*eta)/(rhom*g*alpha_m))**(1/3) #upper bl
    
    Ram= rhom*g*alpha_m*(Tm-Ts)*(r-rc-d0)**3/(kappa*eta) 
    
    return Ram, d0

def Rayleigh_differentiate(t,T,Tb,ncells,dr,model=default):
    """
    

    Parameters
    ----------
    t : float
        time [s]
    T : float
        temperature profile of the body
    Tb : float
        temperature at the base of the convecting region
    ncells: float
        number of cells in body
    dr: float
        spacing between cells [m]
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
    #r = np.linspace(dr,ncells*dr,ncells)
    eta = viscosity(Tb,model)
    g = 4*np.pi*r*rhom*G/3
    h = h0*Al0*XAl*np.exp(-np.log(2)*t/thalf_al)
    #Ra = rhom*alpha_m*r**3*g*(Tb-Ts)/(kappa*eta) #Robuchon & Nimmo 2011
    #Ra = rhom**2*alpha_m*h*r**5/(km*kappa*eta) #Plane layer heated from within (Schubert 2001)
    Ra = rhom**3*alpha_m*h*G*r**6/(km*kappa*eta) #Internally heated sphere (Schubert 2001)
    Ra_crit = Rayleigh_crit(Tb)
    convect = Ra>Ra_crit
    
    return Ra, Ra_crit, convect