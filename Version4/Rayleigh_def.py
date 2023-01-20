#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for calculating the Rayleigh number, expression for d0 comes from substituting in the expression for Rayleigh number
so that length of the domain cancels
Also script for Rayleigh number and critical Rayleigh number for differentiation section
"""
import numpy as np
from viscosity_def import viscosity #import viscosity model

from parameters import gamma, rhom, alpha_m, g, r, rc, kappa, km, Rac, Ts, default, c1, G, E, R, Tref, h0

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

def Rayleigh_differentiate(T,Tb,ncells,dr,model=default):
    """
    

    Parameters
    ----------
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
    #Ra = rhom*alpha_m*r**3*g*(Tb-Ts)/(kappa*eta) #Robuchon & Nimmo 2011
    #Ra = rhom**2*alpha_m*h0*r**5/(km*kappa*eta) #Plane layer heated from within (Schubert 2001)
    Ra = rhom**3*alpha_m*h0*G*r**6/(km*kappa*eta) #Internally heated sphere (Schubert 2001)
    
    Ra_crit = 20.9*((E/(R*Tref**2))*(Tb-Ts))**4 #2.18 Dodds thesis 
    
    #convect = (Ra-Ra_crit) >= -0.0001*Ra 
    convect = Ra>Ra_crit
    return Ra, Ra_crit, convect