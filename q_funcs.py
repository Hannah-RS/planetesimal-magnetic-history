#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Expressions for individual power sources in the core during core solidification
divided by dTcdt from Nimmo (2009)
"""
from parameters import rc, rhoc, Lc, G, gc, rho_exp, rhofe_s
from fe_fes_liquidus import fe_fes_density
from heating import fe_heating
import numpy as np

def qlt(Tc,f,dTl_dP):
    """
    Latent heat release
    Parameters
    ----------
    Tc : float
        core temperature [K]
    f : float
        fractional inner core radius.
    dTl_dP : float
        pressure gradient of the liquidus [K/Pa]

    Returns
    -------
    qlt: float
        Power contribution due to release of latent heat divided by dTc/dt [W]

    """
    qlt = 4*np.pi*(f*rc)**2*Lc/(gc*dTl_dP)
    
    return qlt
    
def qgt(Tc,f,dTl_dP,Xs):
    """
    GPE release
    Parameters
    ----------
    Tc : float
        core temperature [K]
    f : float
        fractional inner core radius
    dTl_dP : float
        pressure gradient of the liquidus [K/Pa]
    Xs : float
        sulfur content [wt %]

    Returns
    -------
    qgt: float
        Power contribution due to release of GPE from release of light elements 
        when inner core solidifies [W]

    """
    rhol = fe_fes_density(Xs)*rho_exp
    drho = rhofe_s - rhol 
    qgt = -8/3*np.pi**2*G*rhol*drho*(f*rc)**4/(rhoc*gc*dTl_dP)
 
    return qgt

def qr(t):
    """
    Radiogenic heating
    
    Parameters
    ----------
    t: float
        time since beginning of simulation (differentiation of asteroid) [s]
    Tc : float
        core temperature [K]

    Returns
    -------
    Power source due to radiogenic heat production [W]

    """
    
    Mc = 4/3*np.pi*rc**3*rhoc #mass of the core [kg]
    h_fe = fe_heating(t) #internal heat generation rate from iron [J /kg /s]
    
    return Mc*h_fe