#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Solidus and liquidus functions from Katz et al. 2003"""
import numpy as np
from fe_fes_liquidus import fe_fes_density

#Constants
G = 6.67e-11 # gravitational constant [N kg^-2 m^2]

@np.vectorize
def pmid_calc(rhom,rhoc,rc,rmid,r):
    """
    Mid mantle pressure assuming incompressible core and mantle [Pa]

    Parameters
    ----------
    rhom : float
        mantle density [kg m^-3]
    rhoc : float
        core density [kg m^-3]
    rc : float
        core radius [m]
    rmid : float
        radius in mantle where you want to calculate the pressure [m]
    r : float
        planetesimal radius [m]

    Returns
    -------
    pressure at rmid [Pa]

    """
    pmid = 2/3*np.pi*G/3*(rhom**2*(r**2-rmid**2)+2*rhom*(rhoc-rhom)*rc**3*((1/rmid)-(1/r))) 
    return pmid

def deltawater(xwater,melt):
    """
    Calculate the temperature reduction on solidus or liquidus due to water.
    Parameters
    ----------
    xwater : float
        Water content of the mantle [wt%]
    melt : float
        Melt fraction of the mantle
    Returns
    -------
    deltawater : float
        Reduction in temperature due to water [K]
    """
    #calculate effect of water - constants from Table 2 of Katz et. al. 2003
    dwater = 0.01 # bulk distribution coefficient between water and melt
    kwater = 43 # [degree C wt %^-gamma-water]
    gamma_water = 0.75 
    deltawater = kwater*(xwater/(dwater+melt*(1-dwater)))**gamma_water # Eqn 16 and 18 of Katz et. al. 2003 [degree C]
    return deltawater

def solidus(r,rc,xwater,Xs,rhom):
    """
    Calculate the mantle solidus temperature for a given water content and radius.
    Parameters
    ----------
    r : float
        Radius of the body [m]
    rc : float
        Core radius [m]
    xwater : float
        Water content of the mantle [wt%]
    Xs : float
        Core sulfur content  [wt%]
    rhom : float
        mantle density [kg m^-3]
    Returns
    -------
    Tms : float
        Mantle solidus temperature [K]
    """
    rmid = rc + (r-rc)/2
    rhoc = fe_fes_density(Xs) #core density [kg m^-3]
    pmid = pmid_calc(rhom,rhoc,rc,rmid,r)/1e9 #mid mantle pressure [GPa]
    a1 = 1085.7 # [degree C]
    a2 = 132.9 # [degree C /GPa]
    a3 = -5.1 # [degree C /GPa^2]

    #calculate effect of water
    dtwater = deltawater(xwater,0) # no melt at solidus
    Tms = a1 + a2*pmid + a3*pmid**2 - dtwater + 273 # Eqn 4 and 11 of Katz et. al. 2003 [K]
    return Tms

def liquidus(r,rc,xwater,Xs,rhom):
    """
    Calculate the mantle liquidus temperature for a given water content and radius.
    Parameters
    ----------
    r : float
        Radius of the body [m]
    rc : float
        Core radius [m]
    xwater : float
        Water content of the mantle [wt%]
    Xs : float
        Core sulfur content  [wt%]
    rhom : float
        mantle density [kg m^-3]
    Returns
    -------
    Tml : float
        Mantle liquidus temperature [K]
    """
    rmid = rc + (r-rc)/2
    rhoc = fe_fes_density(Xs) #core density [kg m^-3]
    pmid = pmid_calc(rhom,rhoc,rc,rmid,r)/1e9 #mid mantle pressure [GPa]
    c1 = 1780.0 # [degree C]
    c2 = 45 # [degree C /GPa]
    c3 = -2.0 # [degree C /GPa^2]
    #calculate effect of water
    dtwater = deltawater(xwater,1) # all melt at liquidus
    Tml = c1 + c2*pmid + c3*pmid**2 - dtwater + 273 # Eqn 10 and 13 of Katz et. al. 2003 [K]
    return Tml