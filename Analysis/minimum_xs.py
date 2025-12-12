#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate minimum Xs for a given planetesimal size and critical melt fraction
"""
#%% Import numpy modules
import numpy as np
import sys
# setting path
sys.path.append('../')
from fe_fes_liquidus import fe_fes_liquidus_bw, central_pressure, fe_fes_density
from parameters import rho_exp
from solidus_calc import solidus, liquidus

#%% Calculation
def min_xs_func(r,rcr,phi,xwater,solidusval='katz'):
    """
    Calculate minimum core sulfur content for a given set of planetesimal parameters
    Parameters
    ----------
    r : float
        Planetesimal radius [m]
    rcr : float
        Core radius fraction
    phi : float
        Critical melt fraction
    xwater : float
        Water content in silicate [wt %]
    solidusval : str, optional
        Solidus model to use, by default 'katz' for Katz et. al. 2003
    Returns
    -------
    minS : float
        Minimum core sulfur content [wt %]
    """
    rc = rcr*r #core radius [m]

    #Create S array
    Xs = np.linspace(0,40,100)
    #Calculate core density
    rhoc = fe_fes_density(Xs)*rho_exp
    rhom = 3000
    # Calculate solidus
    Xsnom = 15 #chosen Xs in Tms and Tml makes 0.01 K difference
    if solidusval == 'katz':
        Tms = solidus(r,rc,xwater,Xsnom,rhom) 
        Tml = liquidus(r,rc,xwater,Xsnom,rhom)
    else:
        Tms = 1400
        Tml = 1800
    #pressure
    
    P = central_pressure(rhom,rhoc,r,rc) #pressure at centre [GPa]
    Tphi = Tms+phi*(Tml-Tms)

    bw = fe_fes_liquidus_bw(Xs,P)   
    #find lowest Xs for a given silicate melting
    minS=Xs[bw<Tphi][0]

    return minS
#%% Test a specific case
#r = np.array([200e3]) #m
#rcr = 0.1 #fractional core radius
#phi = np.array([0.3]) #critical melt fraction
#xwater = np.array([0.07]) #wt %
#minS = min_xs_func(r,rcr,phi,xwater)

#print(f'The minimum initial core sulfur content is {minS:.2f} wt %')
