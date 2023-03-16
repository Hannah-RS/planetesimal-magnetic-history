#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fe-FeS liquidus - linear liquidus approximation for eutectic solidification at 32% S at 1234K
liquid iron melting at 1810 K 

Pressure dependent liquidus from Buono & Walker (2011)
"""

import numpy as np

def fe_fes_liquidus(Xs):
    """
    Fe-FeS liquidus - linear liquidus approximation for eutectic solidification at 32% S at 1234K
    liquid iron melting at 1810 K 
    see Kathryn's thesis pg 180
    Parameters
    ----------
    Xs : float
        wt % sulfur

    Returns
    -------
    Liquidus temp [K]

    """
    return 1810-18*Xs

def fe_fes_liquids_bw(Xs,P):
    """
    

    Parameters
    ----------
    Xs : float
        wt % sulfur
    P : float
        pressure [GPa] at centre 
      

    Returns
    -------
    Liquidus temp [K]

    """
    Mr_s = 32.07 # Pub chem [amu]
    Mr_fe = 55.84 #Pub chem  [amu]
    Xsd = Xs/100 #convert wt % to decimal

    mrr = Mr_fe/Mr_s
    molefrac_fes = Xsd*mrr*(1-Xsd*(1+mrr))**(-1)

    x = molefrac_fes
    
    
    A = -2.4724*P**4 + 28.025*P**3 + 9.1404*P**2 + 581.71*P + 3394.8
    B = 1.7978*P**4 - 6.7881*P**3 - 197.69*P**2 - 271.69*P - 8219.5
    C = -0.1702*P**4 - 9.3959*P**3 + 163.53*P**2 - 319.35*P + 5698.6
    D = -0.2308*P**4 + 7.1*P**3 - 64.118*P**2 + 105.98*P - 1621.9
    E = 0.2302*P**4 - 5.3688*P**3 + 38.124*P**2 - 46.681*P + 1813.8
    
    Tl = A*x**4 + B*x**3 + C*x**2 + D*x**1 + E
    
    return Tl