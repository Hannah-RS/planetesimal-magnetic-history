#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fe-FeS liquidus - linear liquidus approximation for eutectic solidification at 32% S at 1234K
liquid iron melting at 1810 K 

Pressure dependent liquidus from Buono & Walker (2011)
"""

import numpy as np

def fe_fes_liquidus_linear(Xs):
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

def fe_fes_liquidus_bw(Xs,P):
    """
    Fe-FeS liquidus from Buono & Walker (2011) - eqn 29. Valid for sub-eutectic S contents up to 10GPa. 

    Parameters
    ----------
    Xs : float
        wt % sulfur
    P : float
        pressure [GPa]       

    Returns
    -------
    Liquidus temp [K]

    """
    Mr_s = 32.07 # Pub chem [amu]
    Mr_fe = 55.84 #Pub chem  [amu]
    Xsd = Xs/100 #convert wt % to decimal

    mrr = Mr_fe/Mr_s
    x = Xsd*mrr/(1-Xsd) #mole fraction of FeS in Fe,FeS
    
    A = (-2.4724*P**4) + (28.025*P**3) + (9.1404*P**2) + (581.71*P) + 3394.8
    B = (1.7978*P**4) + (-6.7881*P**3) + (-197.69*P**2) + (-271.69*P) + (-8219.5)
    C = (-0.1702*P**4) + (-9.3959*P**3) + (163.53*P**2) + (-319.35*P) + 5698.6
    D = (-0.2308*P**4) + (7.1*P**3) + (-64.118*P**2) + 105.98*P + (-1621.9)
    E = 0.2302*P**4 + (-5.3688*P**3) + 38.124*P**2 + (-46.681*P) + 1813.8
    
    Tl = A*x**4 + B*x**3 + C*x**2 + D*x**1 + E
    
    return Tl

def fe_fes_liquidus_dp(Xs,P):
    """
    dTl/dP using Fe-FeS liquidus from Buono & Walker (2011) - eqn 29. 
    Valid for sub-eutectic S contents up to 10GPa. 

    Parameters
    ----------
    Xs : float
        wt % sulfur
    P : float
        pressure [GPa]       

    Returns
    -------
    dTl/dP 

    """
    Mr_s = 32.07 # Pub chem [amu]
    Mr_fe = 55.84 #Pub chem  [amu]
    Xsd = Xs/100 #convert wt % to decimal

    mrr = Mr_fe/Mr_s
    x = Xsd*mrr/(1-Xsd) #mole fraction of FeS in Fe,FeS
       
    F = 4*(-2.4724*P**3)+3*(28.025*P**2)+2*(9.1404*P)+(581.71)
    G = 4*(1.7978*P**3)+3*(-6.7881*P**2)+2*(-197.69*P)+(-271.69)
    H = 4*(-0.1702*P**3)+3*(-9.3959*P**2)+2*(163.53*P)+(-319.35)
    I = (-4*0.2308*P**3) + 3*7.1*P**2 + (-2*64.118*P) + 105.98
    J = 4*0.2302*P**3 + (-3*5.3688*P**2) + 2*38.124*P + (-46.681)
    
    dTl_dP = F*x**4 + G*x**3 + H*x**2 + I*x**1 + J
    
    return dTl_dP
