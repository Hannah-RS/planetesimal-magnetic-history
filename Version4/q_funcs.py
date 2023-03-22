#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Expressions for heat fluxes divided by dTcdt from Nimmo (2009)
"""
from parameters import rc, rhoc, Lc, D, cpc, drho, G
from heating import Fe_heating
import numpy as np

def Qlt(Tc,f,Delta):
    """
    From Table 1 in Nimmo, F. (2009). Energetics of asteroid dynamos and the role of compositional convection. but divided by dTc/dt

    Parameters
    ----------
    Tc : float
        core temperature.
    f : float
        fractional inner core radius.
    Delta : float
        dTl/dP*rho*cp/(alpha*Tc) from Nimmo (2009)

    Returns
    -------
    Power contribution due to release of latent heat divided by dTc/dt

    """
    if Delta < 1:
        Qlt = -2*np.pi*f*rc*Lc*rhoc*D**2/(Tc*(1-Delta))
    else:
        Qlt = -2*np.pi*f*rc*Lc*rhoc*D**2/(Tc*(Delta-1))
    
    return Qlt


def Qst(Tc):
    """
    From Table 1 in Nimmo, F. (2009). Energetics of asteroid dynamos and the role of compositional convection. but divided by dTc/dt
    

    Parameters
    ----------
    Tc : float
        core temperature

    Returns
    -------
    Power contribution from secular cooling divided by dTc/dt

    """
    
    Mc = 4/3*np.pi*rc**3*rhoc
    
    return -Mc*cpc*(1+2/5*rc**2/D**2)

    
def Qgt(Tc,f,Delta):
    """
    From Table 1 in Nimmo, F. (2009). Energetics of asteroid dynamos and the role of compositional convection. Qg but divided by dTc/dt

    Parameters
    ----------
    Tc : float
        core temperature.
    f : float
        fractional inner core radius
    Delta : float
        dTl/dP*rho*cp/(alpha*Tc) from Nimmo (2009)

    Returns
    -------
    Power contribution due to release of GPE from release of light elements when inner core solidifies

    """

    Mc = 4/3*np.pi*rc**3*rhoc # mass of core [kg]
    
    from F_def import F_calc
    F = F_calc(f)
    
    if Delta < 1:
        Qgt = -3*np.pi*G*rhoc*Mc*F*drho/(1-Delta)*D**2/Tc
    else: 
        Qgt = -3*np.pi*G*rhoc*Mc*F*drho/(Delta-1)*D**2/Tc
    
    return Qgt

def Qr(t):
    """
    From Table 1 in Nimmo, F. (2009). Energetics of asteroid dynamos and the role of compositional convection. Neglecting heat from potassium. 
    For source of radioactive parameters see calculations in yellow folder and refereences in parameters file

    Parameters
    ----------
    t: float
        time since beginning of simulation (differentiation of asteroid) [s]
    Tc : float
        core temperature [K]

    Returns
    -------
    Power source due to radiogenic heat production

    """
    
    Mc = 4/3*np.pi*rc**3*rhoc #mass of the core [kg]
    h_fe = Fe_heating(t) #internal heat generation rate from iron [J /kg /s]
    
    return Mc*h_fe