#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Expressions for heat fluxes divided by dTcdt from Nimmo (2009)
"""
from parameters import rc, rhoc, Lc, Delta, D, cpc, drho, G
import numpy as np

def Qlt(Tc,f):
    """
    From Table 1 in Nimmo, F. (2009). Energetics of asteroid dynamos and the role of compositional convection. but divided by dTc/dt

    Parameters
    ----------
    Tc : float
        core temperature.
    f : float
        fractional inner core radius.

    Returns
    -------
    Power contribution due to release of latent heat divided by dTc/dt

    """

    
    Mc = 4/3*np.pi*rc**3*rhoc
    
    return -3/2*Mc*f*Lc*D**2/(rc**2*Tc*(Delta-1))


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

    
def Qgt(Tc,f):
    """
    From Table 1 in Nimmo, F. (2009). Energetics of asteroid dynamos and the role of compositional convection. Qg but divided by dTc/dt

    Parameters
    ----------
    Tc : float
        core temperature.
    f : float
        fractional inner core radius

    Returns
    -------
    Power contribution due to release of GPE from release of light elements when inner core solidifies

    """

    Mc = 4/3*np.pi*rc**3*rhoc # mass of core [kg]
    
    from F_def import F_calc
    F = F_calc(f)
    
    return -3*np.pi*G*rhoc*Mc*F*drho/(Delta-1)*D**2/Tc