#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Expressions for heat fluxes divided by dTcdt from Nimmo (2009)
"""
from parameters import rc, rhoc, rhofe_l, rhofe_s, rho_eut, Mr_fe, Mr_s, Lc, D, cpc, G, gc
from heating import Fe_heating
import numpy as np

def Qlt(Tc,f,dTl_dP):
    """
    From Table 1 in Nimmo, F. (2009). Energetics of asteroid dynamos and the role of compositional convection. but divided by dTc/dt

    Parameters
    ----------
    Tc : float
        core temperature.
    f : float
        fractional inner core radius.
    dTl_dP : float
        dTl/dP 

    Returns
    -------
    Power contribution due to release of latent heat divided by dTc/dt

    """
    Qlt = 4*np.pi*(f*rc)**2*Lc/(gc*dTl_dP)
    
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

    
def Qgt(Tc,f,dTl_dP,Xs):
    """
    From Table 1 in Nimmo, F. (2009). Energetics of asteroid dynamos and the role of compositional convection. Qg but divided by dTc/dt

    Parameters
    ----------
    Tc : float
        core temperature.
    f : float
        fractional inner core radius
    dTl_dP : float
        dTl/dP pressure gradient of the liquidus
    Xs : float
        sulfur content wt %

    Returns
    -------
    Power contribution due to release of GPE from release of light elements when inner core solidifies

    """
    Xsd = Xs/100
    Mrr = 1+Mr_fe/Mr_s
    rhol = Xsd*Mrr*rho_eut + (1-Xsd*Mrr)*rhofe_l
    drho = rhofe_s - rhol 
    Qgt = 8/3*np.pi**2*G*rhol*drho*(f*rc)**4/(rhoc*gc*dTl_dP)
    
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