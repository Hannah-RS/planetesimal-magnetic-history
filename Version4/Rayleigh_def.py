#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for calculating the Rayleigh number, expression for d0 comes from substituting in the expression for Rayleigh number
so that length of the domain cancels
Also script for Rayleigh number and critical Rayleigh number for differentiation section
"""

from viscosity_def import viscosity #import viscosity model
from heating import Al_heating, AlFe_heating

from parameters import gamma, rhom, alpha_m, g, r, rc, kappa, km, Rac, Ts, default, c1, G, convect_ratio, cpm_p

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
    
def Rayleigh_calc(t,Tb,dTmdt,model=default):
    """
    

    Parameters
    ----------
    t : float
        time [s]
    Tb : float
        mantle base temperature
    dTmdt : float
        rate of temperature change of convecting mantle [K/s]

    Returns
    -------
    Rayleigh number, stagnant lid thickness

    """

    RaH, RaRob = Rayleigh_H(t,Tb,dTmdt,model=model)
    RanoH, d0 = Rayleigh_noH(Tb,model)
    
    if RaH > RanoH:
        Ram = RaH
        d0 = 0.65*(r-rc)*(gamma*(Tb-Ts))**(1.21)*RanoH**(-0.27) #eqn 26 Deschamps & Villela (2021) using average for alid
    else:
        Ram = RanoH
    
    return Ram, d0, RaH, RanoH, RaRob

def Rayleigh_noH(Tb,model=default): 
    """
    

    Parameters
    ----------
    t : float
        time [s]
    Tb : float
        mantle base temperature

    Returns
    -------
    Rayleigh number, stagnant lid thickness

    """
    eta = viscosity(Tb,model)
    d0 = (gamma/c1)**(4/3)*(Tb-Ts)*((Rac*kappa*eta)/(rhom*g*alpha_m))**(1/3) #upper bl
    Ram= rhom*g*alpha_m*abs(Tb-Ts)*(r-rc)**3/(kappa*eta)
    
    return Ram, d0
    
def Rayleigh_H(t,Tb,dTmdt,rcore = rc, model=default,Fe=False):
    """
    Rayleigh number for radiogenic heating

    Parameters
    ----------
    t : float
        time [s]
    Tb : float
        temperature at the base of the convecting region [K]
    dTmdt : float
        rate of temperature change of convecting mantle [K/s]
    rcore: float
        core radius, defaults to value in parameters [m]   
    model : str, optional
        viscosity model The default is default (set in parameters.py).

    Returns
    -------
    Ra : float
        Rayleigh number

    """
    eta = viscosity(Tb,model)

    if Fe == False: #exclude radiogenic heating from Fe
        h = Al_heating(t)
    else:
        h = AlFe_heating(t)
    
    hstar = abs(rhom*h-rhom*cpm_p*dTmdt)
    Ra = rhom**3*alpha_m*h*G*(r-rcore)**6/(km*kappa*eta) #Internally heated sphere (Schubert 2001)
    RaRob = rhom*alpha_m*hstar*g*(r-rcore)**5/(km*kappa*eta)
    
    return Ra, RaRob


def Rayleigh_differentiate(t,Tb,dTmdt,model=default):
    """
    Check for onset of differentiation

    Parameters
    ----------
    t : float
        time [s]
    Tb : float
        temperature at the base of the convecting region
    dTmdt : float
        rate of temperature change of convecting mantle [K/s]
    model : str, optional
        viscosity model The default is default (set in parameters.py).

    Returns
    -------
    RaH : float
        Rayleigh number
    d0H : float
        stagnant lid thickness
    Ra_crit : float
        critical Rayleigh number
    convect: bool
        True if cell convects, false if not

    """
   
    RaH, RaRob = Rayleigh_H(t,Tb,dTmdt,0,model,Fe=True)

    RanoH, d0_noH = Rayleigh_noH(Tb,model)

    d0H = 0.65*r*(gamma*abs(Tb-Ts))**(1.21)*RanoH**(-0.27) #eqn 26 Deschamps & Villela (2021) using average for alid
    
    Ra_crit = Rayleigh_crit(Tb)
    if (d0H/r < convect_ratio) & (RaRob>Ra_crit): #still working on this criteria
        convect = True
    else: 
        convect = False
    
    return RaRob, d0H,  Ra_crit, convect