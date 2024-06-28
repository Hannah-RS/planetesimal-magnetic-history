#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rayleigh numbers
"""
from viscosity_def import viscosity #import viscosity model
from heating import al_heating, alfe_heating

from parameters import frht, rhom, alpha_m, g, r, rc, kappa, km, Ts, default, G, convect_ratio, Rac

def rayleigh_crit(Tb):
    """
    Critical Rayleigh number from eqn 26 of Solomatov 1995

    Parameters
    ----------
    Tb : float
        temperature at base of convecting region [K]

    Returns
    -------
    Ra_crit : float
        critical Rayleigh number

    """
    Ra_crit = 20.9*(frht*(Tb-Ts))**4 
    if Ra_crit < Rac: #if variable critical Ra is less than analytical min ignore
        Ra_crit = Rac
    return Ra_crit
    
def rayleigh_calc(t,Tb,Ur,model=default):
    """
    Calculate all Rayleigh numbers and stagnant lid thickness
    Parameters
    ----------
    t : float
        time [s]
    Tb : float
        mantle base temperature [K]
    Ur : float
        Urey ratio
    model : str
        viscosity model default set in parameters

    Returns
    -------
    Ram : float
        model Rayleigh number (whichever of RanoH and RaH is greater)
    d0 : float
        stagnant lid thickness [m]
    RaH : float
        radiogenic heating Rayleigh number
    RanoH : float
        non-radiogenic (internal) Rayleigh number
    eta : float
        viscosity [Pas]

    """

    RaH = rayleigh_H(t,Tb,model=model)
    RanoH, eta = rayleigh_noH(Tb,model)

    if Ur > 1:
        d0 = 0.667*(r-rc)*(frht*(Tb-Ts))**(1.21)*RanoH**(-0.27) #eqn 26 Deschamps & Villela (2021) 
    else:
        d0 = 0.633*(r-rc)*(frht*(Tb-Ts))**(1.21)*RanoH**(-0.27) #eqn 26 Deschamps & Villela (2021) 
        
    if RaH > RanoH:
        Ram = RaH
    else:
        Ram = RanoH
    
    return Ram, d0, RaH, RanoH, eta

def rayleigh_noH(Tb,model=default): 
    """
    Non-radiogenic Rayleigh number 
    Eqn. 9 of Sanderson et. al. (2024)

    Parameters
    ----------
    t : float
        time [s]
    Tb : float
        mantle base temperature [K]
    model : str
        viscosity model default set in parameters

    Returns
    -------
    Ram : float
        non-radiogenic (internal) Rayleigh number
    eta : float
        viscosity [Pas]
    """
    eta = viscosity(Tb,model)
    Ram= rhom*g*alpha_m*abs(Tb-Ts)*(r-rc)**3/(kappa*eta)
    
    return Ram, eta

def rayleigh_H(t,Tb,rcore = rc, model=default,Fe=False):
    """
    Rayleigh number for radiogenic heating for an internally heated sphere 
    (Schubert 2021)

    Parameters
    ----------
    t : float
        time [s]
    Tb : float
        temperature at the base of the convecting region [K]
    rcore: float
        core radius, defaults to value in parameters [m]   
    model : str, optional
        viscosity model The default is default (set in parameters.py).

    Returns
    -------
    Ra : float
        radiogenic Rayleigh number

    """
    eta = viscosity(Tb,model)

    if Fe == False: #exclude radiogenic heating from Fe
        h = al_heating(t)
    else:
        h = alfe_heating(t)
    
    Ra = rhom**3*alpha_m*h*G*(r-rcore)**6/(km*kappa*eta) #Internally heated sphere (Schubert 2001)

    return Ra


def rayleigh_differentiate(t,Tb,Ur,model=default):
    """
    Check for onset of convection in undifferentiated body

    Parameters
    ----------
    t : float
        time [s]
    Tb : float
        temperature at the base of the convecting region [K]
    Ur : float
        Urey ratio
    model : str, optional
        viscosity model The default is default (set in parameters.py).

    Returns
    -------
    RaH : float
        Rayleigh number
    d0H : float
        stagnant lid thickness [m]
    Ra_crit : float
        critical Rayleigh number
    convect: bool
        True if cell convects, false if not
    eta : float
        viscosity [Pas]
    """
   
    RaH = rayleigh_H(t,Tb,0,model,Fe=True)

    RanoH, eta = rayleigh_noH(Tb,model)
    
    if Ur > 1:
        d0H = 0.667*r*(frht*abs(Tb-Ts))**(1.21)*RanoH**(-0.27) #eqn 26 Deschamps & Villela (2021) 
    else:
        d0H = 0.633*r*(frht*abs(Tb-Ts))**(1.21)*RanoH**(-0.27) #eqn 26 Deschamps & Villela (2021) 
    
    Ra_crit = rayleigh_crit(Tb)
    if d0H/r < convect_ratio:
        convect = True
    else: 
        convect = False

    return RaH, d0H,  Ra_crit, convect, eta