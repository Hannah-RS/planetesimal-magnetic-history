#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate  magnetic Reynolds numbers. For thermal convection uses Bryson (2019) from Nimmo 2009. For compositional
convection it uses three different formalisms and compares them.
1. Nimmo 2009 ucomp with dr/dt
2. Nichols 2021 - calculate flux based Rayleigh number, relate that to a convective power and use the results
from Aubert 2009
3. My own expression for the convective power in the equation from Aubert 2009 (haven't figured this out yet)
"""
import numpy as np
from parameters import G, alpha_c, rc, cpc, Omega, lambda_mag, rhofe_s, rho_eut, rhoc, gc, dr #drho here is delta_rho/rho

def Rem_therm(Fdrive,f,min_unstable):
    """
    thermally drived magnetic Reynolds number 
    Convective velocity from MAC balance (Weiss 2010, eqn 17 in Bryson 2019 supplementary)

    Parameters
    ----------
    Fdrive: float
        FCMB-Fad , heat flux available for driving convection
    f : float
        fractional size of liquid inner core
    min_unstable : float
        index of base of unstable layer in core 0 if whole core convecting

    Returns
    -------
    Magnetic reynolds number for a thermally driven dynamo

    """    
    #calculate utherm (eqn 17 in Bryson 2019)
    l = f*rc - dr*min_unstable #lengthscale over which convection can occur
    utherm = ((2*np.pi*G*alpha_c*rc*Fdrive)/(cpc*Omega))**0.5
    
    return utherm*l/lambda_mag

def Rem_comp(dfdt,f):
    """
    
    Parameters
    ----------
    dfdt : float
        rate of change of inner core radius
    f: float
        fractional inner core radius
    Returns
    -------
    compositional magnetic reynolds number

    """
    ucomp = ucomp_aubert(dfdt,f)
    Re_c = ucomp*rc*f/lambda_mag
    
    return Re_c



def ucomp_aubert(dfdt,f):
    """
    RMS velocity of compositional convection from equation 24 in Aubert 2009

    Parameters
    ----------
    dfdt : float
        rate of change of inner core radius
    f: float
        fractional inner core radius

    Returns
    -------
    ucomp : float
        rms velocity

    """
    d = rc*f
    Raq = rc*gc*abs(dfdt)*(rhofe_s - rho_eut)/(rhoc*Omega**3*d**2)
    p = 3/5*Raq
    c3 = 1.31 # from Figure 10 caption in Aubert 2009
    
    ucomp = c3*p**0.42*Omega*d
    
    return ucomp
