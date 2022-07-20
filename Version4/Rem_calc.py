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
from parameters import G, alpha_c, rc, cpc, Omega, lambda_mag, drho, rhoc, f0, gc #drho here is delta_rho/rho
from aubert import gamma_aubert

def Rem_therm(Fdrive):
    """
    thermally drived magnetic Reynolds number (equation 16 from Bryson et al. 2019)

    Parameters
    ----------
    Fdrive: float
        FCMB-Fad , heat flux available for driving convection

    Returns
    -------
    Magnetic reynolds number for a thermally driven dynamo

    """    
    #calculate utherm (eqn 17 in Bryson 2019)
 
    utherm = ((2*np.pi*G*alpha_c*rc*Fdrive)/(cpc*Omega))**0.5
    
    return utherm*rc/lambda_mag

def Rem_comp(ucomp,f):
    """
    

    Parameters
    ----------
    ucomp : float
        convective velocity from compositional convection
    f : float 
        fractional inner core radius
    Returns
    -------
    compositional magnetic reynolds number

    """
    Re_c = ucomp*rc*(1-f)/lambda_mag
    
    return Re_c

def ucomp_nimmo(t,f):
    """
    compositionally drived magnetic Reynolds number 
    Rem= ucomp*l/lambda
    (ucomp is equation 1 from Nimmo 2009 which is equivalent to equation 16 from Bryson et al. 2019 but am using Fb = delta rho dri/dt as not
     solidifying at the eutectic)

    Parameters
    ----------
    t: float
        time
    f: float
        fractional inner core radius

    Returns
    -------
    convective velocity for a solidifying inner core

    """
    
    #calculate dfdt
    "two options here - use analytical expressions in already written functions or numerically take an array of t and f and approximate"
    #calculate numerically
    dfdt = np.gradient(f,t) #this might do weird things when f is constant so use a mask to replace those values
    fmask = np.ma.masked_values(f== f0,f,True)
    dfdt[fmask]=0 #set values with initial inner core size = 0 as core is not growing
    
    #calculate ucomp (eqn 1 in Nimmo 2009 Fb = delta rho dri/dt)
    
    
    ucomp = 0.85*Omega**(-1/5)*rc*(1-f)**(1/5)*(4/3*np.pi*G*drho*rhoc*dfdt)**(2/5)
    
    return ucomp

def ucomp_aubert(p,f):
    """
    RMS velocity of compositional convection from equation 24 in Aubert 2009

    Parameters
    ----------
    p : float
        convective power density
    f: float
        fractional inner core radius

    Returns
    -------
    ucomp : float
        rms velocity

    """
    c3 = 1.31 # from Figure 10 caption in Aubert 2009
    d = rc*(1-f)
    ucomp = c3*p**0.42*Omega*d
    
    return ucomp

def p_nichols(t,f):
    """
    Expression for convective power density based on buoyancy flux Rayleigh number from Nichols 2021 (pg 8-9)
    Parameters
    ----------
    t: float
        time
    f: float
        fractional inner core radius

    Returns
    -------
    Convective power density for a compositionally driven dynamo

    """
    dfdt = np.gradient(f,t) #this might do weird things when f is constant so use a mask to replace those values
    fmask = np.ma.masked_values(f== f0,f,True)
    dfdt[fmask]=0 #set values with initial inner core size = 0 as core is not growing
    
    Qb = 4*np.pi*f**2*rc**2*drho*rhoc*rc*dfdt #buoyancy flux
    d = (1-f)*rc
    Raq = gc*Qb/(4*np.pi*Omega**3*d**4*rhoc) #  flux based Rayleigh number corrected version of eqn 15 Nichols 2021 
    
    p = gamma_aubert(f)*Raq#convective power density
    
    return p