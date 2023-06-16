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
from parameters import G, alpha_c, rc, cpc, Omega, lambda_mag, rhofe_s, rhoc, gc, dr, rho_exp, Bp_frac, mu0, r, fohm, cu, cb
from fe_fes_liquidus import fe_fes_density

def u_therm(Fdrive,f,min_unstable):
    """
    thermally driven convective velocities
    Convective velocity from MAC balance (Weiss 2010, eqn 17 in Bryson 2019 supplementary) and CIA balance (Christensen 2009)

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
    umac : float
        convective velocity for a thermally driven dynamo for MAC balance 
    ucia : float
        convective velocity  for a thermally driven dynamo for CIA balance 
    """    
    #calculate utherm (eqn 17 in Bryson 2019)
    l = f*rc - dr*min_unstable #lengthscale over which convection can occur
    umac = ((2*np.pi*G*alpha_c*rc*Fdrive)/(cpc*Omega))**0.5
    ucia = (Fdrive*alpha_c*gc/(rhoc*cpc))**(2/5)*(l/Omega)**(1/5)

    return umac, ucia

def Rem_therm(Fdrive,f,min_unstable):
    """
    thermally drived magnetic Reynolds number 
    Convective velocity from MAC balance and CIA balance

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
    Rem_mac : float
        Magnetic reynolds number for a thermally driven dynamo for MAC balance
    Rem_cia : float
        Magnetic reynolds number for a thermally driven dynamo for CIA balance
    """    
    #calculate utherm (eqn 17 in Bryson 2019)
    l = f*rc - dr*min_unstable #lengthscale over which convection can occur
    umac, ucia = u_therm(Fdrive,f,min_unstable)
    Rem_mac = umac*l/lambda_mag
    Rem_cia = ucia*l/lambda_mag
    return Rem_mac, Rem_cia

def Rem_comp(dfdt,f, Xs):
    """
    
    Parameters
    ----------
    dfdt : float
        rate of change of inner core radius
    f: float
        fractional inner core radius
    Xs : float
        sulfur content of core [wt %]
    Returns
    -------
    Rec : float
        compositional magnetic reynolds number
    Bcomp : float
        compositional magnetic field strength [T]

    """
    ucomp, p = ucomp_aubert(dfdt,f,Xs)
    Bcomp = Bp_frac*(rc/r)**3*cb*p**0.34*(rhoc*mu0*fohm)**0.5*Omega*f*rc
    Re_c = ucomp*rc*f/lambda_mag
    return Re_c, Bcomp

def ucomp_aubert(dfdt,f,Xs):
    """
    RMS velocity of compositional convection from equation 24 in Aubert 2009

    Parameters
    ----------
    dfdt : float
        rate of change of inner core radius
    f: float
        fractional inner core radius
    Xs : float
        sulfur content of core [wt %]

    Returns
    -------
    ucomp : float
        rms velocity

    """
    rhol = fe_fes_density(Xs)*rho_exp
    Raq = rc*gc*abs(dfdt)*(rhofe_s - rhol)/(rhol*Omega**3*(f*rc)**2)
    p = 3/5*Raq #convective power per unit volume
    
    ucomp = cu*p**0.42*Omega*f*rc
    
    return ucomp, p

def B_flux_therm(Fdrive,f,min_unstable):
    """
    Magnetic field strength from flux based scaling for thermal flux from Christensen 2009
    
    Parameters
    ----------
    Fdrive: float
        FCMB-Fad , heat flux available for driving convection [Wm^-2]
    f : float
        fractional size of liquid inner core
    min_unstable : float
        index of base of unstable layer in core 0 if whole core convecting

    Returns
    -------
    Bflux_ml : float
        magnetic field strength using mixing length scaling for convective velocity [T]
    Bflux_mac : float
        magnetic field strength using MAC scaling for convective velocity [T]
    Bflux_cia : float
        magnetic field strength using CIA scaling for convective velocity [T]

    """
    l = f*rc - dr*min_unstable #lengthscale over which convection can occur
    Bflux_ml = Bp_frac*(rc/r)**3*(2*mu0*fohm)**0.5*((4*np.pi*f*rc*l*alpha_c*G*rhoc**(3/2)*Fdrive)/(3*cpc))**(1/3)
    Bflux_mac = Bp_frac*(rc/r)**3*(2*mu0*fohm*l)**0.5*((4*np.pi*f*rc*alpha_c*G*rhoc**2*Omega*Fdrive)/(3*cpc))**(1/4)
    Bflux_cia = Bp_frac*(rc/r)**3*(2*mu0*fohm*rhoc*Omega**(1/5)*l**(4/5))**0.5*((4*np.pi*f*rc*alpha_c*G*Fdrive)/(3*cpc))**(3/10)
    
    return Bflux_ml, Bflux_mac, Bflux_cia
