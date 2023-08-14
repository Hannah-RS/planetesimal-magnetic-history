#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate  magnetic Reynolds numbers and magnetic field strengths.
There are two implementations:
    1. A unified buoyancy flux is calculated, which is then used to calculate a convective power per unit volume. 
    This is then combined with scaling laws for magnetic field strength and convective velocity.
    2. A separate implementation for thermal and compositional dynamos which also calculates the CIA reynolds number.
    This implementation will be removed in future but is currently kept in for comparison.
    For thermal convection uses Bryson (2019) from Nimmo 2009. For compositional
    convection it uses three different formalisms and compares them. Nichols 2021 - calculate flux based Rayleigh number, relate that to a convective power and use the results
    from Aubert 2009
"""
import numpy as np
from parameters import G, alpha_c, rc, cpc, Omega, lambda_mag, rhofe_s, rhoc, gc, dr, rho_exp, Bp_frac, mu0, r, fohm, cu, cb
from parameters import Xs_eutectic, kc
from fe_fes_liquidus import fe_fes_density

#new implementation
def conv_power(f,dfdt,Xs,Tcore,Fcmb,solid):
    """
    Convective power per unit volume for combined thermal and buoyancy flux adapted from Nichols 2021 
    and Ruckriemen 2015

    Parameters
    ----------
    f : float
        fractional inner core radius
    dfdt : float
        rate of change of inner core radius
    Xs : float
        core sulfur content [wt %]
    Tcore : float
        temperature array of core [K]
    Fcmb : float
        CMB heat flux [Wm^-2]
    solid : bool
        whether the core is solidifying or not

    Returns
    -------
    p : float
        convective power per unit volume as defined in equation 17 in Aubert 2009
    """
    
    if solid == True:
        #compositional buoyancy
        rhol = fe_fes_density(Xs)*rho_exp
        if Xs >= Xs_eutectic:
            drho = 0
        else: 
            drho = rhofe_s - rhol
        comp = drho*rc*dfdt
        #thermal buoyancy
        nic = round(f*rc/dr)
        Ficb = -kc*(Tcore[nic]-Tcore[nic-1])/dr #worried if this will return anything if index is wrong
        Fad = kc*gc*alpha_c*Tcore[nic-1]/cpc
         
    else: #boundary of convecting region is at CMB
        comp = 0
        Ficb = Fcmb
        Fad = kc*alpha_c*gc*Tcore[-2]/cpc
    
    therm = alpha_c/cpc*(Ficb-Fad)
    #convert total buoyancy to convective power per unit volume
    buoy = 4*np.pi*f**2*rc**2*(therm+comp)
    Raq = gc*buoy/(4*np.pi*rhoc*Omega**3*(f*rc)**4)
    p = 3/5*Raq #convective power per unit volume
    
    return p

def Rem_b(f,dfdt,Xs,Tcore,Fcmb,solid):
    """
    Calculation of dipole field strength on the surface and magnetic Reynolds number
    Based on scaling laws of Aubert 2009, Davidson 2013 and Davies 2022
    
    Parameters
    ----------
    f : float
        fractional inner core radius
    dfdt : float
        rate of change of inner core radius
    Xs : float
        core sulfur content [wt %]
    Tcore : float
        temperature array of core [K]
    Fcmb : float
        CMB heat flux [Wm^-2]
    solid : bool
        whether the core is solidifying or not

    Returns
    -------
    Rem : float
        magnetic Reynolds number
    Bdip_surf : float
        RMS dipole magnetic field strength at the surface [T]
    """
    p = conv_power(f, dfdt, Xs, Tcore, Fcmb, solid)
    if p<0: #no dynamo
        Rem = 0
        Bdip_surf = 0 
    else:
        uconv = cu*p**0.42*Omega*f*rc
        Rem = uconv*f*rc/lambda_mag
        Bdip_cmb = 0.23*p**0.31*(fohm*mu0*rhoc)*Omega*f*rc
        Bdip_surf = Bdip_cmb*((f*rc)/r)**3
        
    return Rem, Bdip_surf

#old implementation
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
