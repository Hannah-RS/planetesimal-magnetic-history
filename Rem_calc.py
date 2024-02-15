#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate  magnetic Reynolds numbers and magnetic field strengths.
    A unified buoyancy flux is calculated, which is then used to calculate a 
    convective power per unit volume. 
    This is then combined with scaling laws for magnetic field strength 
    and convective velocity.
"""
import numpy as np
from parameters import alpha_c, rc, cpc, Omega, lambda_mag, rhofe_s, rhoc, gc, \
    dr, rho_exp, mu0, r, fohm, cu, cb, Xs_eutectic, kc, Lc, icfrac
from fe_fes_liquidus import fe_fes_density

@np.vectorize
def r1(x,f):
    """
    outer solid shell inner radius

    Parameters
    ----------
    x : float
        fraction of mass in outer shell
    f : float
        fractional inner core radius for concentric inward

    Returns
    -------
    r1 : float
        inner radius of outer solid shell

    """
    r1 = rc*(1-(1-x)*(1-f**3))**(1/3)
    return r1

@np.vectorize
def r2(x,f):
    """
    inner solid core radius

    Parameters
    ----------
    x : float
        fraction of mass in outer shell
    f : float
        fractional inner core radius for concentric inward

    Returns
    -------
    r2 : float
        inner solid core radius

    """
    r2 = rc*(x*(1-f**3))**(1/3)
    return r2

def conv_power(f,dfdt,Xs,Tcore,Fcmb,solid):
    """
    Convective power per unit volume for combined thermal and buoyancy flux 
    adapted from Buffett 1996 (11) and Ruckriemen 2015 (16)

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
    comp : float
            buoyancy flux from solidifying core [kg/s]
    therm : float
            buoyancy flux from superadiabatic heat flux [kg/s]
    """
    
    if solid == True:
        #compositional buoyancy
        rhol = fe_fes_density(Xs)*rho_exp
        if Xs >= Xs_eutectic:
            drho = 0
        else: 
            drho = rhofe_s - rhol
        late = (alpha_c*rhol*Lc)/cpc #latent heat release at ICB
        comp = (late-drho)*rc*dfdt
        if comp <0:
            raise ValueError('comp<0',late/drho)
        #thermal buoyancy
        nic = round(r1(icfrac,f)/dr) #r1/dr
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
    
    return p, comp, therm

def Rem_b(f,dfdt,Xs,Tcore,Fcmb,solid,min_unstable):
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
    min_unstable : int
        minimum index of unstable layer
    Returns
    -------
    Rem : float
        magnetic Reynolds number
    Bdip_surf : float
        RMS dipole magnetic field strength at the surface [T]
    comp : float
            buoyancy flux from solidifying core [kg/s]
    therm : float
            buoyancy flux from superadiabatic heat flux [kg/s]
    """
    r2val = r2(icfrac,f)
    if r2val > min_unstable*dr:
        l = r1(icfrac,f) - r2val
    else:
        l = r1(icfrac,f) - min_unstable*dr 
    p, comp, therm = conv_power(f, dfdt, Xs, Tcore, Fcmb, solid)
    if p<0: #no dynamo
        Rem = 0
        Bdip_surf = 0 
    else:
        uconv = cu*p**0.42*Omega*l
        Rem = uconv*l/lambda_mag
        Bdip_cmb = cb*p**0.31*(fohm*mu0*rhoc)**0.5*Omega*l
        Bdip_surf = Bdip_cmb*((f*rc)/r)**3
        
    return Rem, Bdip_surf, comp, therm