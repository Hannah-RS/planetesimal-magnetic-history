#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate different heat fluxes 
"""
def Flux_calc(Tm,Tc,cond=False,Tcmb =0):
    """
    

    Parameters
    ----------
    Tm : float
        mantle temperature (in conductive regime this is the base temp)
    Tc : float
        core temperature
    cond: bool
        whether mantle is conductive, default is false
    Tcmb: float,
        CMB temperature, not required if mantle is convecting but is required if mantle is conducting

    Returns
    -------
    Fs: float
        surface heat flux
    Fcmb: float
        CMB heat flux
    Fdrive: float
        heat flux avilable for driving the dynamo
    """
    from parameters import kc, alpha_c, rc, rhoc, G, cpc, Ts, km, gamma, dr, r
    import numpy as np
  
    if cond == False: #mantle is convecting
        #calculate surface flux as Fcmb is a scaled version
        #calculate lid thickness
        from Rayleigh_def import Rayleigh_calc
        Ra, d0 = Rayleigh_calc(Tm)
        
        #surface flux
        Fs = km*(Tm-Ts)/d0
        #CMB heat flux
        Fcmb = Fs*gamma*(Tc-Tm)
        
    else: #mantle is conductive Fcmb comes from eqn 26 in Dodds 2020
        from core_bl_calc import core_bl
        delta_c = core_bl(Tc,Tcmb)
    
        #Fcmb = kc*(Tc-Tcmb)/delta_c
        Fcmb = km*(Tc-Tm)/dr
        Fs = km*(Tm-Ts)/(r-rc) #approximate surface flux as a linear profile across the mantle
        
    
    #calculate adiabatic heat flux
    gc = G*4/3*np.pi*rc*rhoc #g at CMB
    Fad = kc*Tc*alpha_c*gc/cpc
    
    return Fs, Fcmb, Fcmb - Fad