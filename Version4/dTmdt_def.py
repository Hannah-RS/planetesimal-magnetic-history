#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Expression for dTm/dt from rearranging eqn 25 in Dodds (2020) and dTdt for an undifferentiated body
"""
from parameters import rhom, As, Ts, Acmb, km, gamma, Myr, cpa, rhoa, V, r, rc, Vm
import numpy as np
from Rayleigh_def import Rayleigh_calc
from heating import Al_heating, AlFe_heating
from cp_func import cp_calc_int

def dTmdt_calc(t,Tconv,d0,Fs,Fcmb):
    """

    Parameters
    ----------
    t : float
        time, s
    Tconv : float
        convective temp [K]
    Fs: float
        surface heat flux [W m^-2]
    Fcmb : float
        CMB heat flux

    Returns
    -------
    dTmdt : float
            rate of change of mantle temperature 

    """
    if (r-d0) < rc: #i.e. lid thickness is less than mantle thickness
        Vocean = 4/3*np.pi*((r-d0)**3-rc**3)
    else:
        Vocean = Vm #put filler here as the output of this function won't be used
    #calculate radiogenic heating 
    h = Al_heating(t)
    rad = h*rhom*Vocean #radiogenic heating contribution
    cp = cp_calc_int(Tconv,False)
    
    return 1/(rhom*cp*Vocean)*(rad-Fs*As+Fcmb*Acmb)

def dTadt_calc(t,Tconv,Fs): #not sure if this is called anywhere
    """

    Parameters
    ----------
    t : float
        time, s
    Tconv : float
        convective temp [K]
    Fs: float
        surface heat flux [W m^-2]

    Returns
    -------
    dTadt : float
            rate of change of body temperature 

    """
    
    #calculate radiogenic heating 
    h = AlFe_heating(t)
    rad = h*rhoa*V #radiogenic heating contribution
    cp = cp_calc_int(Tconv, True)
    
    return 1/(rhoa*cpa*V)*(rad-Fs*As)