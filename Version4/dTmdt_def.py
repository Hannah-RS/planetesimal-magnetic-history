#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Expression for dTm/dt from rearranging eqn 25 in Dodds (2020) and dTdt for an undifferentiated body
"""
from parameters import rhom, cpm_p, Vm, As, Ts, Acmb, km, gamma, Myr, cpa, rhoa, V
import numpy as np
from Rayleigh_def import Rayleigh_calc
from heating import Al_heating, AlFe_heating

def dTmdt_calc(t,Fs,Fcmb):
    """

    Parameters
    ----------
    t : float
        time, s
    Fs: float
        surface heat flux [W m^-2]
    Fcmb : float
        CMB heat flux

    Returns
    -------
    dTmdt : float
            rate of change of mantle temperature 

    """
    
    #calculate radiogenic heating 
    h = Al_heating(t)
    rad = h*rhom*Vm #radiogenic heating contribution
    
    
    return 1/(rhom*cpm_p*Vm)*(rad-Fs*As+Fcmb*Acmb)

def dTadt_calc(t,Fs):
    """

    Parameters
    ----------
    t : float
        time, s
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
    
    
    return 1/(rhoa*cpa*V)*(rad-Fs*As)