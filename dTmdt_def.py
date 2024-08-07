#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from parameters import rhom, As, Acmb, rhoa, V, r, rc, Vm, convect_ratio
from heating import al_heating, alfe_heating
from cp_func import cp_calc_int

def dTmdt_calc(t,Tconv,d0,Flid,Fcmb):
    """
    Temperature change of convecting portion of the mantle.
    Uses Eqn. 14 of Sanderson et. al. (2024)
    
    Parameters
    ----------
    t : float
        time, [s]
    Tconv : float
        convective temp [K]
    d0 : float
        stagnant lid thickness [m]
    Flid: float
        heat flux through top of stagnant lid [W m^-2]
    Fcmb : float
        CMB heat flux [W m^-2]

    Returns
    -------
    dTmdt : float
            rate of change of mantle temperature [K/s]

    """
    if (r-d0) < rc: #i.e. lid thickness is less than mantle thickness
        Vocean = 4/3*np.pi*((r-d0)**3-rc**3)
        Alid = 4*np.pi*(r-d0)**2
    else:
        Vocean = Vm #put filler here as the output of this function won't be used
        Alid = As
    #calculate radiogenic heating 
    h = al_heating(t)
    rad = h*rhom*Vocean #radiogenic heating contribution
    cp = cp_calc_int(Tconv,False)
    return 1/(rhom*cp*Vocean)*(rad-Flid*Alid+Fcmb*Acmb)

def dTadt_calc(t,Tconv,d0,Flid,eutectic): 
    """
    Temperature change of convecting portion of undifferentiated body.
    Uses Eqn. 6 of Sanderson et. al. (2024)
    
    Parameters
    ----------
    t : float
        time, [s]
    Tconv : float
        convective temp [K]
    d0 : float
        stagnant lid thickness [m]
    Flid: float
        heat flux through top of stagnant lid [W m^-2]
    eutectic : bool
        if the initial sulfur content is eutectic composition

    Returns
    -------
    dTadt : float
            rate of change of body temperature [K/s]

    """
    if d0/r < convect_ratio: #i.e. lid thickness is less than mantle thickness
        Vocean = 4/3*np.pi*(r-d0)**3
        Alid = 4*np.pi*(r-d0)**2
    else:
        Vocean = V #put filler here as the output of this function won't be used
        Alid = As
    #calculate radiogenic heating 
    h = alfe_heating(t)
    rad = h*rhoa*Vocean #radiogenic heating contribution
    cp = cp_calc_int(Tconv,True,eutectic)
   
    return 1/(rhoa*cp*Vocean)*(rad-Flid*Alid)