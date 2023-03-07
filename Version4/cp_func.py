#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function for choosing correct specific heat capacity
"""
from parameters import Tml, Tms, Ts_fe, Tl_fe, cpa, cpa_fe, cpa_fesi, cpa_si, cpm, cpm_p
import numpy as np

def cp_calc_int(T,differentiate):
    """
    
    Calculate specific heat capacity for a single temperature
    
    Parameters
    ----------
    T : float
        temperature [K]
    differentiate : bool
        is the body differentiated

    Returns
    -------
    cp : float
        effective specific heat capacity

    """
    
    if differentiate == True:
       if T < Ts_fe:
           cp = cpa
       elif T >= Ts_fe and T <= Tl_fe and T < Tms:
           cp = cpa_fe
       elif T >= Ts_fe and T < Tl_fe and T >= Tms and T < Tml:
           cp = cpa_fesi
       elif T >= Tl_fe and T < Tms:
           cp = cpa
       elif T >= Tl_fe and T >= Tms and T < Tml:
           cp = cpa_si
       elif T >= Tl_fe and T >= Tml:
           cp = cpa
       else:
           print(f"T={T},Ts_fe={Ts_fe},Tl_fe={Tl_fe}")
           raise ValueError ('Cp scenario not coded')
        
    else:
        if T <= Tml and T >= Tms: #mantle is melting
            cp = cpm_p
        else: #mantle is solid or melted
            cp = cpm
    return cp

def cp_calc_arr(Tarr,differentiate):
    """
    
    Calculate specific heat capacity for a temperature array
    
    Parameters
    ----------
    T : float
        temperature array [K]
    differentiate : bool
        is the body differentiated

    Returns
    -------
    cp : float
        array of heat capacities

    """
    cp = np.zeros([len(Tarr)])
    if differentiate == True:
        for i, T in enumerate(Tarr):
           if T < Ts_fe:
               cp[i] = cpa
           elif T >= Ts_fe and T <= Tl_fe and T < Tms:
               cp[i] = cpa_fe
           elif T >= Ts_fe and T < Tl_fe and T >= Tms and T < Tml:
               cp[i] = cpa_fesi
           elif T >= Tl_fe and T < Tms:
               cp[i] = cpa
           elif T >= Tl_fe and T >= Tms and T < Tml:
               cp[i] = cpa_si
           elif T >= Tl_fe and T >= Tml:
               cp[i] = cpa
           else:
               print(f"T={T},Ts_fe={Ts_fe},Tl_fe={Tl_fe}")
               raise ValueError ('Cp scenario not coded')
        
    else:
        for i, T in enumerate(Tarr):
            if T <= Tml and T >= Tms: #mantle is melting
                cp[i] = cpm_p
            else: #mantle is solid or melted
                cp[i] = cpm
    return cp

