#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function for differentiating an asteroid. Based on the process in Dodds et. al. (2021)
"""
def differentiation(Tint,tacc,dr):
    """
    

    Parameters
    ----------
    Tint : float
        array of initial temperatures
    tacc: float
        accretion time after CAIs [s]
    dr : float
        cell spacing [m]

    Returns
    -------
    Tdiff: float
        final temperature profile after differentiation
    Tdiff_profile: float
        temperature profiles at each step in differentiation
    k_profile: float
        thermal conductivity profiles at each step in differentiation. Can be kc 
        (core), km (mantle), ka (undifferentiated)
    t_diff: float
        array of timesteps during differentiation

    """
    
    return Tdiff, Tdiff_profile, k_profile, t_diff