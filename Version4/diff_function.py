#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function for differentiating an asteroid. Based on the process in Dodds et. al. (2021)
"""
from stencil import cond_stencil_general
from scipy import sparse as sp
import numpy as np
from parameters import kc, km, ka, r, dr, h0, Al0, XAl, thalf_al
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
    sparse_mat = sp.dia_matrix(cond_stencil_general)
    ncells = int(r/dr)
    
    #Initial step
    k_profile = np.ones([ncells])*ka #don't know how long differentiation will last so append at each step
    rho_profile = np.ones([ncells])*rho_a
    heating = np.ones([ncells]) #1 if radiogenic heating i.e. mantle or undifferentiated, 0 if core
    
    Tk = k_profile*Tint
    H = h0*Al0*XAl*np.exp(-np.log(2)*t/thalf_al)
    
    #Calculate rhs 1/r^2dt/dr(1/r^2dt/dr)
    rhs = sparse_mat.dot(Tk) + H*heating*rho_profile
    
    #Below here do if statement for dtdt vs dxfedt
    
    
    #Write in appending part
    
    # Write the return part
    #In here to test the function passes
    Tdiff = 1600
    Tdiff_profile = 1600
    k_profile = 0
    t_diff = 0 
    return Tdiff, Tdiff_profile, k_profile, t_diff