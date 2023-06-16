#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from parameters import de_fe, de_k, thalf_fe, thalf_k, f60, ppm_k , rc, rhoc, D
import numpy as np

def Qr(t,Tc):
    """
    From Table 1 in Nimmo, F. (2009). Energetics of asteroid dynamos and the role of compositional convection. 
    For source of radioactive parameters see calculations in yellow folder and refereences in parameters file

    Parameters
    ----------
    t: float
        time since beginning of simulation (differentiation of asteroid) [s]
    Tc : float
        core temperature [K]

    Returns
    -------
    Power source due to radiogenic heat production

    """
    
    Mc = 4/3*np.pi*rc**3*rhoc #mass of the core [kg]
    h_fe = de_fe*f60*np.exp(-t/thalf_fe)/thalf_fe #internal heat generation rate from iron [J /kg /s]
    h_k = de_k*ppm_k*np.exp(-t/thalf_k)/thalf_k #internal heat generation rate from potassium
    
    return Mc*(h_fe+h_k)
    