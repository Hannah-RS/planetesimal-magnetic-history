#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for calculating the Rayleigh number, expression for d0 comes from substituting in the expression for Rayleigh number
so that length of the domain cancels
"""
import numpy as np
from viscosity_def import viscosity #import viscosity model

from parameters import gamma, rhom, alpha_m, g, r, rc, kappa, Rac, Ts

def Rayleigh_calc(Tm):
    """
    

    Parameters
    ----------
    Tm : float
        mantle temperature

    Returns
    -------
    Rayleigh number, stagnant lid thickness

    """
         
    eta = viscosity(Tm)
    d0 = (gamma/8)**(4/3)*(Tm-Ts)*((Rac*kappa*eta)/(rhom*g*alpha_m))**(1/3) #upper bl
    drh = d0/(gamma*(Tm-Ts)) #lower bl
    Ram= rhom*g*alpha_m*(Tm-Ts)*(r-rc-(d0-drh))**3/(kappa*eta) 
    
    return Ram, d0
    