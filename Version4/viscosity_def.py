#!/usr/bin/env python3
# -*- coding: utf-8 -*-
def viscosity(Tm):
    """
    Viscosity model from Dodds et al (2021)
    
    Parameters
    ----------
    Tm : float
        mantle temperature

    Returns
    -------
    Silicate viscosity

    """
    import numpy as np
    
    log10_eta = 64 - Tm/29 - 5*np.tanh((Tm-1625)/15)
    
    return 10**log10_eta