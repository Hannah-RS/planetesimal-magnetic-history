#!/usr/bin/env python3
# -*- coding: utf-8 -*-
def Er(Tc):
    """
    From Table 1 in Nimmo, F. (2009). Energetics of asteroid dynamos and the role of compositional convection. 

    Parameters
    ----------
    Tc : float
        core temperature

    Returns
    -------
    Entropy source due to radiogenic heat production

    """
    from parameters import h, rc, rho, D
    import numpy as np
    
    Mc = 4/3*np.pi*rc**3*rho #mass of the core [kg]
    
    return Mc*h/Tc*2/5*rc**2/(D**2)
    