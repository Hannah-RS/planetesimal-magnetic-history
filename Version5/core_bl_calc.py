#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function for calculating thickness of core boundary layer
"""
def core_bl(Tc,Tcmb):
    """
    

    Parameters
    ----------
    Tc : float
        core temperature [K]
    Tcmb : float
        cmb temperature [K]

    Returns
    -------
    delta_c: float
        core boundary layer thickness [m]
    """
    
    from parameters import eta_c, rhoc, kappa_c, alpha_c, gc
    import numpy as np
    a = (kappa_c*eta_c)/(rhoc*alpha_c*gc)**(1/3)
    b = np.abs(Tc-Tcmb) #use absolute value as temperature could go the other way, just temp difference matters
    c = b**(-1/3)
    #return ((kappa_c*eta_c)/(rhoc*alpha_c*gc*(Tc-Tcmb)))**(1/3)
    return a*c