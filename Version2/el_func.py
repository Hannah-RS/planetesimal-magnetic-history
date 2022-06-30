#!/usr/bin/env python3
# -*- coding: utf-8 -*-
def Elt(Tc,f):
    """
    From Table 1 in Nimmo, F. (2009). Energetics of asteroid dynamos and the role of compositional convection. but divided by dTc/dt

    Parameters
    ----------
    Tc : float
        core temperature.
    f : float
        fractional inner core radius.

    Returns
    -------
    Entropy contribution due to release of latent heat divided by dTc/dt

    """
    from parameters import rc, rho, Lh, Delta
    import numpy as np
    
    Mc = 4/3*np.pi*rc**3*rho
    
    return -3/2*Mc*f*(1-f**2)*Lh/(Tc**2*(Delta-1))