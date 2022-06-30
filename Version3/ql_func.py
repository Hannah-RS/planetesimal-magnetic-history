#!/usr/bin/env python3
# -*- coding: utf-8 -*-
def Qlt(Tc,f):
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
    Power contribution due to release of latent heat divided by dTc/dt

    """
    from parameters import rc, rho, Lh, Delta, D
    import numpy as np
    
    Mc = 4/3*np.pi*rc**3*rho
    
    return -3/2*Mc*f*Lh*D**2/(rc**2*Tc*(Delta-1))