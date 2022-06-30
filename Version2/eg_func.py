#!/usr/bin/env python3
# -*- coding: utf-8 -*-
def Egt(Tc,f):
    """
    From Table 1 in Nimmo, F. (2009). Energetics of asteroid dynamos and the role of compositional convection. but divided by dTc/dt

    Parameters
    ----------
    Tc : float
        core temperature.
    f : float
        fractional inner core radius

    Returns
    -------
    Entropy contribution due to release of GPE from release of light elements when inner core solidifies

    """
    from parameters import rho, rc, drho, D, Delta, G
    import numpy as np
    
    Mc = 4/3*np.pi*rc**3*rho # mass of core [kg]
    
    from F_def import F_calc
    F = F_calc(f)
    
    return - 3*np.pi*G*rho*Mc*F*drho*1/(Delta-1)*D**2/Tc**2