#!/usr/bin/env python3
# -*- coding: utf-8 -*-
def Qst(Tc):
    """
    From Table 1 in Nimmo, F. (2009). Energetics of asteroid dynamos and the role of compositional convection. but divided by dTc/dt
    

    Parameters
    ----------
    Tc : float
        core temperature

    Returns
    -------
    Power contribution from secular cooling divided by dTc/dt

    """
    from parameters import cpc, rc, rhoc, D
    import numpy as np
    
    Mc = 4/3*np.pi*rc**3*rhoc
    
    return -Mc*cpc*(1+2/5*rc**2/D**2)