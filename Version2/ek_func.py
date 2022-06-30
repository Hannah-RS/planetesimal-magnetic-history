#!/usr/bin/env python3
# -*- coding: utf-8 -*-
def Ek():
    """
    From Table 1 in Nimmo, F. (2009). Energetics of asteroid dynamos and the role of compositional convection.

    Returns
    -------
    float, Entropy due to conductive losses

    """
    from parameters import rc, rho, k, D
    import numpy as np
    Mc = 4/3*np.pi*rc**3*rho # core mass [kg]
    
    return Mc*12*k*rc**2/(5*rho*D**4)