#!/usr/bin/env python3
# -*- coding: utf-8 -*-
def Rem_calc(f,dfdt):
    """
    Calculate magnetic Reynolds number based on results of solver.py
    Uses equation 1 from Nimmo 2009, but uses d= rc- ri= rc(1-f)
    
    Rem= \Omega^{-1/5}\beta r_c^2(1-f)^{1/5}/eta*(4/3 pi G rho drho df/dt)^2/5

    Parameters
    ----------
    f : float
        fractional inner core radius.
    dfdt : float
        rate of change of fractional inner core radius.
        

    Returns
    -------
    Rem: float
        magnetic Reynolds number

    """
    import numpy as np
    from parameters import Omega, beta, rc, eta, G, drho, rho
    
    #calculate bracketed term first
    brac = 4/3*np.pi*G*rho*drho*dfdt
    
    return Omega**(-1/5)*beta*rc**2*(1-f)**(1/5)/eta*brac**(2/5)