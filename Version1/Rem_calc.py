#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Magnetic Reynolds number based on Nimmo 2009
Uses constants from that paper, but varies core size and core solidification rate

Parameters
----------
drdt: float, rate of change of inner core radius

Returns
--------
Rem: float, magnetic Reynolds number
"""
def Rem(rc,drdt):
    import numpy as np
    beta=0.85 # a constant
    omega= 3e-4 # rotation angular frequency [s^-1] here 6hr period
    eta=2 # magnetic difusivity [m^2 s^-1]
    drho= 0.05 #delta rho/rho
    rho=7019 # density [kgm^-3]
    G= 6.67e-11 #gravitational constant [N m^2 kg^-2]
    
    return beta*omega**(-1/5)*rc**(8/5)*(4*np.pi*G*rho*drho*drdt/3)**(2/5)/eta