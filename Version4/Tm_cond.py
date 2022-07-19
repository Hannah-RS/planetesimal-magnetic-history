#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Expression for dTm/dt from rearranging 13 in supplementary materials of Bryson (2019)
Calculates rate of change of conductive profile for the whole body and new temperature profile
"""
def Tm_cond_calc(dt,T):
    """

    Parameters
    ----------
    dt : float
        time step [s]
    T: array
        temperature profile, first value is r=0, last value is surface

    Returns
    -------
    Tm_new :  array
            new temperature profile

    """
    # set up
    import numpy as np
    import scipy.sparse as sp
         
    # calculate dTdt
    #import stencil
    r_mat = sp.load_npz('Stencil.npz')

        
    #calculate dTdt for conduction
    dTdt = r_mat.dot(T)
    
    Tnew = dTdt*dt + T
        
    return Tnew