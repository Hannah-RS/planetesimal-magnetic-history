#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Expression for dT/dt from rearranging 13 in supplementary materials of Bryson (2019)
Calculates rate of change of conductive profile for the whole body and new temperature profile
"""
# set up
import numpy as np
import scipy.sparse as sp
    
def T_cond_calc(dt,T,sparse_mat):
    """

    Parameters
    ----------
    dt : float
        time step [s]
    T: array
        temperature profile, first value is r=0, last value is surface
    sparse_mat: sparse matrix
        stencil for conductive temperature evolution

    Returns
    -------
    Tm_new :  array
            new temperature profile

    """

    #calculate dTdt for conduction

    dTdt = sparse_mat.dot(T)
    
    Tnew = dTdt*dt + T
        
    return Tnew