#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Expression for dT/dt from rearranging 13 in supplementary materials of Bryson (2019)
Calculates rate of change of conductive profile for the whole body and new temperature profile
"""
# set up
import numpy as np
import scipy.sparse as sp
from parameters import cpc
from cp_func import cp_calc_arr
from heating import Al_heating, Fe_heating
    
def Tm_cond_calc(t,dt,T,sparse_mat):
    """
    Conductive thermal evolution of the mantle
    Parameters
    ----------
    t: float
        time [s]
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
    cp = cp_calc_arr(T,False)   
    h = Al_heating(t)
    dTdt = dTdt + h/cp

    Tnew = dTdt*dt + T
    Tnew[-1] = T[-1]   #top cell of mantle is always at surface temp so need to overwrite so radiogenic heating term doesn't heat it up 
    return Tnew

def Tc_cond_calc(t,dt,T,sparse_mat,radio=True):
    """
    Conductive thermal evolution of the core
    
    Parameters
    ----------
    t: float
        time [s]
    dt : float
        time step [s]
    T: array
        temperature profile, first value is r=0, last value is surface
    sparse_mat: sparse matrix
        stencil for conductive temperature evolution

    Returns
    -------
    Tc_new :  array
            new temperature profile

    """

    #calculate dTdt for conduction

    dTdt = sparse_mat.dot(T)
    
    if radio == True:
        
        h = Fe_heating(t)
        dTdt = dTdt + h/cpc
    else: 
        pass
    
    Tcnew = dTdt*dt + T
        
    return Tcnew