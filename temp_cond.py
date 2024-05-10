#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# set up
from parameters import cpc
from cp_func import cp_calc_arr
from heating import al_heating, fe_heating
    
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
        temperature profile, first value is r=rc, last value is surface [K]
    sparse_mat: sparse matrix
        stencil for conductive temperature evolution

    Returns
    -------
    Tm_new :  array
            new mantle temperature profile [K]

    """

    #calculate dTdt for conduction

    cpdTdt = sparse_mat.dot(T) #lhs of 16 divided by rhom
    cp = cp_calc_arr(T,False)   
    h = al_heating(t)
    dTdt = cpdTdt/cp + h/cp

    Tnew = dTdt*dt + T
    Tnew[0] = Tnew[0] - h*dt/cp[0] #no radiogenic heating at cmb
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
        temperature profile, first value is r=0, last value is rc [K]
    sparse_mat: sparse matrix
        stencil for conductive temperature evolution

    Returns
    -------
    Tc_new :  array
            new core temperature profile [K]

    """

    #calculate dTdt for conduction

    dTdt = sparse_mat.dot(T)
    
    if radio == True:
        
        h = fe_heating(t)
        dTdt = dTdt + h/cpc
    else: 
        pass
    
    Tcnew = dTdt*dt + T
        
    return Tcnew