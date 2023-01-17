#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function for differentiating an asteroid. Based on the process in Dodds et. al. (2021)
"""
from stencil import cond_stencil_general
from fe_fes_liquidus import fe_fes_liquidus
from scipy import sparse as sp
import numpy as np
from parameters import kc, km, ka, h0, Al0, XAl, thalf_al, rhoa, rhoc, rhom, Xs_0, cpa, Lc
def differentiation(Tint,tacc,r,dr,dt):
    """
    

    Parameters
    ----------
    Tint : float
        array of initial temperatures
    tacc: float
        accretion time after CAIs [s]
    r : float
        asteroid size [m]
    dr : float
        cell spacing [m]
    dt : float
        timesetp [s]

    Returns
    -------
    Tdiff: float
        final temperature profile after differentiation [K]
    Tdiff_profile: float
        temperature profiles at each step in differentiation [K]
    k_profile: float
        thermal conductivity profiles at each step in differentiation. Can be kc 
        (core), km (mantle), ka (undifferentiated) [W /m /K]
    Xfe: float
        proportion of iron each cell which is melted in differentiation [0 to 1]
    rho_profile: float
        density of each cell at each step in differentiation [kg m^-3]
    t_diff: float
        array of timesteps during differentiation [s]

    """
    sparse_mat = sp.dia_matrix(cond_stencil_general(r,dr))
    ncells = int(r/dr)
    Tliquidus = fe_fes_liquidus(Xs_0)

    #Initial step
    # Create arrays - column is one timestep
    k_profile = np.ones([ncells,1])*ka #don't know how long differentiation will last so append at each step
    rho_profile = np.ones([ncells,1])*rhoa
    heating = np.ones([ncells,1]) #1 if radiogenic heating i.e. mantle or undifferentiated, 0 if core
    Xfe = np.zeros([ncells,1]) #fraction of iron melted
    T = np.ones([ncells,1]) #temperature
    
    t = np.asarray([tacc])

    Tk = k_profile[:,0]*Tint
    H = h0*Al0*XAl*np.exp(-np.log(2)*t/thalf_al)

    #Calculate rhs 1/r^2dt/dr(r^2dt/dr)
    rhs = sparse_mat.dot(Tk) + H*heating[:,0]*rho_profile[:,0]
    
    #Calculate temperature change or melt change
    if np.any(Tint >= Tliquidus):
        dXfedt = rhs/(rhoc*Lc)
        dTdt = rhs/(rhoa*cpa)
        
        #no temp change where iron is melting
        T[Xfe[:,0] < 1,0] = Tint[Xfe[:,0] < 1]
        Xfe[Xfe[:,0] < 1,0] = dXfedt[Xfe[:,0] < 1]*dt
        
        #iron already melted increase the temperature
        T[Xfe[:,0] >= 1,0] = Tint[Xfe[:,0] >= 1,0] + dTdt[Xfe[:,0] >= 1]*dt
        Xfe[Xfe[:,0] >= 1,0] = Xfe[Xfe[:,0] >= 1,0]
        
    else: #no nodes are melting
        dTdt = rhs/(rhoa*cpa)
        T[:,0] = Tint + dt*dTdt
    
    # Add composition check and movement here
    #overwrite heating array?
     
     
    #Now loop
    i = 1
    #while np.any(rho_profile[:-1]==rhoa): #whilst any part except the top cell is not differentiated
    while i < 100:
        # make all the arrays one column bigger
        app_array = np.ones([ncells,1])
        T = np.append(T,app_array,1)
        k_profile = np.append(k_profile,app_array,1)
        rho_profile = np.append(rho_profile, app_array, 1)
        heating = np.append(heating, app_array,1) 
        Xfe = np.append(Xfe, app_array, 1)

        t = np.append(t,t[i-1]+dt)
        
        #Calculate radiogenic heating
        H = h0*Al0*XAl*np.exp(-np.log(2)*t[i]/thalf_al)
 
        #Calculate rhs 1/r^2dt/dr(r^2dt/dr)
        Tk = k_profile[:,i-1]*T[:,i-1]
        rhs = sparse_mat.dot(Tk) + H*heating[:,i-1]*rho_profile[:,i-1]
 
        #Calculate temperature change or melt change
        if np.any(T[:,i-1] >= Tliquidus):
            dXfedt = rhs/(rhoc*Lc)
            dTdt = rhs/(rhoa*cpa)
        
            #no temp change where iron is melting
            T[Xfe[:,i-1] < 1,i] = T[Xfe[:,i-1] < 1,i-1]
            Xfe[Xfe[:,i-1] < 1,i] = Xfe[Xfe[:,i-1] < 1,i-1] + dXfedt[Xfe[:,i-1] < 1]*dt
        
            #iron already melted increase the temperature
            T[Xfe[:,i-1] >= 1,i] = T[Xfe[:,i-1] >= 1,i-1] + dTdt[Xfe[:,i-1] >= 1]*dt
            Xfe[Xfe[:,0] >= 1,i] = Xfe[Xfe[:,0] >= 1,0]

        else: #no nodes are melting
            dTdt = rhs/(rhoa*cpa)
            T[:,i] = T[:,i-1] + dt*dTdt


        # Add composition check and movement here
        #overwrite heating array?
        i +=1


    #relabel for returning
    Tdiff = T
    Tdiff_profile = T[:,-1]
    t_diff = t
    
    return Tdiff, Tdiff_profile, k_profile, Xfe, rho_profile, t_diff