#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

def find_593_depth(tdatalow,tdataup,t,temp):
    """
    Find the depth at which the temperature is 593 K for a given time
    Parameters
    ----------
    tdatalow : float
        Lower bound of the Ar-Ar age [Myr]
    tdataup : float
        Upper bound of the Ar-Ar age [Myr]
    t : array
        Model output time [Myr]
    temp : array
        Model mantle temperature profile [K]
    Returns
    -------
    rind : array
        Indices within the mantle at which the temperature is 593 K at
        the given Ar-Ar age. 6 indices, 3 upper and 3 lower bounds
    """
    rind = np.zeros([3,2],dtype=int)
    for i, time in enumerate(tdatalow):
        tind = np.where(t>=time)[0][0]
        rind[i,0] = np.where(temp[tind,:]<=593)[0][0] #find deepest depth, closest to 593 (smallest r)
    for i, time in enumerate(tdataup):
        tind = np.where(t>=time)[0][0]
        rind[i,1] = np.where(temp[tind,:]<=593)[0][0] #find deepest depth, closest to 593 (smallest r)   
    return rind

def find_623_cool(rind,temp,tempdt,cr_low,cr_up,rplot):
    """
    Find the cooling rate at 623K for a given depth
    Parameters
    ----------
    rind : array
        Indices within the mantle at which the temperature is 593 K at
        the given Ar-Ar age. 6 indices, 3 upper and 3 lower bounds
    temp : array
        Model mantle temperature profile [K]
    tempdt : array
        Model mantle cooling rate [K/Myr]
    cr_low : array
        Lower bound  on cooling rate at 623K [K/Myr]
    cr_up : array
        Upper bound on cooling rate at 623K [K/Myr]
    rplot : array
        Range of radii within the mantle
    Returns
    -------
    f2 : bool
        True if all cooling rates are within the bounds
    depth : array
        depths which satisfy the constraints (3,2) for 3 meteorites, lower and
        upper bounds [km]
    """
    #find time indices of 623K contour using depth midpoints
    rdiff = rind[:,0] - rind[:,1] +1
    rmidval = (rind[:,0] + rind[:,1])/2
    dTdt_check = np.zeros([3],dtype=bool)
    depth = np.zeros([3,2])
    for j in range(3): # iterate over meteorites
        tind = np.zeros([rdiff[j]],dtype=int)
        dTdt = np.zeros([rdiff[j]])
        for i in range(rdiff[j]):
            tind[i] = np.where(temp[:,rind[j,1]+i]<=623)[0][0] #find time when each depth hits 623K
            dTdt[i] = tempdt[tind[i],rind[j,1]+i] #cooling rate at 623K for each depth at correct time
        if np.any(dTdt[dTdt <= cr_up[j]]>= cr_low[j]):
            dTdt_check[j] = True
            #refine depths to those which work
            rind[j,1] = rind[j,1] + np.where(dTdt[dTdt <= cr_up[j]]>= cr_low[j])[0][0]
            rind[j,0] = rind[j,1] + np.where(dTdt[dTdt <= cr_up[j]]>= cr_low[j])[0][-1]
    if np.all(dTdt_check == True):
        f2 = True
        #find depths
        for i in range(3):
            depth[i,0] = rplot[-1] - rplot[rind[i,0]]
            depth[i,1] = rplot[-1] - rplot[rind[i,1]]
    else:
        f2 = False
    return f2, depth