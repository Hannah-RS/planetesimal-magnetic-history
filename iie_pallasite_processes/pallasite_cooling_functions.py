#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cooling rate functions for the pallasite workflow
"""
import numpy as np
import pandas as pd

def find_temp_depth(tempcool,tempdtlow,tempdtup,t,temp,tempdt):
    """
    Find the depth at which cooling rates at input temp match data 
    Parameters
    ----------
    tempcool : float
        Temperature at which to find cooling rate [K]
    tempdtlow : array
        Lower bound of the cooling rate (Yang et. al. 2010) [K/Myr]
    tempdtup : array
        Upper bound of the cooling rate (Yang et. al. 2010) [K/Myr]
    t : array
        Model output time [Myr]
    temp : array
        Model mantle temperature profile [K]
    tempdt : array
        Model mantle cooling rate [K/Myr]
    Returns
    -------
    rind : array
        Indices within the mantle at which the cooling rate at 800K match the Yang et. al. 2010
        cooling rate for a given meteorite
    """
    # ignore initial slow cooling when 26Al heating is strong
    minage = t[t>7][0]  #10 26Al half-lives
    rind = np.zeros([5,2],dtype=int)
    contour = np.zeros([len(t[t>minage])],dtype=int)
    cr = np.zeros([len(t[t>minage])])
    # find temp contour
    for i, tval in enumerate(t[t>minage]):
        iadd = len(t[t<=minage]) #add index to account for the time before minage
        contour[i] = np.where(temp[i+iadd,:]<=tempcool)[0][0] #find deepest depth (smallest r)
        cr[i] = tempdt[i+iadd,contour[i]] #find cooling rate at that depth and temp
    i=0
    for dTdtlow, dTdtup in zip(tempdtlow,tempdtup): #find shallowest index where cooling rate is greater than lower bound
        rvallow = contour[np.where(cr>=dTdtlow)[0]]
        rvalup = contour[np.where(cr<=dTdtup)[0]] 
        rvals = np.intersect1d(rvallow,rvalup) #find depths where cooling rate is within bounds
        if len(rvals)>0:
            rind[i,0] = np.max(rvals) #find shallowest depth
            rind[i,1] = np.min(rvals) #find deepest depth
        i+=1 
    return rind

def check_623_cool(rind,tempdtlow,tempdtup,t,temp,tempdt):
    """
    Check if cooling rates at 623K match Maurel et. al. 2019
    Parameters
    ----------
    rind : array
        Indices within the mantle at which the cooling rate at 800K match the Yang et. al. 2010
        cooling rate for a given meteorite
    tempdtlow : array
        Lower bound of the cooling rate (Yang et. al. 2010) [K/Myr]
    tempdtup : array
        Upper bound of the cooling rate (Yang et. al. 2010) [K/Myr]
    t : array
        Model output time [Myr]
    temp : array
        Model mantle temperature profile [K]
    tempdt : array
        Model mantle cooling rate [K/Myr]
    Returns
    -------
    f : bool
        True if the cooling rates at 623K match Maurel et. al. 2019
    rind: array
        Indices within the mantle at which the Maurel et. al. 2019 and Yang et. al. 2010 cooling rates are satisfied
    """
    if np.all(rind==0): #if no depths were found in previous function
        f = False
        rind = np.zeros([5,2],dtype=int)
    else:
        fcheck = [] #counter for cooling rate test
        rmid = (rind[:,0] + rind[:,1])/2 #find depth midpoints
        for rval, dTdtlow, dTdtup in zip(rmid,tempdtlow,tempdtup):
            if np.any(temp[:,int(rval)]<=623):
                tind = np.where(temp[:,int(rval)]<=623)[0][0] #find first time cool below 623K
                cr = tempdt[tind,int(rval)] #find cooling rate at that depth and temp
                if (cr >= dTdtlow) & (cr <= dTdtup): #check if cooling rate is within bounds
                    fcheck.append(True) 
                else:
                    fcheck.append(False)
            else:
                fcheck.append(False)
        if all(fcheck) == True:
            f = True
            #refine depth estimate to satisfying both cooling constraints
            rind2 = find_temp_depth(623,tempdtlow,tempdtup, t, temp, tempdt)
            #conditionally replace rind
            rind[rind2[:,0]<rind[:,0],0] = rind2[rind2[:,0]<rind[:,0],0] #shallower max depths
            rind[rind2[:,1]>rind[:,1],1] = rind2[rind2[:,1]>rind[:,1],1] #deeper min depths
        else:
            f = False
    return f, rind
