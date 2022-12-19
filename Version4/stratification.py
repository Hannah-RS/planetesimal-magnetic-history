#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for calculating the average temperature over a volume
"""
import numpy as np
from parameters import rhoc

def volume_average(Tprofile, unstable_ind, dr):
    """
    

    Parameters
    ----------
    Tprofile : float
        array of temperatures
    unstable_ind : int
        indices of positions of instability
    dr: float
        grid spacing

    Returns
    -------
    Tave: float
        volume average temperature

    """
        
    
    min_unstable = min(unstable_ind)
    max_unstable = max(unstable_ind)
    r = unstable_ind*dr
    
    if min_unstable==0: #if central node is unstable that is same as one node up being unstable
        min_unstable = 1
    else:
        r = np.insert(r,0,(min_unstable-1)*dr) # add one radius smaller
    
    r = np.append(r,(max_unstable+1)*dr) # add one radius larger     
    m = rhoc*4/3*np.pi*r**3
    mlayer = np.diff(m)
    T = (Tprofile[min_unstable:max_unstable+2]+Tprofile[min_unstable-1:max_unstable+1])*mlayer/2 #use temp of midpoints
    Tave = np.sum(T)/np.sum(mlayer)
    
    return Tave

# # test the function - should return 1
# #single unstable point
# Tprofile = np.ones([10])
# unstable_ind = [2]
# if volume_average(Tprofile,unstable_ind,1)==1:
#     pass
# else: raiseValueError('Test failed')

# #multiple points
# Tprofile = np.ones([10])
# unstable_ind = [2,3]
# if volume_average(Tprofile,unstable_ind,1)==1:
#     pass
# else: raiseValueError('Test failed')

# #including the base
# Tprofile = np.ones([10])
# unstable_ind = [0,1]
# if volume_average(Tprofile,unstable_ind,1)==1:
#     pass
# else: raiseValueError('Test failed')

