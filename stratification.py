#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for calculating the average temperature over a volume
"""
import numpy as np
from parameters import rhoc, rc

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
    rcore = np.arange(0,rc,dr)    
    min_unstable = int(min(unstable_ind)) 
    if min_unstable ==0: #amend unstable indices so layer mass works
        #min_unstable = 1
        unstable_ind = unstable_ind[1:]
    max_unstable = int(max(unstable_ind))
    r = rcore[unstable_ind]
    r = np.append(r,(max_unstable+1)*dr) # add one radius larger 
    m = rhoc*4/3*np.pi*r**3
    mlayer = np.diff(m)
    Tmid = [(a + b) / 2 for a, b in zip(Tprofile[:-1], Tprofile[1:])]
    
    if min_unstable == 0: #append bottom mass
        mlayer = np.insert(mlayer,0,m[0])
    else:
        pass
    T = Tmid[min_unstable:]*mlayer #multiply midpoint temps by layer masses
    Tave = np.sum(T)/np.sum(mlayer)
   
    return Tave

# # test the function - should return 1
#single unstable point
# rc=5
# Tprofile = np.array([2,3,4,5,6,7])
# unstable_ind = [4]
# if volume_average(Tprofile,unstable_ind,1)==6.5:
#     pass
# else: raise ValueError('Test failed')

# # #multiple values
# rc=5
# Tprofile = np.array([2,3,4,5,6,7])
# unstable_ind = [2,3,4]
# if round(volume_average(Tprofile,unstable_ind,1),5)==5.85897:
#     pass
# else: raise ValueError('Test failed')

# # multiple values over whole core
# rc=5
# Tprofile = np.array([2,3,4,5,6,7])
# unstable_ind = [0,1,2,3,4]
# if volume_average(Tprofile,unstable_ind,1)==5.7:
#     pass
# else: raise ValueError('Test failed')