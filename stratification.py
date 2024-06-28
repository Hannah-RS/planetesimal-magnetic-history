#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from parameters import rhoc, rc

def volume_average(Tprofile, unstable_ind, dr):
    """
    Calculate average temperature over a volume. 
    Used in calculating removal of thermal stratification
    Parameters
    ----------
    Tprofile : float
        array of temperatures [K]
    unstable_ind : int
        indices of positions of instability
    dr: float
        grid spacing [m]

    Returns
    -------
    Tave: float
        volume average temperature [K]

    """
    rcore = np.arange(0,rc+dr,dr,dtype='float64') 
    min_unstable = int(min(unstable_ind)) 
    if min_unstable ==0: #amend unstable indices so layer mass works
        unstable_ind = unstable_ind[1:]
    r = rcore[min_unstable:] #all cells above the minimum unstable are unstable
    m = rhoc*4/3*np.pi*r**3
    mlayer = np.diff(m)
    Tmid = [(a + b) / 2 for a, b in zip(Tprofile[:-1], Tprofile[1:])]
    T = Tmid[min_unstable:]*mlayer #multiply midpoint temps by layer masses
    Tave = np.sum(T)/np.sum(mlayer)
    
    return Tave

# # test the function - should return 1
# #single unstable point
# rc=5
# Tprofile = np.array([2,3,4,5,6,7])
# unstable_ind = [4]
# if volume_average(Tprofile,unstable_ind,1)==6.5:
#     print('Test 1 passed')
# else: raise ValueError('Test 1 failed')

# # #multiple values
# rc=5
# Tprofile = np.array([2,3,4,5,6,7])
# unstable_ind = [2,3,4]
# if round(volume_average(Tprofile,unstable_ind,1),5)==5.85897:
#     print('Test 2 passed')
# else: raise ValueError('Test 2 failed')

# # multiple values over whole core
# rc=5
# Tprofile = np.array([2,3,4,5,6,7])
# unstable_ind = [0,1,2,3,4]
# if volume_average(Tprofile,unstable_ind,1)==5.7:
#     print('Test 3 passed')
# else: raise ValueError('Test 3 failed')

# #lower indices are unstable
# rc=5
# Tprofile = np.array([1,1,1,1,1,1])
# unstable_ind = [1,2,3]
# if round(volume_average(Tprofile,unstable_ind,1),5)==1:
#     print('Test 4 passed')
# else: raise ValueError('Test 4 failed')

# #lower indices including base are  unstable
# rc=5
# Tprofile = np.array([1,1,1,1,1,1])
# unstable_ind = [3,4]
# if round(volume_average(Tprofile,unstable_ind,1),5)==1:
#     print('Test 5 passed')
# else: raise ValueError('Test 5 failed')