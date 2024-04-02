#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function for creating the stencil for the solution of the diffusion equation
Requires no inputs just run once at the beginning of model for a given body to create it then use multiple times
Returns r_mat which is the required stencil based on Equation 13 in supplementary materials of Bryson (2019) but with Ts fixed and dT/dr = 0 at r=0
Three stencils: one for core, one for mantle and one without k which can be multiplied in later
"""
from parameters import n_cells
def cond_stencil_mantle(r,rc,dr,krho):
    """
    

    Parameters
    ----------
    r : float
            radius of body [m]
    rc : float
        core radius [m]
    dr : float
        cell width [m]
    krho : float
        thermal conductivity/density of silicate [W kg^-1 K^-1 m^2]
    Returns
    -------
    r_mat : matrix (2d np array)
        stencil used for conductive dTdt for mantle

    """
    
    import numpy as np
    
    
    #create a temperature array same length as mantle
    nmcells = round((n_cells-3)*(r-rc)/r)+2 #number of cells needed to span mantle plus one extra for CMB

    # top value is surface, bottom row is CMB
    rarr = np.arange(rc-dr,r+dr,dr) # create values of cell boundaries
    rmid = (rarr[:-1]+rarr[1:])/2 #find midpoints of cells
    # create matrix for radial steps
    r_mat = np.zeros([nmcells,nmcells])

    
    #first row uses forward differences
    # r_mat[0,0] = 1 - 2*dr/rmid[0]
    # r_mat[0,1] = 2*(dr/rmid[0]-1)
    # r_mat[0,2] = 1
    #last row is zeros as temperature of surface doesn't change
    
    for i in range(1,nmcells-1): #fill matrix but ignore first and last rows as they will need different values for bcs
    
        r_mat[i,i-1] = 1 - dr/rmid[i]
        r_mat[i,i] = -2
        r_mat[i,i+1] = 1 + dr/rmid[i]
    
    # multiply array by kappa_m

    r_mat_m = r_mat*krho/dr**2
       
    return r_mat_m

def cond_stencil_core(r,rc,dr,kappa_c):
    """

    Parameters
    ----------
    r : float
            radius of body [m]
    rc : float
        core radius [m]
    dr : float
        cell width [m]
    kappa_c : float
        thermal diffusivity of core [m^2 /s]
    Returns
    -------
    r_mat : matrix (2d np array)
        stencil used for conductive dTdt for core

    """
    
    import numpy as np
    
    
    #create a temperature array same length as core
    nccells = round((n_cells-3)*(rc/r))+2 #number of cells needed to span core 
    # top value is CMB, bottom row is centre
    rarr = np.arange(0,rc+dr,dr) # create values of cell boundaries
    rmid = (rarr[:-1]+rarr[1:])/2 #find midpoints of cells
    
    # create matrix for radial steps
    r_mat = np.zeros([nccells,nccells])
    
    #first row is set by assuming dT/dr=0 at centre 
    r_mat[0,0] = -1
    r_mat[0,1] = 1
    #last row is zeros as temperature of CMB doesn't change (determined separately from flux condition)
    
    for i in range(1,nccells-1): #fill matrix but ignore first and last rows as they will need different values for bcs
    
        r_mat[i,i-1] = 1 - dr/rmid[i]
        r_mat[i,i] = -2
        r_mat[i,i+1] = 1 + dr/rmid[i]
    
    # multiply left half of array by kappa_c as multiplies core temps
    #if n_cells is odd then extend manle slightly
    r_mat_c = r_mat*kappa_c/dr**2
    
    return r_mat_c

def cond_stencil_general(r,dr):
    """

    Parameters
    ----------
    r : float
            radius of body [m]
    dr : float
        cell width [m]

    Returns
    -------
    r_mat_g : matrix (2d np array)
        stencil used for conductive dTdt 

    """
    
    import numpy as np

    # top value is CMB, bottom row is centre
    rarr = np.arange(0,r+dr,dr) # create values of cell boundaries
    rmid = (rarr[:-1]+rarr[1:])/2 #find midpoints of cells
    
    # create matrix for radial steps
    r_mat = np.zeros([n_cells,n_cells])
    
    #first row is set by assuming dT/dr=0 at centre 
    r_mat[0,0] = -1
    r_mat[0,1] = 1
    #last row is zeros as temperature of CMB doesn't change (determined separately from flux condition)
    
    for i in range(1,n_cells-1): #fill matrix but ignore first and last rows as they will need different values for bcs
    
        r_mat[i,i-1] = 1 - dr/rmid[i]
        r_mat[i,i] = -2
        r_mat[i,i+1] = 1 + dr/rmid[i]
    r_mat_g = r_mat/dr**2
    
    return r_mat_g