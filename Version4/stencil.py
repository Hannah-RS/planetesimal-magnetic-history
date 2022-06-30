#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function for creating the stencil for the solution of the diffusion equation
Requires no inputs just run once at the beginning of model for a given body to create it then use multiple times
Returns r_mat which is the required stencil based on Equation 13 in supplementary materials of Bryson (2019) but with Ts fixed and dT/dr = 0 at r=0
Left half of array has kappa_c, right half has kappa_m
"""
def cond_stencil(r,rc,dr,kappa_m,kappa_c):
    """
    

    Parameters
    ----------
    r : float
            radius of body [m]
    rc : float
        core radius [m]
    dr : float
        cell width [m]
    kappa_m : float
        thermal diffusivity of silicate [m^2 /s]
    kappa_c: float
        thermal diffusivity of core [m^2 /s]
    Returns
    -------
    r_mat : matrix (2d np array)
        stencil used for conductive dTdt

    """
    
    import numpy as np
    
    
    #create a temperature array same length as body
    n_cells = int(r/dr) #number of cells needed to span body

    # top value is surface, bottom row is CMB
    rarr = np.arange(0,r,500) # create values of cell boundaries
    rmid = (rarr[:-1]+rarr[1:])/2 #find midpoints of cells
    
    # create matrix for radial steps
    "Could make sparse later?"
    r_mat = np.zeros([n_cells,n_cells])
    
    #first row is set by assuming dT/dr=0 at centre 
    r_mat[0,0] = -1
    r_mat[0,1] = 1
    #last row is zeros as temperature of surface doesn't change
    
    for i in range(1,n_cells-1): #fill matrix but ignore first and last rows as they will need different values for bcs
    
        r_mat[i,i-1] = 1 - dr/rmid[i]
        r_mat[i,i] = -2
        r_mat[i,i+1] = 1 + dr/rmid[i]
    
    # multiply left half of array by kappa_c as multiplies core temps
    #if n_cells is odd then extend manle slightly
    r_mat_c = r_mat[:,:int(n_cells/2)]*kappa_c/dr**2
    r_mat_m = r_mat[:,int(n_cells/2):]*kappa_m/dr**2
    
    out = np.hstack([r_mat_c,r_mat_m])
    
    return out

