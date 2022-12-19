#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function for creating the stencil for the solution of the diffusion equation
Requires no inputs just run once at the beginning of model for a given body to create it then use multiple times
Returns r_mat which is the required stencil based on Equation 13 in supplementary materials of Bryson (2019) but with Ts fixed and dT/dr = 0 at r=0
Left half of array has kappa_c, right half has kappa_m
"""
def cond_stencil_mantle(r,rc,dr,kappa_m):
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
    Returns
    -------
    r_mat : matrix (2d np array)
        stencil used for conductive dTdt for mantle

    """
    
    import numpy as np
    
    
    #create a temperature array same length as body
    n_cells = int(r/(2*dr))+1 #number of cells needed to span half body plus one extra for CMB


    # top value is surface, bottom row is CMB
    rarr = np.arange(r/2-dr,r,dr) # create values of cell boundaries
    rmid = (rarr[:-1]+rarr[1:])/2 #find midpoints of cells
    
    # create matrix for radial steps
    r_mat = np.zeros([n_cells,n_cells])

    
    #first row is zero as assuming CMB temp is constant (this is determined separately from flux balance)
    #last row is zeros as temperature of surface doesn't change
    
    for i in range(1,n_cells-1): #fill matrix but ignore first and last rows as they will need different values for bcs
    
        r_mat[i,i-1] = 1 - dr/rmid[i]
        r_mat[i,i] = -2
        r_mat[i,i+1] = 1 + dr/rmid[i]
    
    # multiply array by kappa_m

    r_mat_m = r_mat*kappa_m/dr**2
       
    return r_mat_m

def cond_stencil_core(r,rc,dr,kappa_c):
    """
    Assume core radius is half of body radius

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
    
    
    #create a temperature array same length as half
    n_cells = int(r/(2*dr)) #number of cells needed to span half body 

    # top value is CMB, bottom row is centre
    rarr = np.arange(0,r/2,dr) # create values of cell boundaries
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
    
    # multiply left half of array by kappa_c as multiplies core temps
    #if n_cells is odd then extend manle slightly
    r_mat_c = r_mat*kappa_c/dr**2
    
    return r_mat_c