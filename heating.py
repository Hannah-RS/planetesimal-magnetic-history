#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Radiogenic heating equations
"""
import numpy as np
from parameters import h0Al, Al0, XAl_a, XAl_d, thalf_al, h0Fe, Fe0, XFe_a, XFe_d, thalf_fe

def Al_heating(t):
    """
    Radiogenic heating by 26Al in a silicate mantle

    Parameters
    ----------
    t : float
        time after CAIs [s]

    Returns
    -------
    h : float
        internal heating rate [W /kg]

    """
    h = h0Al*Al0*XAl_d*np.exp(-np.log(2)*t/thalf_al)
    
    return h

def Fe_heating(t):
    """
    Radiogenic heating by 60Fe in an iron core

    Parameters
    ----------
    t : float
        time after CAIs [s]

    Returns
    -------
    h : float
        internal heating rate [W /kg]

    """
    h = h0Fe*Fe0*XFe_d*np.exp(-np.log(2)*t/thalf_fe)
    
    return h

def AlFe_heating(t):
    """
    Radiogenic heating by 26Al and 60Fe (undifferentiated material)

    Parameters
    ----------
    t : float
        time after CAIs [s]

    Returns
    -------
    h : float
        internal heating rate [W /kg]

    """
    Al_heat= h0Al*Al0*XAl_a*np.exp(-np.log(2)*t/thalf_al)
    Fe_heat = h0Fe*Fe0*XFe_a*np.exp(-np.log(2)*t/thalf_fe)
    
    h = Fe_heat + Al_heat

    return h
