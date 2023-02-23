#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Radiogenic heating equations
"""
import numpy as np
from parameters import h0, Al0, XAl, thalf_al

def Al_heating(t):
    """
    Radiogenic heating by 26Al

    Parameters
    ----------
    t : float
        time after CAIs [s]

    Returns
    -------
    h : float
        internal heating rate [W /kg]

    """
    h = h0*Al0*XAl*np.exp(-np.log(2)*t/thalf_al)
    
    return h