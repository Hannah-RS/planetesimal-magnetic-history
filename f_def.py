#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Definition of f as a function of Tc

:param Tc: float or array of floats, CMB temp [K]
:param Tc0: float, intial CMB temp [K]
:param f0: float, initial inner core size as fraction of core radius 
:param B: float, B=D^2/(2rc^2(\Delta-1))

:returns f: float, inner core size as a fraction of core radius
"""
import numpy as np

def f_calc(Tc,Tc0,f0,B):

    #print(Tc)
    out = np.sqrt(-2*B*np.log(Tc/Tc0)+f0**2)

    return out