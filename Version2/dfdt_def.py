#!/usr/bin/env python3
# -*- coding: utf-8 -*-
def dfdt_calc(f,Tc,dTcdt):
    """
    Calculate df/dt based on Equation 7 in Nimmo 2009

    Parameters
    ----------
    f : float
        fractional inner core radius
    Tc : float
        core temperature
    Tc : float
        rate of change of core temperature
    Returns
    -------
    dfdt

    """
    from parameters import B
    
    return -B*dTcdt/(Tc*f)