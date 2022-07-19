#!/usr/bin/env python3
# -*- coding: utf-8 -*-
def dfdt_calc(t,y,Qr=True,Qs=True,Qg=True,Ql=True):
    """
    Calculate df/dt based on Equation 7 in Nimmo 2009

    Parameters
    ----------
    t: float
        time 
    y: vector consisting of Tm, Tc, f
    f : float
        fractional inner core radius
    Tm: float
        mantle temperature
    Tc : float
        core temperature
    Qr : boolean, optional
        power contribution from radiogenic heating. The default is True.
    Qs : boolean, optional
        power contribution from secular cooling The default is True.
    Qg : boolean, optional
        power contribution from compositional buoyancy. The default is True.
    Ql : boolean, optional
        power contribution from latent heat. The default is True.
    Returns
    -------
    dfdt

    """
    #separate out y
    Tm = y[0]
    Tc = y[1]
    f = y[2]
    
    
    from parameters import Tsolidus
    if Tc > Tsolidus:
        out = 0 #core not solidifying - based on 3-9 wt% S and phase diagram in Scheinberg 2016 
    else:  
        from parameters import B
        from dTcdt_def import dTcdt_calc
        out = -B*dTcdt_calc(t,y,Qr,Qs,Qg,Ql)/(Tc*f) 
    
    return out