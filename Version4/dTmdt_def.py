#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Expression for dTm/dt from rearranging eqn 9 in supplementary materials of Bryson (2019)
"""
from parameters import rhom, cpm_p, Vm, As, Ts, h0, Al0, XAl, thalf_al, Acmb, km, gamma, Myr
import numpy as np
from Rayleigh_def import Rayleigh_calc

def dTmdt_calc(t,Tm,Tc,Qr=True,Qs=True,Qg=True,Ql=True):
    """

    Parameters
    ----------
    t : float
        time, s
    Tm: float
        base mantle temperature
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
    dTmdt : float
            rate of change of mantle temperature 

    """

    
    #calculate lid thickness
    
    Ra, d0 = Rayleigh_calc(Tm)
    
    #calculate radiogenic heating 
    h = h0*Al0*XAl*np.exp(-np.log(2)*t/thalf_al) 
    rad = h*rhom*Vm #radiogenic heating contribution
    
    #surface flux
    Fs = km*(Tm-Ts)/d0
    
    #CMB heat flux
    Fcmb = Fs*gamma*(Tc-Tm)
    
    return 1/(rhom*cpm_p*Vm)*(rad-Fs*As+Fcmb*Acmb)