#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
When not solidifying:
    dTc/dt  equation 27 of Dodds et al. (2020)
    
When solidiying:
    Uses expressions for Ql and Qg from Nimmo (2009)
    
Neglects heating due to radioacitivity in core    
Can toggle on and off solidification

"""
#import constants and parameters    
from parameters import gamma, Ts, Acmb, dr, kc, rc, rhoc, cpc
from cmb_bl import delta_c
from q_funcs import Qlt, Qgt
import numpy as np

def dTcdt_calc(Fcmb,Tcore,f,solidification = False):
    """

    Parameters
    ----------
    Fcmb: float
        CMB heat flux [W m^-2]
    Tcore : array
        core temperature array [K]
    f : float
        fractional inner core radius 
    solidification: bool
        is the core solidifying, default is False

    Returns
    -------
    dTcdt : float
            rate of change of core temperature (if negative core is cooling)

    """
    #calculate f3 - eqn 28 in Dodds (2020)
    nic_cells = round(f*rc/dr)
    f3 = -kc*(Tcore[nic_cells]-Tcore[nic_cells-1])/dr
    
    #calculate Vconv - volume of cmb boundary layer has negligiblee volume 
    Vconv = 4/3*np.pi*rc**3*(1-f**3)
    Qst = rhoc*cpc*Vconv
    
    Tc = Tcore[nic_cells] #take temperature just above ICB as core convective temp
    
    if solidification == True: #core solidifying so consider buoyancy, latent heat
        dTcdt = (f3*4*np.pi*(f*rc)**2-Fcmb*Acmb)/(Qst+Qlt(Tc,f)+Qgt(Tc,f))
   
    else:
        dTcdt = (f3*4*np.pi*(f*rc)**2-Fcmb*Acmb)/Qst
     
 
        
    return dTcdt