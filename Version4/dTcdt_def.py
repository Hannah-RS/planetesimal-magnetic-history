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
from q_funcs import Qlt, Qgt, Qr
import numpy as np

def dTcdt_calc(t,Fcmb,Tcore,f,solidification = False, stratification = [False,0]):
    """

    Parameters
    ----------
    t : float
        time [s]
    Fcmb: float
        CMB heat flux [W m^-2]
    Tcore : array
        core temperature array [K]
    f : float
        fractional inner core radius 
    solidification: bool
        is the core solidifying, default is False
    stratification: list
        is the unstably core stratified: bool
        index of lowest unstable cell 
    Returns
    -------
    dTcdt : float
            rate of change of core temperature (if negative core is cooling)

    """
    
    if stratification[0] == True:
        stratification[1] = min_unstable_ind
        f3 = -kc*(Tcore[min_unstable_ind]-Tcore[min_unstable_ind-1])/dr #calculate f3 - eqn 28 in Dodds (2020)
        rconv = min_unstable_ind*dr
        Vconv = 4/3*np.pi*(rc**3-rconv**3) #calculate Vconv - volume of cmb boundary layer has negligible volume 
        Tc = Tcore[min_unstable_ind]


    else:
        f3 = 0 #if there is an inner core treat as isothermal
        rstrat = 0
        Vconv = 4/3*np.pi*rc**3*(1-f**3)
        #nic_cells = round(f*rc/dr)
        Tc = Tcore[-2] #take temperature just below CMB as core convective temp
    
    Qst = rhoc*cpc*Vconv
    Qrad = Qr(t)
       
    if solidification == True: #core solidifying so consider buoyancy, latent heat
        dTcdt = (f3*4*np.pi*(rstrat)**2-Fcmb*Acmb+Qrad)/(Qst+Qlt(Tc,f)+Qgt(Tc,f))
   
    else:
        dTcdt = (f3*4*np.pi*(rstrat)**2-Fcmb*Acmb+Qrad)/Qst

 
        
    return dTcdt