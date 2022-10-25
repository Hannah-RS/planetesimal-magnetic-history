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

def dTcdt_calc(Fcmb,Tcore,f,solidification = False, stratification = False):
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
    
    if stratification == True:
        dTdr = np.gradient(Tcore,dr)
        lnb = np.max(np.where(dTdr <= 0))
        f3 = -kc*(Tcore[lnb]-Tcore[lnb-1])/dr #calculate f3 - eqn 28 in Dodds (2020)
        rstrat = lnb*dr
        Vconv = 4/3*np.pi*(rc**3-rstrat**3) #calculate Vconv - volume of cmb boundary layer has negligible volume 
        Tc = Tcore[lnb]
    else:
        f3 = 0 #if there is an inner core treat as isothermal
        rstrat = 0
        Vconv = 4/3*np.pi*rc**3*(1-f**3)
        nic_cells = round(f*rc/dr)
        Tc = Tcore[nic_cells] #take temperature just above ICB as core convective temp
    
    Qst = rhoc*cpc*Vconv
       
    if solidification == True: #core solidifying so consider buoyancy, latent heat
        dTcdt = (f3*4*np.pi*(rstrat)**2-Fcmb*Acmb)/(Qst+Qlt(Tc,f)+Qgt(Tc,f))
   
    else:
        dTcdt = (f3*4*np.pi*(rstrat)**2-Fcmb*Acmb)/Qst
     
 
        
    return dTcdt