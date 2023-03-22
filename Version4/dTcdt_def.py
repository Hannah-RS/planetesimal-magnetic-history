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
from parameters import Acmb, dr, kc, rc, rhoc, cpc, alpha_c, Vc, Xs_0, D, Pc
from q_funcs import Qr, Qlt, Qgt
from fe_fes_liquidus import fe_fes_liquidus_dp
import numpy as np

def dTcdt_calc(t,Fcmb,Tcore,f,Xs=Xs_0,stratification = [False,0]):
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
    Xs : float
        sulfur content [wt %]
    stratification: list
        is the unstably core stratified: bool
        index of lowest unstable cell 
    Returns
    -------
    dTcdt : float
            rate of change of core temperature (if negative core is cooling)

    """
    
    if stratification[0] == True:
        min_unstable_ind = int(stratification[1])
        f3 = -kc*(Tcore[min_unstable_ind]-Tcore[min_unstable_ind-1])/dr #calculate f3 - eqn 28 in Dodds (2020)
        rstrat = min_unstable_ind*dr
        Vconv = 4/3*np.pi*(rc**3-rstrat**3) #calculate Vconv - volume of cmb boundary layer has negligible volume 
        Tc = Tcore[min_unstable_ind]


    else:
        f3 = 0 #if there is an inner core treat as isothermal
        rstrat = 0
        Vconv = 4/3*np.pi*rc**3*(1-f**3)
        #nic_cells = round(f*rc/dr)
        Tc = Tcore[-2] #take temperature just below CMB as core convective temp
    
    Qst = rhoc*cpc*Vconv
    Qrad = Qr(t)
       
    dTcdt = (f3*4*np.pi*(rstrat)**2-Fcmb*Acmb+Qrad)/Qst
        
    return dTcdt

def dTcdt_calc_solid(t,Fcmb,Tcore,f,Xs,dt):
    """
    Temp change for solidification
    
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
    Xs : float
        sulfur content [wt %]
    dt : float
        timestep [s]
    Returns
    -------
    dTcdt : float
            rate of change of core temperature (if negative core is cooling)
    f_new : float
        new fractional inner core radius

    """
    dTl_dP = fe_fes_liquidus_dp(Xs, Pc)
    Delta = dTl_dP*(rhoc*cpc)/(alpha_c*Tcore[0])
    Qst = rhoc*cpc*Vc
    Qrad = Qr(t)
    Ql = Qlt(Tcore[0],f,Delta)
    Qg = Qgt(Tcore[0],f,Delta)

    dTcdt = (Qrad-Fcmb*Acmb)/(Qst-Ql-Qg)
    if Delta < 1:
        dfdt = -D**2/(2*Tcore[0]*f*rc**2*(1-Delta))*dTcdt
    else:
        dfdt = -D**2/(2*Tcore[0]*f*rc**2*(Delta-1))*dTcdt
    f_new = f+dfdt*dt
          
    return dTcdt, f_new