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
from parameters import Acmb, dr, kc, rc, rhoc, cpc, Vc, Xs_0, Pc, gc
from q_funcs import Qr, Qlt, Qgt
from fe_fes_liquidus import fe_fes_liquidus_dp
import numpy as np
from Rem_calc import Rem_comp

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
    else:
        f3 = 0 #if there is an inner core treat as isothermal
        rstrat = 0
        Vconv = 4/3*np.pi*rc**3*(f**3)
        #nic_cells = round(f*rc/dr)
        
    
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
    Rem_c : float
        compositional magnetic Reynolds number
    Bcomp : float
        flux based magnetic field strength [T]

    """
    dTl_dP = fe_fes_liquidus_dp(Xs, Pc)
    Qst = rhoc*cpc*Vc
    Qrad = Qr(t)
    Ql = Qlt(Tcore[0],f,dTl_dP)
    #Qg = Qgt(Tcore[0],f,dTl_dP,Xs) #exclude as makes neglible difference

    dTcdt = (Fcmb*Acmb-Qrad)/(Qst+Ql)
    dfdt = - dTcdt/(rhoc*gc*dTl_dP*rc)

    f_new = f+dfdt*dt
    Rem_c, Bcomp = Rem_comp(dfdt,f,Xs) 
   
    return dTcdt, f_new, Rem_c, Bcomp