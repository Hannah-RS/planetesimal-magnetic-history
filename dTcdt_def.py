#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#import constants and parameters    
from parameters import Acmb, dr, kc, rc, rhoc, cpc, Vc, Xs_0, Pc, gc
from q_funcs import qlt, qr, qgt
from fe_fes_liquidus import fe_fes_liquidus_dp
import numpy as np
from rem_calc import rem_b
from heating import fe_heating

def dTcdt_calc(t,Fcmb,Tcore,f,Xs=Xs_0,stratification = [False,0]):
    """
    Core temperature change prior to solidification for convecting portion of the core.
    Includes the possibility of core thermal stratification.
    Uses Eqn. 20 of Sanderson et. al. (2024)
    Also calls function to calculate magnetic Reynolds number and field strengths.
    
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
            rate of change of core temperature (if negative core is cooling) [K/s]
    Rem : float
        magnetic Reynolds number
    Bdip_cmb : float
        dipole magnetic field strength at the surface [T]
    comp : float
            buoyancy flux from solidifying core [kg/s]
    therm : float
            buoyancy flux from superadiabatic heat flux [kg/s]
    qcore : np.ndarray
            heat sources in the core, qr=radiogenic, qs=secular cooling
            ql= latent heat release, qg=gpe release [qr, qs, ql, qg] [W]
            
    """
    
    if stratification[0] == True:
        min_unstable_ind = int(stratification[1])
        f3 = -kc*(Tcore[min_unstable_ind]-Tcore[min_unstable_ind-1])/dr #calculate f3 - eqn 28 in Dodds (2020)
        rstrat = min_unstable_ind*dr
        Vconv = 4/3*np.pi*(rc**3-rstrat**3) #calculate Vconv - volume of cmb boundary layer has negligible volume 
    else:
        f3 = 0 #if there is an inner core treat as isothermal
        min_unstable_ind = 0
        rstrat = 0
        Vconv = 4/3*np.pi*rc**3*(f**3)
        
    
    qst = rhoc*cpc*Vconv #secular cooling
    qrad = rhoc*Vconv*fe_heating(t) #radiogenic heating
       
    dTcdt = (f3*4*np.pi*(rstrat)**2-Fcmb*Acmb+qrad)/qst
    #calculate magnetic field
    Rem, Bdip_cmb, comp, therm = rem_b(f, 0, Xs, Tcore, Fcmb, False,min_unstable_ind) #dfdt = 0 for  non-solidifying   
    #calculate fluxes to return
    qcore =np.array([qrad,qst*dTcdt,0,0])
    
    return dTcdt, Rem, Bdip_cmb, comp, therm, qcore

def dTcdt_calc_solid(t,Fcmb,Tcore,f,Xs,dt):
    """
    Core temperature change during solidification.
    Uses Eqn. 25 of Sanderson et. al. (2024)
    Also calls function to calculate magnetic Reynolds number and field strengths.
    
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
            rate of change of core temperature (if negative core is cooling) [K/s]
    f_new : float
        new fractional inner core radius
    Rem : float
        magnetic Reynolds number
    Bdip_cmb : float
        dipole magnetic field strength at the surface [T]
    comp : float
            buoyancy flux from solidifying core [kg/s]
    therm : float
            buoyancy flux from superadiabatic heat flux [kg/s]
    qcore : np.ndarray
            heat sources in the core, qr=radiogenic, qs=secular cooling
            ql= latent heat release, qg=gpe release [qr, qs, ql, qg] [W]
    """
    dTl_dP = fe_fes_liquidus_dp(Xs, Pc)
    qst = rhoc*Vc*cpc
    qrad = qr(t)
    ql = qlt(Tcore[0],f,dTl_dP)
    #qg = qgt(Tcore[0],f,dTl_dP,Xs) #exclude as makes neglible difference
    qg=0
    
    dTcdt = (qrad-Fcmb*Acmb)/(qst-ql+qg)
    dfdt = - dTcdt/(rhoc*gc*dTl_dP*rc)
    f_new = f+dfdt*dt
    Rem, Bdip_cmb, comp, therm = rem_b(f, dfdt, Xs, Tcore, Fcmb, True,0) #if core is solidifying there is no thermal stratification
    #calculate heat fluxes to return
    qcore =np.array([qrad,qst*dTcdt,ql*dTcdt,qg*dTcdt])
    
    return dTcdt, f_new, Rem, Bdip_cmb, comp, therm, qcore