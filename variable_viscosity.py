#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from parameters import alpha_n, Tml, Tms, rcmf, eta0, etal, w, beta, xwater, Dw, mmr 

def wtperc_to_h(xwater,mr):
    """
    Convert water concentrations in wt % to H/1e6 Si
    Parameters:
        xwater: float
            Concentration of water in wt %
        mr: float
            Molar mass of water bearing mineral e.g. forsterite = 140 g/mol
    Returns:
        ch: float
            Concentration of H in 1e6 Si

    """
    mwater = 18.01528 #molar mass of water
    ch = xwater*1000*2*mr/mwater
    return ch

def ch_sol(ch_tot,phi,D):
    """
    Water concentration in H/1e6 Si in solid phase
    Parameters:
        ch_tot: float
            Total water concentration in H/1e6 Si
        phi: float
            Melt fraction
        D: float
            Partition coefficient
    Returns:
        ch_s: float
            Water concentration in H/1e6 Si in solid phase
    """
    ch_s = ch_tot/((1-phi)+phi/D)
    if type(phi) != np.float64: #dealing with an array, ignore negative values
        ch_s[phi<0] = ch_tot
    return ch_s

def part_a(T,Tms,ch_tot):
    """
    Exponential (Frank-Kamenetskii) viscosity variation with no melt fraction component
    Eqn. 1a from Sanderson et. al. (2024)
    
    Parameters
    ----------
    T : float
        temperature [K]
    Tms : float
        solidus temperature [K]
    ch_tot: float
        total water concentration [H / 10^6 Si]
    Returns
    -------
    eta : float
        viscosity [Pas]

    """
    if ch_tot > 0: #subsolidus - no melt partitioning
        eta = (10/ch_tot)*eta0*np.exp(-beta*(T-Tms)) # 1 or 10
    else:
        eta = eta0*np.exp(-beta*(T-Tms))
    return eta

def part_b(T,Tms,Tml,ch_tot):
    """
    Exponential (Frank-Kamenetskii) viscosity variation with melt weakening
    Eqn. 1b from Sanderson et. al. (2024)
    
    Parameters
    ----------
    T : float
        temperature [K]
    Tms : float
        solidus temperature [K]
    Tml : float
        liquidus temperature [K]
    ch_tot: float
        total water concentration [H / 10^6 Si]
    Returns
    -------
    eta : float
        viscosity [Pas]

    """
    if ch_tot > 0: #calculate melt partitioning
        phi = (T-Tms)/(Tml-Tms)
        ch_s = ch_sol(ch_tot, phi, Dw)
        eta = (10/ch_s)*eta0*np.exp(-(beta+alpha_n/(Tml-Tms))*(T-Tms)) #1 or 10?
    else:
        eta = eta0*np.exp(-(beta+alpha_n/(Tml-Tms))*(T-Tms))
    return eta

def part_c(T,Tms,Tml,Trcmf,ch_tot):
    """
    Linear in log-space approximation to initial decrease after critical melt fraction
    Eqn. 1c from Sanderson et. al. (2024)

    Parameters
    ----------
    T : float
        temperature [K]
    Tms : float
        solidus temperature [K]
    Tml : float
        liquidus temperature [K]
    T: float
        temperature of critical melt fraction [K]
    ch_tot: float
        total water concentration [H / 10^6 Si]
    Returns
    -------
    eta : float
        viscosity [Pas]
    """
    etaa = np.log10(part_b(Trcmf,Tms,Tml,ch_tot))
    etab = np.log10(part_d(Trcmf+w,Tms,Tml))
    log_eta = etaa + (etab-etaa)/w*(T-Trcmf)
    eta = 10**log_eta
    return eta

def part_d(T,Tms,Tml):
    """
    Krieger-Dougherty viscosity variation for supercritical melt fractions
    Eqn. 1d from Sanderson et. al. (2024)
    
    Parameters
    ----------
    T : float
        temperature [K]
    Tms : float
        solidus temperature [K]
    Tml : float
        liquidus temperature [K]
    Returns
    -------
    eta : float
        viscosity [Pas]

    """
    phi = (T-Tms)/(Tml-Tms)
    #Krieger-Dougherty
    if np.any((phi-rcmf) <0):
        print(phi[(phi-rcmf)<0])
    eta = etal*((phi-rcmf)/(1-rcmf))**(-2.5*(1-rcmf))
    return eta

def eta_calc(T,Tms=Tms,Tml=Tml,xwater=xwater):
    """
    Piecewise mantle viscosity 
    Eqn 1. from Sanderson et. al. (2024)
    Modified to include water content Sanderson et. al. (2025)
    Parameters
    ----------
    T : float
        temperature [K]
    Tms : float
        solidus temperature [K], default in parameters file
    Tml : float
        liquidus temperature [K], default in parameters file
    xwater : float
        water content [wt%], default in parameters file
    Returns
    -------
    eta : float
        viscosity [Pas]
    """
    Trcmf = Tms + rcmf*(Tml-Tms) #calculate Trcmf
    if xwater > 0:
        ch_tot = wtperc_to_h(xwater,mmr)   
    else:
        ch_tot = 0
    if type(T) == np.ndarray:
        phi = (T-Tms)/(Tml-Tms)
        eta1 = part_a(T[phi<=0],Tms,ch_tot)
        eta2 = part_b(T[((phi>0)&(phi<rcmf))],Tms,Tml,ch_tot)
        eta3 = part_c(T[(phi>=rcmf)&(T<(Trcmf+w))],Tms,Tml,Trcmf,ch_tot)
        eta4 = part_d(T[T>=(Trcmf+w)],Tms,Tml)
        eta = np.concatenate((eta1,eta2,eta3,eta4),dtype='float64')
       
    else:
        phi = (T-Tms)/(Tml-Tms)
        if phi<=0:
            eta = part_a(T,Tms,ch_tot)
        elif (phi>0) & (phi <= rcmf):
            eta = part_b(T,Tms,Tml,ch_tot)
        elif (phi>rcmf) & (T<(Trcmf+w)):
            eta = part_c(T,Tms,Tml,Trcmf,ch_tot)
        elif T>=(Trcmf+w):
            eta= part_d(T,Tms,Tml)  
        else:
            print(T,'Scenario not considered')
    return eta

