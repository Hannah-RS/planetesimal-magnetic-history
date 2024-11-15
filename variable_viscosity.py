#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from parameters import alpha_n, Tml, Tms, rcmf, eta0, etal, w, beta 

def part_a(T,Tms):
    """
    Exponential (Frank-Kamenetskii) viscosity variation with no melt fraction component
    Eqn. 1a from Sanderson et. al. (2024)
    
    Parameters
    ----------
    T : float
        temperature [K]
    Tms : float
        solidus temperature [K]
    Returns
    -------
    eta : float
        viscosity [Pas]

    """
    eta = eta0*np.exp(-beta*(T-Tms))
    return eta

def part_b(T,Tms,Tml):
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
    Returns
    -------
    eta : float
        viscosity [Pas]

    """
    eta = eta0*np.exp(-(beta+alpha_n/(Tml-Tms))*(T-Tms))
    return eta

def part_c(T,Tms,Tml,Trcmf):
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
    Returns
    -------
    eta : float
        viscosity [Pas]

    """
    etaa = np.log10(part_b(Trcmf,Tms,Tml))
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

def eta_calc(T,Tms=Tms,Tml=Tml):
    """
    Piecewise mantle viscosity 
    Eqn 1. from Sanderson et. al. (2024)
    Parameters
    ----------
    T : float
        temperature [K]
    Tms : float
        solidus temperature [K], default in parameters file
    Tml : float
        liquidus temperature [K], default in parameters file
    Returns
    -------
    eta : float
        viscosity [Pas]
    """
    Trcmf = Tms + rcmf*(Tml-Tms) #calculate Trcmf
    if type(T) == np.ndarray:
        phi = (T-Tms)/(Tml-Tms)
        eta1 = part_a(T[phi<=0],Tms)
        eta2 = part_b(T[((phi>0)&(phi<rcmf))],Tms,Tml)
        eta3 = part_c(T[(phi>=rcmf)&(T<(Trcmf+w))],Tms,Tml,Trcmf)
        eta4 = part_d(T[T>=(Trcmf+w)],Tms,Tml)
        eta = np.concatenate((eta1,eta2,eta3,eta4),dtype='float64')
    else:
        phi = (T-Tms)/(Tml-Tms)
        if phi<=0:
            eta = part_a(T,Tms)
        elif (phi>0) & (phi <= rcmf):
            eta = part_b(T,Tms,Tml)
        elif (phi>rcmf) & (T<(Trcmf+w)):
            eta = part_c(T,Tms,Tml,Trcmf)
        elif T>=(Trcmf+w):
            eta= part_d(T,Tms,Tml)  
        else:
            print(T,'Scenario not considered')
    return eta   

