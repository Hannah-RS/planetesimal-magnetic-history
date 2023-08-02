#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plots for creating variable viscosity profiles
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../')
from parameters import alpha_n, eta0, Tml, Tms, rcmf

def part_a(T,eta0,g):
    """
    Exponential (Frank-Kamenetskii) viscosity variation with no melt fraction component
    
    Parameters
    ----------
    T : float
        temperature [K]
    eta0 : float
        reference viscosity [Pas]
    g : float
        variation of viscosity with temperature [K^-1]

    Returns
    -------
    eta : float
        viscosity [Pas]

    """
    eta = eta0*np.exp(-g*(T-Tms))
    return eta

def part_b(T,eta0,g):
    """
    Exponential (Frank-Kamenetskii) viscosity variation
    
    Parameters
    ----------
    T : float
        temperature [K]
    eta0 : float
        reference viscosity [Pas]
    g : float
        variation of viscosity with temperature [K^-1]

    Returns
    -------
    eta : float
        viscosity [Pas]

    """
    eta = eta0*np.exp(-(g+alpha_n/(Tml-Tms))*(T-Tms))
    return eta

def part_c(T,Trcmf,etaa,etab):
    """
    Linear in log-space approximation to initial decrease after rcmf

    Parameters
    ----------
    T : float
        temperature [K]
    Trcmf : float
        temperature of the rheologically critical melt fraction [K]
    etaa : float
        log10 of viscosity at Trcmf  (defined by part_b)
    etab : float
        log10 of viscosity at Trcmf+w  (defined by part_d)

    Returns
    -------
    eta : float
        viscosity [Pas]

    """
    log_eta = etaa + (etab-etaa)/w*(T-Trcmf)
    eta = 10**log_eta
    return eta

def part_d(T,rcmf,eta0,g):
    """
    Krieger-Dougherty viscosity variation for supercritical melt fractions

    Parameters
    ----------
    T : float
        temperature [K]
    rcmf : float
        rheologically critical melt fraction
    eta0 : float
        reference viscosity [Pas]
    g : float
        variation of viscosity with temperature [K^-1]

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

def eta_calc(T,rcmf,Trcmf,eta0,g):
    """
    Piecewise viscosity calculation
    Parameters
    ----------
    T : float
        temperature [K]
    rcmf : float
        rheologically critical melt fraction
    Trcmf : float
        temperature of the rheologically critical melt fraction [K]
    eta0 : float
        reference viscosity [Pas]
    g : float
        variation of500 viscosity with temperature [K^-1]

    Returns
    -------
    eta : float
        viscosity [Pas]
    """
    if type(T) == np.ndarray:
        phi = (T-Tms)/(Tml-Tms)
        eta1 = part_a(T[phi<=0],eta0,g)
        eta2 = part_b(T[((phi>0)&(phi<=rcmf))],eta0,g)
        eta3 = part_c(T[(phi>rcmf)&(T<(Trcmf+w))],Trcmf,etaa,etab)
        eta4 = part_d(T[T>=(Trcmf+w)],rcmf,eta0,g)
        eta = np.concatenate((eta1,eta2,eta3,eta4))
    else:
        phi = (T-Tms)/(Tml-Tms)
        if phi<0:
            eta = part_a(T,eta0,g)
        elif (phi>0) & (phi <= rcmf):
            eta = part_b(T,eta0,g)
        elif (phi>rcmf) & (T<(Trcmf+w)):
            eta = part_c(T,Trcmf,etaa,etab)
        elif T>=(Trcmf+w):
            eta= part_d(T,rcmf,eta0,g)  
        else:
            print(T,'Scenario not considered')
    return eta

#things that will be determined in the parameters file
w = 10 #width of linear decrease region [K]
etal = 100 # liquid viscosity [Pas]
Trcmf = rcmf*(Tml-Tms)+Tms
from parameters import E, R, Tref
g = E/(R*Tref**2)
eta0 = 1e21

#calculate at the top of the script on import so can be cached
etaa = np.log10(part_b(Trcmf,eta0,g))
etab = np.log10(part_d(Trcmf+w,rcmf,eta0,g))
if etab > etaa:
    raise ValueError(f'rcmf = {rcmf}, eta0 ={eta0}, g={g} is an invalid viscosity model - exponential decrease is too steep.')
#testing
T = np.linspace(1200,1800,200)

# phi = (T-Tms)/(Tml-Tms)
# eta = eta_calc(1600,rcmf,eta0,g)
# eta = eta_calc(T,rcmf,eta0,g)
# #eta = part_b(T,rcmf,eta0,g)
# plt.figure()
# plt.scatter(T,eta)
# plt.xlabel('T/K')if etab > etaa:
# plt.ylabel('$\\eta$ \Pas')
# plt.yscale('log')

#for parameter list
n=10
m=10
p = 10
rcmf = np.linspace(0.2,0.5,n)
Trcmf = rcmf*(Tml-Tms)+Tms
eta0 = np.logspace(14,21,m)
g = np.linspace(0.005,0.08,p)

mono = np.zeros([n,m,p])
#calculate viscosity
for i in range(n):
    for j in range(m):
        for k in range(p):
            etaa = np.log10(part_b(Trcmf[i],eta0[j],g[k]))
            etab = np.log10(part_d(Trcmf[i]+w,rcmf[i],eta0[j],g[k]))
            if etab > etaa:
                raise ValueError(f'rcmf = {rcmf[i]:.2g}, eta0 ={eta0[j]:.2e}, g={g[k]:.3f} is an invalid viscosity model - exponential decrease is too steep.')
            eta = eta_calc(T,rcmf[i],Trcmf[i],eta0[j],g[k])
            eta_diff = np.diff(eta) #calculate differences with sucessive elements
            # plt.figure()
            # plt.plot(T,eta)
            # plt.xlabel('T/K')
            # plt.ylabel('$\\eta$ \Pas')
            # plt.yscale('log')
            # plt.title(f'rcmf={rcmf[i]}, eta0={eta0[j]:.1e}, g={g[k]}')
            if np.all(eta_diff<=0):
                mono[i,j,k] = 1 #function is monotonically decreasing yay!
    
if np.all(mono==1):
    print('All the versions are monotonic - yay!')
else:
    print('Sad times')
    
eta = eta_calc(T,rcmf[6],Trcmf[6],eta0[0],g[8])
plt.figure()
plt.scatter(T,eta)
plt.xlabel('T/K')
plt.ylabel('$\\eta$ \Pas')
plt.yscale('log')