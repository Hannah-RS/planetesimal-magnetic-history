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

etal = 100 # liquid viscosity [Pas]

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

def part_c(T,rcmf,eta0,g):
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
    eta = etal*((phi-rcmf)/(1-rcmf))**(-2.5*(1-rcmf))
    return eta

def eta_calc(T,rcmf,eta0,g):
    """
    Piecewise viscosity calculation
    Parameters
    ----------
    T : float
        temperature [K]
    rcmf : float
        rheologically critical melt fraction
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
        eta2 = part_b(T[(phi>0)&(phi<=rcmf)],eta0,g)
        eta3 = part_c(T[phi>rcmf],rcmf,eta0,g)
        eta = np.append(eta1,eta2)
        eta = np.append(eta,eta3)
    else:
        phi = (T-Tms)/(Tml-Tms)
        if phi<0:
            eta = part_a(T,eta0,g)
        elif (phi>0) & (phi <= rcmf):
            eta = part_b(T,eta0,g)
        else:
            eta= part_c(T,rcmf,eta0,g)       
    return eta
    
from parameters import E, R, Tref
#eta0 = 1e14
T = np.linspace(1200,1800,100)
g = E/(R*Tref**2)

eta = eta_calc(1600,rcmf,eta0,g)
eta = eta_calc(T,rcmf,eta0,g)
#eta = part_b(T,rcmf,eta0,g)
plt.figure()
plt.plot(T,eta)
plt.xlabel('T/K')
plt.ylabel('$\\eta$ \Pas')
plt.yscale('log')

#for parameter list
n=6
m=6
p = 6
rcmf = np.linspace(0.2,0.5,n)
eta0 = np.logspace(14,21,m)
g = np.linspace(0.005,0.05,p)

mono = np.zeros([n,m,p])
#calculate viscosity
for i in range(n):
    for j in range(m):
        for k in range(p):
            eta = eta_calc(T,rcmf[i],eta0[j],g[k])
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