#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check volume conservation of core is valid during solidification
"""
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('../')

from fe_fes_liquidus import fe_fes_density
from parameters import rho_exp, rhofe_s 
rc = 200e3
f0 = 0.999
Xs0_val = [10, 15, 20, 27.1] #lowest possible initial conc

#make plot for different initial sulfur contents
#if assume constant volume how does mass change
plt.figure()
for Xs_0 in Xs0_val:
    ffinal = (Xs_0/32)**(1/3)
    m0 = 4/3*np.pi*rc**3*fe_fes_density(Xs_0)*rho_exp

    f = np.linspace(f0,ffinal,100)
    Xs = f**(-3)*Xs_0

    m = 4/3*np.pi*(f*rc)**3*fe_fes_density(Xs)*rho_exp+4/3*np.pi*(rc**3-(f*rc)**3)*rhofe_s

    plt.plot(f,(m-m0)/m0,label=f'{Xs_0} wt%')

plt.xlabel('$r_i$/$r_c$ up to eutectic')
plt.ylabel('$\Delta M$/M$_{initial}$')
plt.legend()
plt.xlim([f0,0.6])

#if assume constant mass how does volume change
plt.figure()
for Xs_0 in Xs0_val:
    ffinal = (Xs_0/32)**(1/3)
    m0 = 4/3*np.pi*rc**3*fe_fes_density(Xs_0)*rho_exp

    f = np.linspace(f0,ffinal,100)
    Xs = f**(-3)*Xs_0

    rc2 = ((f*rc)**3*(1-fe_fes_density(Xs)*rho_exp/rhofe_s)+rc**3*fe_fes_density(Xs_0)/rhofe_s)**(1/3)

    plt.plot(f,(rc2-rc)/rc,label=f'{Xs_0} wt%')

plt.xlabel('$r_i$/r$_{c,0}$ up to eutectic')
plt.ylabel('$\Delta r$/r$_{c,0}$')
plt.legend()
plt.xlim([f0,0.6])
