#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 10:25:30 2022

@author: exet5460
"""
from flux_definitions_dep import f1_conductive, f1_convective, f2_conductive, f2_convective, flux_balance
import numpy as np
import matplotlib.pyplot as plt

Tm = 1601.4546508457229
Tc =1599.9999999999995

Tcmb = np.linspace(Tc+1e-5,Tm-1e-5,100)

f1_cond = f1_conductive(Tm,Tc,Tcmb)
f1_conv = f1_convective(Tm,Tc,Tcmb)
f2_cond = f2_conductive(Tm,Tc,Tcmb)
f2_conv = f2_convective(Tm,Tc,Tcmb)

plt.figure()
#plt.semilogy(Tcmb,abs(f1_cond),label='f1 cond')
plt.semilogy(Tcmb,abs(f1_conv),label='f1 conv')
#plt.semilogy(Tcmb,abs(f2_cond),label='f2 cond')
plt.semilogy(Tcmb,abs(f2_conv),label='f2 conv')
plt.legend()
plt.xlabel('CMB temperature/K')
plt.ylabel('Absolute value of heat flux')
#plt.savefig('Plots/Heat_flux_convergence.png')

#compare intersection with solver
import scipy.optimize as sco
Tcmb_actual = sco.root_scalar(flux_balance,args=(Tc,Tm,f1_convective,f2_convective),x0=Tc+1e-5,x1=Tm-1e-5).root
