#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for comparing different viscosity profiles
"""
import numpy as np
import matplotlib.pyplot as plt
from viscosity_def import viscosity
from parameters import Tml, Tms

# create an array of mantle temperatures
n = 100 #number of temperature points
Tm = np.linspace(1400,1700, n)
phi = (Tm-Tms)/(Tml-Tms)


eta_Dodds = viscosity(Tm, 'Dodds')
eta_Bryson = viscosity(Tm, 'Bryson')
eta_Steren = viscosity(Tm, 'Sterenborg')
#eta_Arr = viscosity(Tm, 'Arrhenius')
#eta_Dodds2 = viscosity(Tm, 'new')

#plot models
plt.figure()
plt.semilogy(Tm, eta_Dodds, label='Dodds et al.')
plt.semilogy(Tm, eta_Bryson, label='Bryson et al. - FK')
plt.semilogy(Tm, eta_Steren, label='Sterenborg and Crowley')
#plt.semilogy(Tm, eta_Arr, label='Arrhenius', linestyle = '--', color ='red')
#plt.semilogy(Tm, eta_Dodds2, label='new')
plt.xlabel('Tm/K')
plt.ylabel('Viscosity/Pas')
plt.legend()
#plt.savefig('Plots/viscosity_comparison_T.png')

#melt fraction version
#plot models
plt.figure()
plt.semilogy(phi, eta_Dodds, label='Dodds et al.')
plt.semilogy(phi, eta_Bryson, label='Bryson et al.')
plt.semilogy(phi, eta_Steren, label='Sterenborg and Crowley')
#plt.semilogy(phi, eta_Arr, label='Arrhenius', linestyle = '--', color ='red')
#plt.semilogy(phi, eta_Dodds2, label='new')
plt.xlabel('Melt fraction')
plt.ylabel('Viscosity/Pas')
plt.legend()
plt.savefig('Plots/viscosity_comparison_phi.png')