#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Comparison of viscosity profile from Bryson et al. (2019) using FK approximation with Arrhenius expression 
from Robuchon & Nimmo (2011)
The constants used are the same in both expressions, this script just explores the functional form
"""
import numpy as np
import matplotlib.pyplot as plt
from viscosity_def import viscosity
from parameters import Tml, Tms

# create an array of mantle temperatures
n = 100 #number of temperature points
Tm = np.linspace(1400,1600, n)
phi = (Tm-Tms)/(Tml-Tms)


eta_Bryson = viscosity(Tm, 'Bryson')
eta_Robuchon = viscosity(Tm, 'Robuchon-Bryson')
eta_Dodds = viscosity(Tm, 'Dodds')


#plot models
plt.figure()
plt.semilogy(Tm, eta_Bryson, label='Bryson et al. - FK')
plt.semilogy(Tm, eta_Robuchon, label='Bryson et al. - Arrhenius')
plt.semilogy(Tm, eta_Dodds, label='Dodds et al.')
plt.xlabel('Tm/K')
plt.ylabel('Viscosity/Pas')
plt.legend()
plt.savefig('Plots/viscosity_comparison_FK.png')

#melt fraction version
#plot models
plt.figure()
plt.semilogy(phi, eta_Bryson, label='Bryson et al. - FK')
plt.semilogy(phi, eta_Robuchon, label='Bryson et al. - Arrhenius')
plt.semilogy(phi, eta_Dodds, label='Dodds et al.')
plt.xlabel('Melt fraction')
plt.ylabel('Viscosity/Pas')
plt.legend()
plt.savefig('Plots/viscosity_comparison_FKphi.png')