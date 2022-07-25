#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for comparing different viscosity profiles
"""
import numpy as np
import matplotlib.pyplot as plt
from viscosity_def import viscosity

# create an array of mantle temperatures
n = 100 #number of temperature points
Tm = np.linspace(1400,1700, n)

eta_Dodds = viscosity(Tm, 'Dodds')
eta_Bryson = viscosity(Tm, 'Bryson')
#eta_Arr = viscosity(Tm, 'Arrhenius')
eta_Dodds2 = viscosity(Tm, 'new')

#plot models
plt.figure()
plt.semilogy(Tm, eta_Dodds, label='Dodds')
plt.semilogy(Tm, eta_Bryson, label='Bryson')
#plt.semilogy(Tm, eta_Arr, label='Arrhenius', linestyle = '--', color ='red')
plt.semilogy(Tm, eta_Dodds2, label='new')
plt.xlabel('Tm/K')
plt.ylabel('Viscosity/Pas')
plt.legend()
plt.savefig('Plots/viscosity_investigation.png')