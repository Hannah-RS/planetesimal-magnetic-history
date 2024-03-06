#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot viscosity model for chosen set of parameters (in parameters file)
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
# setting path
sys.path.append('../')
from viscosity_def import viscosity
from parameters import Tml, Tms
import pandas as pd

# create an array of mantle temperatures
n = 100 #number of temperature points
Tm = np.linspace(1200,1700, n)
phi = (Tm-Tms)/(Tml-Tms)
#calculate viscosity
eta = viscosity(Tm,'vary')

#%% Make plot
fig, ax1 = plt.subplots(1,1)
ax2 = ax1.twiny()
ax1.semilogy(Tm, eta,color='black')
ax2.semilogy(phi,eta, color='black')
ax2.set_xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7],[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7])
ax1.set_ylabel('Viscosity/Pas')
ax1.set_xlabel('T$_m$/K')
ax2.set_xlabel('Melt fraction, $\phi$',)
plt.savefig('../Plots/Icarus_paper/viscosity_profile.pdf',bbox_inches='tight',dpi=500)