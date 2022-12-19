#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for comparing different viscosity profiles
"""
import numpy as np
import matplotlib.pyplot as plt
from viscosity_def import viscosity
from parameters import Tml, Tms
import pandas as pd

# create an array of mantle temperatures
n = 100 #number of temperature points
Tm = np.linspace(1400,1700, n)
phi = (Tm-Tms)/(Tml-Tms)

# import scott and kohlstedt data
sk_data=pd.read_excel('Results/scott_and_kohlstedt_06.xlsx')
sk_diff = sk_data[sk_data['regime']=='diffusion'] # get diffusion points
#hk_diff = sk_diff[sk_diff['source']=='hk1'] # get hirth and kohlstedt points (3)
#sk_diff = sk_diff[sk_diff['source']=='sk'] #get scott and kohlstedt points
sk_gbs = sk_data[sk_data['regime']=='GBS'] #get gbs points
#hk_gbs = sk_gbs[sk_gbs['source']=='hk2'] # get hirth and kohlstedt (4)
#sk_gbs = sk_gbs[sk_gbs['source']=='sk'] # get scott and kohlstedt points

eta_Dodds = viscosity(Tm, 'Dodds')
eta_Bryson = viscosity(Tm, 'Bryson')
eta_Steren = viscosity(Tm, 'Sterenborg')
#eta_Arr = viscosity(Tm, 'Arrhenius')
#eta_Dodds2 = viscosity(Tm, 'new')

#plot models
plt.figure()
plt.semilogy(Tm, eta_Dodds, label='Dodds et al. (2021)',color='cornflowerblue')
plt.semilogy(Tm, eta_Bryson, label='Bryson et al. (2019)',linestyle='--',color='green')
plt.semilogy(Tm, eta_Steren, label='Sterenborg and Crowley (2013)',linestyle='dotted',color='navy')
#plt.semilogy(Tm, eta_Arr, label='Arrhenius', linestyle = '--', color ='red')
#plt.semilogy(Tm, eta_Dodds2, label='new')
plt.xlabel('Tm/K')
plt.ylabel('Viscosity/Pas')
plt.legend()
#plt.savefig('Plots/viscosity_comparison_T.png')

#melt fraction version
#plot models
plt.figure()
plt.semilogy(phi, eta_Dodds, label='Dodds et al. (2021)',color='cornflowerblue')
plt.semilogy(phi, eta_Bryson, label='Bryson et al. (2019)',linestyle='--',color='green')
plt.semilogy(phi, eta_Steren, label='Sterenborg and Crowley (2013)',linestyle='dotted',color='navy')
#plt.semilogy(phi, eta_Arr, label='Arrhenius', linestyle = '--', color ='red')
#plt.semilogy(phi, eta_Dodds2, label='new')
#plt.xlabel('Melt fraction')

plt.ylabel('Viscosity/Pas')
plt.legend()
#plt.savefig('Plots/viscosity_comparison_phi.png')

########################  plot models with two axes  ###########################
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
ax1.semilogy(phi, eta_Steren, label='Sterenborg and Crowley (2013)',linestyle='dotted',color='navy')
ax2.semilogy(Tm, eta_Dodds, label='Dodds et al. (2021)',color='cornflowerblue')
ax2.semilogy(Tm, eta_Bryson, label='Bryson et al. (2019)',linestyle='--',color='green')

ax1.set_ylabel('Viscosity/Pas')
ax2.set_xlabel('T$_m$/K')
ax1.set_xlabel('Melt fraction, $\phi$')
fig.legend(loc=(0.15,0.15))
plt.savefig('Plots/viscosity_comparison.png')

######################### Version with all experiments ########################
plt.figure()
plt.semilogy(phi, eta_Dodds,color='cornflowerblue')
plt.semilogy(phi, eta_Bryson,linestyle='--',color='green')
plt.semilogy(phi, eta_Steren,linestyle='dotted',color='navy')
plt.scatter(sk_diff['melt fraction'],sk_diff['viscosity (Pas)'],marker = 'x',color = 'black',label='diffusion creep - Scott & Kohlstedt (2006)')
#plt.scatter(hk_diff['melt fraction'],hk_diff['viscosity (Pas)'],marker = '+',color = 'black')
plt.scatter(sk_gbs['melt fraction'],sk_gbs['viscosity (Pas)'],marker = 'x',color = 'red',label='dislocation creep - Scott & Kohlstedt (2006)')
#plt.scatter(hk_gbs['melt fraction'],hk_gbs['viscosity (Pas)'],marker = '+',color = 'red')
plt.arrow(x=0.3, y=1e4, dx = 0, dy = 1e7,head_width=0.02,head_length=5e7,color='black')#,text='upper bound') #start here
plt.xlabel('Melt fraction')
plt.ylabel('Viscosity/Pas')
plt.legend()
plt.savefig('Plots/experiments.pdf')

################# plot Bryson model sequentially  #############################
plt.figure()
plt.semilogy(phi[phi<0.5],eta_Bryson[phi<0.5],linestyle='--')
plt.xlabel('Melt fraction')
plt.ylabel('Viscosity/Pas')
plt.ylim([1,1e22])
plt.xlim([0,0.8])
#plt.savefig('Plots/Bryson_1.png')

plt.figure()
plt.semilogy(phi[phi<0.63],eta_Bryson[phi<0.63],linestyle='--')
plt.semilogy(phi[phi<0.5],eta_Bryson[phi<0.5],color='grey')
plt.xlabel('Melt fraction')
plt.ylabel('Viscosity/Pas')
plt.ylim([1,1e22])
plt.xlim([0,0.8])
#plt.savefig('Plots/Bryson_2.png')

plt.figure()
plt.semilogy(phi,eta_Bryson,linestyle='--')
plt.semilogy(phi[phi<0.63],eta_Bryson[phi<0.63],color='grey')
plt.xlabel('Melt fraction')
plt.ylabel('Viscosity/Pas')
plt.ylim([1,1e22])
plt.xlim([0,0.8])
#plt.savefig('Plots/Bryson_3.png')
