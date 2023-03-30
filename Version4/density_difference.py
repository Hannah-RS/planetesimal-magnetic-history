#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 11:10:18 2023

@author: exet5460
"""
import matplotlib.pyplot as plt
import numpy as np
from parameters import Mr_fe, Mr_s, rho_eut, rhofe_l, rhofe_s

Xs = np.linspace(0,32,100)
Xsd = Xs/100
Mrr = 1+Mr_fe/Mr_s
rhol = Xsd*Mrr*rho_eut + (1-Xsd*Mrr)*rhofe_l
drho1 = (rhofe_s - rhol)/rhol
drho2 = (rhofe_s - rho_eut)/rhol

ticks = np.arange(0,0.65,0.05)
plt.figure()
plt.plot(Xs,drho1,label='$\\Delta \\rho = \\rho_{Fe,solid}-\\rho_c$')
plt.plot(Xs,drho2,label='$\\Delta \\rho =\\rho_{Fe,solid}-\\rho_{FeS,eut}$')
plt.hlines(0.05,0,32,linestyle='dashed',color='black',label='Nimmo (2009)')
plt.xlabel('Xs/wt%')
plt.ylabel('$\\frac{\\Delta \\rho}{\\rho_c} $ ')
plt.yticks(ticks,minor=False)
plt.legend(bbox_to_anchor=(0.97,0.4))
plt.savefig('Plots/density_difference.png')
