#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 11:10:18 2023

@author: exet5460
"""
import matplotlib.pyplot as plt
import numpy as np
from parameters import Mr_fe, Mr_s, rho_eut, rhofe_l, rhofe_s

Xs = np.linspace(25,32,100)
Xsd = Xs/100
Mrr = 1+Mr_fe/Mr_s
rhol = Xsd*Mrr*rho_eut + (1-Xsd*Mrr)*rhofe_l
drho1 = (rhofe_s - rhol)/rhol
drho2 = (rhofe_s - rho_eut)/rhol

plt.figure()
plt.plot(Xs,drho1,label='$\\Delta \\rho = \\rho_{Fe,sol}-\\rho_c$')
plt.plot(Xs,drho2,label='$\\Delta \\rho =\\rho_{Fe,sol}-\\rho_{FeS,eut}$')
plt.xlabel('Xs/wt%')
plt.ylabel('$\\frac{\\Delta \\rho}{\\rho_c} $ ')
plt.legend()
plt.savefig('density_difference.png')
