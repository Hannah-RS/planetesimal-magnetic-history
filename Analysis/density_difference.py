#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot fractional density between bulk liquid core 
and solidified iron as a function of S content
"""
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../')
from parameters import Mr_fe, Mr_s, rho_eut, rhofe_l, rhofe_s, Tl_fe,\
    alpha_c, rho_exp, Xs_eutectic
from fe_fes_liquidus import fe_fes_density, weight_perc_to_at_frac

Xs = np.linspace(0.1,Xs_eutectic,100)
Xsd = Xs/100
at = weight_perc_to_at_frac(Xs)
Mrr = 1+Mr_fe/Mr_s

rhom19 = fe_fes_density(Xs)*rho_exp

drho = rhofe_s-rhom19


plt.figure()
plt.plot(Xs,drho,label='$\\Delta \\rho =\\rho_{Fe,solid}-\\rho_c$',color='black')
plt.plot(Xs,0.05*rhom19,label='$\\Delta \\rho =$const.',color='black',linestyle='dashed')
plt.xlabel('Core sulfur content /wt%')
plt.ylabel('$\\Delta \\rho $ / kg $m^{-3}$ ')
plt.xlim([0,33])
plt.legend()
plt.savefig('../Plots/density_difference.png',dpi=500)

#%% Absolute core densities 
# plt.figure()
# plt.plot(Xs,rhom19*rho_exp,label='Morard 2019')
# plt.scatter(min(Xs),7800,label='solid Fe - Bryson 2015',color='black',marker='x')
# plt.scatter(32,4992,label='eutectic FeS - Dodds 2022',color='green',marker='x')
# plt.scatter(min(Xs),6980,label='liquid Fe - Dodds 2022',color='blue',marker='x')
# plt.legend(loc='lower left')
# plt.xlabel('Xs/ wt %')
# plt.ylabel('$\\rho / kg m^{-3}$')
# #plt.savefig('Plots/core_density.png')

# plt.figure()
# plt.plot(at,rhom19,label='Morard 2019')
# plt.hlines(7800,min(at),max(at),label='solid Fe - Bryson 2015',linestyle='dashed',color='black')
# plt.hlines(4992,max(at)-0.05,max(at),label='eutectic FeS - Dodds 2022',linestyle='dashed',color='green')
# plt.hlines(6980,min(at),max(at),label='liquid Fe - Dodds 2022',linestyle='dashed',color='red')
# plt.legend(loc='lower left')
# plt.xlabel('atom % S')
# plt.ylabel('$\\rho kg m^{-3}$')