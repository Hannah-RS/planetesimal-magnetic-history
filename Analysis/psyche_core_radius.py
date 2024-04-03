#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate possible core-mantle ratios for different core sulfur contents for Psyche
Assume constant bulk density
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../')

from fe_fes_liquidus import fe_fes_density

#%% Define functions
def weight_perc_fe(Xs0,rhoc,rhom,rhoa):
    """

    Parameters
    ----------
    Xs0 : float
        initial core sulfur content [wt %]
    rhoc : float
        core density [kg m^-3]
    rhom : float
        mantle density [kg m^-3]
    rhoa : float
        bulk density [kg m^-3]
        
    Returns
    -------
    XFe : float
        weight percent iron in accreted material

    """
    XFe = (rhoc/rhoa)*((rhoa-rhom)/(rhoc-rhom))*(1-0.01*Xs0)
    return XFe

def core_radius(Xs0,rhom,rhoa):
    """

    Parameters
    ----------
    Xs0 : float
        initial core sulfur content [wt %]
    rhom : float
        mantle density [kg m^-3]
    rhoa : float
        bulk density [kg m^-3]

    Returns
    -------
    rc : float
        core radius as fraction of planetesimal radius

    """
    rhoc = fe_fes_density(Xs0)
    XFe = weight_perc_fe(Xs0,rhoc,rhom,rhoa)
    rc = ((XFe/(1-0.01*Xs0))*(rhoa/rhoc))**(1/3)
    return rc

#%% Set parameters
rhoms = np.array([2500,3000,3500]) # mantle density [kg m^-3] Johnson et. al. (2020)
rhoa = 4000 # bulk density [kg m^-3]

#%% Create arrays for output
Xs0 = np.linspace(0,33,33)

#%% Plot result for rc
colors = ['navy','mediumblue','skyblue']
lines = ['-','-.',':']
fig, ax = plt.subplots(1,1)
ax.grid(which='both',alpha=0.5)
for rhom, color, ls in zip(rhoms,colors,lines):
    rc = core_radius(Xs0,rhom,rhoa)
    ax.plot(Xs0,rc,label=f'$\\rho_m =${rhom} $kg m^{{-3}}$',color=color,linestyle=ls)
ax.set(ylabel='Core radius fraction - $\\frac{r_c}{r}$',xlabel='Initial core sulfur content/ wt %')
ax.legend()
ax.set_title('bulk density 4000 kg$m^{-3}$')
fig.savefig('../Plots/Psyche/core_radius_frac.png',dpi=500)
#%% Plot results for XFe
rhoc = fe_fes_density(Xs0)
Xfe = weight_perc_fe(Xs0,rhoc,rhoms[1],rhoa)
maxfe = weight_perc_fe(fe_fes_density(27.1),27.1)
fig, ax = plt.subplots(1,1)
ax.plot(Xs0,Xfe)
ax.set(ylabel='Iron abundance/ wt %',xlabel='Initial core sulfur content/ wt %')
ax.hlines(maxfe,0,27.1,linestyle='dashed',color='black')
ax.vlines(27.1,0.45,maxfe,linestyle='dashed',color='black')
ax.text(0,0.565,f'{maxfe:.2f}')
ax.text(27.2,0.45,'27.1')