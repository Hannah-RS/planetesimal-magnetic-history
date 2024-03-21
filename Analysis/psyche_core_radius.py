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
def weight_perc_fe(rhoc,Xs0):
    """

    Parameters
    ----------
    rhoc : float
        core density [kg m^-3]
    Xs0 : float
        initial core sulfur content [wt %]

    Returns
    -------
    XFe : float
        weight percent iron in accreted material

    """
    XFe = (rhoc/rhoa)*((rhoa-rhom)/(rhoc-rhom))*(1-0.01*Xs0)
    return XFe

def core_radius(Xs0):
    """

    Parameters
    ----------
    Xs0 : float
        initial core sulfur content [wt %]

    Returns
    -------
    rc : float
        core radius as fraction of planetesimal radius

    """
    rhoc = fe_fes_density(Xs0)
    XFe = weight_perc_fe(rhoc,Xs0)
    rc = ((XFe/(1-0.01*Xs0))*(rhoa/rhoc))**(1/3)
    return rc

#%% Set parameters
rhom = 3000 # mantle density [kg m^-3] Bryson et. al. (2019)
rhoa = 4000 # bulk density [kg m^-3]

#%% Create arrays for output
Xs0 = np.linspace(0,33,33)

rc = core_radius(Xs0)

#%% Plot result for rc
fig, ax = plt.subplots(1,1)
ax.plot(Xs0,rc)
ax.set(ylabel='$\\frac{r_c}{r}$',xlabel='Initial core sulfur content/ wt %')
ax.hlines(core_radius(27.1),0,27.1,linestyle='dashed',color='black')
ax.vlines(27.1,core_radius(0),core_radius(27.1),linestyle='dashed',color='black')
ax.text(0,0.89,f'{core_radius(27.1):.2f}')
ax.text(27.2,0.63,'27.1')

#%% Plot results for XFe
rhoc = fe_fes_density(Xs0)
Xfe = weight_perc_fe(rhoc,Xs0)
maxfe = weight_perc_fe(fe_fes_density(27.1),27.1)
fig, ax = plt.subplots(1,1)
ax.plot(Xs0,Xfe)
ax.set(ylabel='Iron abundance/ wt %',xlabel='Initial core sulfur content/ wt %')
ax.hlines(maxfe,0,27.1,linestyle='dashed',color='black')
ax.vlines(27.1,0.45,maxfe,linestyle='dashed',color='black')
ax.text(0,0.565,f'{maxfe:.2f}')
ax.text(27.2,0.45,'27.1')