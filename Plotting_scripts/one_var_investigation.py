#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make scatter plots to compare effects of one variables on the following parameters:
Dynamo start
Dynamo stop
Duration
Ng l10
Ngl100
Peak B
End of convection
fconv_T
Diff t
Core solid
Peak mantle T
Peak core T
No filtering is applied to the dataset
"""
from scatter_function import make_scatter, make_sub_scatter
import sys
# setting path
sys.path.append('../')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

folder = 'Fullrun2'
data = pd.read_csv(f'../Results_combined/{folder}/all_sucess_info.csv',delimiter=',',skiprows=[1],header=0,index_col=False)
data['r']=data['r']/1e3 #rescale to km
save = True

varlist =['r','Xs_0','rcmf','frht','etal','eta0','Fe0']
varlabel =['radius /km','$X_{S,0}$ /wt%','rcmf','frht','liquid viscosity /Pas','reference viscosity /Pas','$^{60}Fe/^{56}Fe$']
######################## Dynamo start, stop, duration #########################
# This plot is kinda gross
ylist = ['_on','_off','_dur']
ylab = ['Dynamo onset /Myr','Dynamo cessation /Myr','Dynamo duration /Myr']
names = ['mac','cia','comp']
fig, ax = plt.subplots(nrows=len(ylist),ncols=len(varlist),sharey='row',sharex='col',figsize=[15,15]) 
for j, lab in enumerate(ylist):
    i=0
    for var, label in zip(varlist,varlabel):
        for name in names:
            ax[j,i].scatter(data[var],data[name+ylist[j]],label=name)
        if (i==4) | (i==5):
            ax[j,i].set_xscale('log')
        ax[j,i].set_xlabel(label)
        i=i+1
    ax[j,0].set_ylabel(ylab[j])
ax[0,0].legend()
if save == True:
    plt.savefig(f'../Plots/{folder}dynamo_duration_1var.png')
    
######################### Intermittence ##################################
ylist = ['ngl10','ngl100','ngg100']
ylab = ['Gaps <10 Myr','10Myr < gap < 100 Myr','gap >100Myr']
names = ['mac_','cia_','comp_']
for name in names:
    data[name+'ngg100']=data[name+'n']-data[name+'ngl100']-data[name+'ngl10']-1
fig, ax = plt.subplots(nrows=len(ylist),ncols=len(varlist),sharey='row',sharex='col',figsize=[15,15]) 
for j, lab in enumerate(ylist):
    i=0
    for var, label in zip(varlist,varlabel):
        for name in names:
            ax[j,i].scatter(data[var],data[name+lab],label=name)
        if (i==4) | (i==5):
            ax[j,i].set_xscale('log')
        ax[j,i].set_xlabel(label)
        i=i+1
    ax[j,0].set_ylabel(ylab[j])
ax[0,0].legend()
if save == True:
    plt.savefig(f'../Plots/{folder}dynamo_intermittence_1var.png')
    
####################### Peak field strength ###############################
ylist = ['maxB_mac','maxB_cia','maxB_comp']
ylab = ['MAC /$\\mu$T','CIA /$\\mu$T','comp /$\\mu$T']

fig, ax = plt.subplots(nrows=len(ylist),ncols=len(varlist),sharey='row',sharex='col',figsize=[15,15]) 
for j, lab in enumerate(ylist):
    i=0
    for var, label in zip(varlist,varlabel):
        ax[j,i].scatter(data[var],data[lab]/1e-6,label=name)
        if (i==4) | (i==5):
            ax[j,i].set_xscale('log')
        ax[j,i].set_xlabel(label)
        i=i+1
    ax[j,0].set_ylabel(ylab[j])
ax[0,0].legend()
if save == True:
    plt.savefig(f'../Plots/{folder}maxB_1var.png')
######################## Convective shut off ################################
varlist =['r','Xs_0','rcmf','frht','etal','eta0','Fe0']
varlabel =['radius /km','$X_{S,0}$ /wt%','rcmf','frht','liquid viscosity /Pas','reference viscosity /Pas','$^{60}Fe/^{56}Fe$']

fig, ax = plt.subplots(nrows=3,ncols=len(varlist),sharey='row',sharex='col',figsize=[15,15]) 
#first time of conduction
i=0
for var, label in zip(varlist,varlabel):
    ax[0,i].scatter(data[var],data['lconv_t'])
    if (i==4) | (i==5):
        ax[0,i].set_xscale('log')
    ax[0,i].set_xlabel(label)
    i=i+1
ax[0,0].set_ylabel('First time of conduction /Myr ')
#last time of convection
i=0
for var, label in zip(varlist,varlabel):
    ax[1,i].scatter(data[var],data['fcond_t'])
    if (i==4) | (i==5):
        ax[1,i].set_xscale('log')
    ax[1,i].set_xlabel(label)
    i=i+1
ax[1,0].set_ylabel(' Last time of convection /Myr ')
#duration of buffering
i=0
for var, label in zip(varlist,varlabel):
    ax[2,i].scatter(data[var],data['fcond_t']-data['lconv_t'])
    if (i==4) | (i==5):
        ax[2,i].set_xscale('log')
    ax[2,i].set_xlabel(label)
    i=i+1
ax[2,0].set_ylabel('Duration of buffering /Myr')

########################Differentiation time #################################
ylist = ['diff_time','tsolid']
ylab = ['Differentiation time /Myr','Solidification time /Myr']

fig, ax = plt.subplots(nrows=len(ylist),ncols=len(varlist),sharey='row',sharex='col',figsize=[15,15]) 
for j, lab in enumerate(ylist):
    i=0
    for var, label in zip(varlist,varlabel):
        ax[j,i].scatter(data[var],data[lab],label=name)
        if (i==4) | (i==5):
            ax[j,i].set_xscale('log')
        ax[j,i].set_xlabel(label)
        i=i+1
    ax[j,0].set_ylabel(ylab[j])
ax[0,0].legend()
if save == True:
    plt.savefig(f'../Plots/{folder}difft_1var.png')
    
######################## Peak mantle and core temperature #################################
ylist = ['peakT','peak_coreT']
ylab = ['Peak mantle temp /K','Peak core temp /K']

fig, ax = plt.subplots(nrows=len(ylist),ncols=len(varlist),sharey='row',sharex='col',figsize=[15,15]) 
for j, lab in enumerate(ylist):
    i=0
    for var, label in zip(varlist,varlabel):
        ax[j,i].scatter(data[var],data[lab],label=name)
        if (i==4) | (i==5):
            ax[j,i].set_xscale('log')
        ax[j,i].set_xlabel(label)
        i=i+1
    ax[j,0].set_ylabel(ylab[j])
ax[0,0].legend()
if save == True:
    plt.savefig(f'../Plots/{folder}peakT_1var.png')