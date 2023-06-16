#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot of flux and temperature evolution with time - rolled back version for old runs for transfer
The first plot from make_plots.py
"""
#import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import ScalarFormatter
import seaborn as sns

#choose your run
run=33
conduction = False #did the mantle start conducting before the core solidified

#scale time to Myr
from parameters import Myr, r, Tm0, Tsolidus

#import data from npz file
npzfile = np.load('Results/run_{}.npz'.format(run))
Tm = npzfile['Tm_base'] 

Tc= npzfile['Tc'] 
t = npzfile['t'] #time in s
Flux = npzfile['Flux']


Fs = Flux[0]
Fcmb = Flux[1]
Fad = Flux[2]
Frad = Flux[3]

t_plot = t/Myr

#import label info - read in from correct row in csv
run_info = pd.read_csv('run_info3.csv',delimiter=',')

row = run_info[run_info['run']==run]
r = row.iloc[0,1] #radius [m]
Tsolidus = row.iloc[0,2] #temperature for onset of solidification [K]
Tm0 = row.iloc[0,3] # initial mantle and core temp [K]
tstart = row.iloc[0,4]
tend = row.iloc[0,5] # max possible time of simulation [Myr]
tstep = row.iloc[0,6] #max timestep [Myr]
tsolid = row.iloc[0,7] #time at which solidifcation finishes [Myr] if tsolid == tend then core may not have finished solidifying
if conduction == True:
    cond_i= int(row.iloc[0,8]) #index in array at which mantle started conducting
    cond_t = t[cond_i]/Myr #time at which mantle switched to conduction
dr=row.iloc[0,10] #cell spacing
dt =row.iloc[0,11] #T_profile output frequency

################### Main Plot #########################################
# make  log-log plot - similar to Bryson 2019
with sns.plotting_context('talk'):
    plt.figure(tight_layout=True,figsize=[10,7])
    plt.suptitle('Thermal evolution of a {:.0f}km asteroid'.format(r/1e3, Tm0, Tsolidus,run))
    
    xmin=tstart
    
    #temperatures as function of time

    plt.semilogx(t_plot,Tm,label='mantle',color='red')
    plt.semilogx(t_plot,Tc,label='core',color='black')
    if conduction == True:
        plt.vlines(cond_t,ymin=min(Tm),ymax=1600,color='black',linestyle='--',label='conduction')
    #plt.xlim([5,10])
    #plt.ylim([1580,1610])
    #plt.xlabel('Time/ Myr')
    #plt.xlim([xmin,500])  #use these limits when comparing runs
    #plt.ylim([1400,1650]) #use these limits when comparing runs
    plt.ylabel('T/K')
    plt.locator_params('y',nbins=4)
    #plt.ylim([450,1700])
    ax=plt.gca()
    ax.xaxis.set_major_formatter(ScalarFormatter())
    plt.legend(loc='upper right')
    plt.xlabel('Time/ Myr')
   # plt.savefig('Plots/transfer_run{}.png'.format(run),dpi=600)
    # #fluxes as function of time
    # plt.subplot(2,1,2)
    # plt.loglog(t_plot,Fs,label='$F_s$',color='green')
    # plt.loglog(t_plot,Fcmb,label='$F_{CMB}$',color='navy',linestyle='dashed')
    # plt.loglog(t_plot,Fad,label='$F_{ad}$',color='orange',linestyle='dotted')
    # plt.loglog(t_plot,Frad,label='$F_{rad}$',color='purple',linestyle='dashdot')
    # plt.xlabel('Time/ Myr')
    # plt.ylim([1e-3,1e2])   
    # #plt.xlim([5,10])
    # # plt.xlim([xmin,500])   #use these limits when comparing runs
    # plt.ylabel('Flux/ W$m^{-2}$')
    # plt.legend(loc='upper right',ncol=2)
    # ax=plt.gca()
    # ax.xaxis.set_major_formatter(ScalarFormatter())
    #plt.savefig('Plots/transfer_run{}.png'.format(run),dpi=600)
    
