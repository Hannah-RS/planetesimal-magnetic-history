#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for plotting the results of the integration - temperature and fractional inner core radius with time
"""
#import modules
import matplotlib.pyplot as plt
import numpy as np

#choose your run
run=36
Qcmb=100 #GW

#scale time to Myr
from parameters import Myr

#import data from npz file
npzfile = np.load('Results/run_{}.npz'.format(run))
Tc = npzfile['Tc'] #core temp in K
t = npzfile['t'] #time in s
f = npzfile['f'] #fractional inner core radius
dfdt = npzfile['dfdt']*Myr # rate of change of fractional inner core radius per Myr
dTcdt = npzfile['dTcdt']*Myr # rate of change in core temperature per Myr 
Rem = npzfile['Rem'] # magnetic Reynolds number


t_plot = t/Myr


#import label info
from parameters import rc

# make plot
plt.figure(tight_layout=True)
plt.suptitle('Thermal evolution of a {:.0f}km asteroid \n Qcmb={}GW'.format(rc/1e3,Qcmb))
plt.subplot(2,2,1)
plt.plot(t_plot, Tc)
plt.ylabel('$T_C$/K')
plt.xlabel('t/Myr')
plt.subplot(2,2,2)
plt.plot(t_plot,f)
plt.xlabel('t/Myr')
plt.ylabel('f')
plt.subplot(2,2,3)
plt.plot(t_plot[:-1],Rem) #last point truncated due to double scalar error
plt.xlabel('t/Myr')
plt.ylabel('Rem')
plt.hlines(y=40,xmin=0,xmax=max(t_plot),color='k',linestyle='--')
plt.yscale('log')
plt.subplot(2,2,4)
plt.plot(t_plot,dfdt)
plt.xlabel('t/Myr')
plt.ylabel('df/dt /Myr')
plt.yscale('log')
plt.savefig('Plots/run{}.png'.format(run),dpi=300) #think of a more systematic naming system later