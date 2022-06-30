#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot the output of runs 11, 12, 13 which turn off different terms
Run 11 - No latent heat or secular cooling
Run 12 - No latent heat or GPE release
Run 13 - No secular cooling or GPE release
Run 14 - No secular cooling
Run 15 - Mo GPE
Run 16 - all terms
"""
#import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


#scale time to Myr
from parameters import Myr

def make_subplot(var,label,log=False,line=False,Rem=False):
    """
    

    Parameters
    ----------
    var : string
        variable to be plotted, options are f, dfdt, Tc, dTcdt, Rem
    label: string, 
         label for y axis

    Returns
    -------
    A subplot with three lines one for each run of chosen variable against time

    """
    

    def results(var,run):
        npzfile = np.load('Results/run_{}.npz'.format(run))
        var_val= npzfile[var] #values of chosen variable 
        if Rem == False:
            t = npzfile['t']/Myr  #time in Myr
        elif Rem == True: #for magnetic Reynolds number there is one less point
            t = npzfile['t'][:-1]/Myr  #time in Myr
        else: raise ValueError('Rem must be boolean')
        return t, var_val

    #have to import all separately as length of time arrays is different due to truncation
    t1, var1 = results(var,11)
    t2, var2 = results(var,12)
    t3, var3 = results(var,13)
    t4, var4 = results(var,14)
    t5, var5 = results(var,15)
    t6, var6 = results(var,16)

    plt.plot(t1, var1,label='Qg')
    plt.plot(t2, var2,label='Qs')
    plt.plot(t3, var3,label='Ql')
    plt.plot(t4, var4,label='Ql + Qg')
    plt.plot(t5, var5,label='Ql + Qs')
    plt.plot(t6, var6,label='Ql + Qs + Qg')
    plt.ylabel('{}'.format(label))
    plt.xlabel('t/Myr')
    plt.legend()
    if log == True:
        plt.yscale('log')
    else: pass
    if line == True:
        plt.hlines(40,xmin=0,xmax=max(t3),linestyle='--',color='k')
    else: pass

    return None

# make plot - 4 subplots of chosen variable against t, one plot for each varying parameter
plt.figure(tight_layout=True,figsize=[20,10])
plt.suptitle('Thermal evolution of an asteroid \n including different terms')

plt.subplot(2,2,1)
make_subplot('f','fractional inner core radius')


plt.subplot(2,2,2)
make_subplot('Tc','Core temperature/K')

plt.subplot(2,2,3)
make_subplot('dfdt','df/dt',log=True)

plt.subplot(2,2,4)
make_subplot('Rem','magnetic Reynolds number',log=True,line=True,Rem=True)

plt.savefig('Plots/single_term.png',dpi=300) #think of a more systematic naming system later