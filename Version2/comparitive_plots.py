#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for plotting the results of multiple runs with variable Ohmic dissipation, core radius, density contrast, adiabatic gradient
Run 2 in the csv is the default parameters
Issue: arrays in each run are different lengths as truncate when f=1 so can't import it all into arrays and use loops :-('
"""
#import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#import varying parameters
vary_p=pd.read_csv('vary_parameters.csv',delimiter=',')
run_num=len(vary_p) #how many runs did you do

#scale time to Myr
from parameters import Myr

#import data on chosen variable from each run 
#import data from npz file
var = 'Rem' #  f, dfdt, Tc, dTcdt or Rem
label = 'Rem' #label for y axis
log = True #conditional log scale
line =True #add horizontal line for Rem=40

def results(var,run):
    npzfile = np.load('Results/run_{}.npz'.format(run))
    var_val= npzfile[var] #values of chosen variable 
    t = npzfile['t'][:-1]/Myr  #time in Myr
    return t, var_val

#have to import all separately as length of time arrays is different due to truncation
t1, var1 = results(var,1)
t2, var2 = results(var,2)
t3, var3 = results(var,3)
t4, var4 = results(var,4)
t5, var5 = results(var,5)
t6, var6 = results(var,6)
t7, var7 = results(var,7)
t8, var8 = results(var,8)
t9, var9 = results(var,9)



# make plot - 4 subplots of chosen variable against t, one plot for each varying parameter
plt.figure(tight_layout=True,figsize=[20,10])
plt.suptitle('Thermal evolution of an asteroid \n magnetic Reynolds number')

plt.subplot(2,2,1)
plt.plot(t1, var1,label='$\Phi_v$={} W/($Km^3$)'.format(vary_p.loc[0,'phiv']))
plt.plot(t2, var2,label='$\Phi_v$={} W/($Km^3$)'.format(vary_p.loc[1,'phiv']))
plt.plot(t3, var3,label='$\Phi_v$={} W/($Km^3$)'.format(vary_p.loc[2,'phiv']))
plt.ylabel('{}'.format(label))
plt.xlabel('t/Myr')
if log == True:
    plt.yscale('log')
else: pass
if line == True:
    plt.hlines(40,xmin=0,xmax=max(t3),linestyle='--',color='k')
else: pass
plt.legend()

plt.subplot(2,2,2)
plt.plot(t4, var4,label='$r_c$={:.0f} km'.format(vary_p.loc[3,'rc']/1e3))
plt.plot(t2, var2,label='$r_c$={:.0f} km'.format(vary_p.loc[1,'rc']/1e3))
plt.plot(t5, var5,label='$r_c$={:.0f} km'.format(vary_p.loc[4,'rc']/1e3))
plt.ylabel('{}'.format(label))
plt.xlabel('t/Myr')
plt.legend()
if log == True:
    plt.yscale('log')
else: pass
if line == True:
    plt.hlines(40,xmin=0,xmax=max(t5),linestyle='--',color='k')
else: pass

plt.subplot(2,2,3)
plt.plot(t6, var6,label='$\delta \\rho /\\rho$={}'.format(vary_p.loc[5,'drho']))
plt.plot(t2, var2,label='$\delta \\rho /\\rho$={}'.format(vary_p.loc[1,'drho']))
plt.plot(t7, var7,label='$\delta \\rho /\\rho$={}'.format(vary_p.loc[6,'drho']))
plt.ylabel('{}'.format(label))
plt.xlabel('t/Myr')
plt.legend()
if log == True:
    plt.yscale('log')
else: pass
if line == True:
    plt.hlines(40,xmin=0,xmax=max(t7),linestyle='--',color='k')
else: pass

plt.subplot(2,2,4)
plt.plot(t8, var8,label='$\Delta$={}'.format(vary_p.loc[7,'Delta']))
plt.plot(t2, var2,label='$\Delta$={}'.format(vary_p.loc[1,'Delta']))
plt.plot(t9, var9,label='$\Delta$={}'.format(vary_p.loc[8,'Delta']))
plt.ylabel('{}'.format(label))
plt.xlabel('t/Myr')
plt.legend()
if log == True:
    plt.yscale('log')
else: pass
if line == True:
    plt.hlines(40,xmin=0,xmax=max(t9),linestyle='--',color='k')
else: pass

plt.savefig('Plots/{}.png'.format(var),dpi=300) #think of a more systematic naming system later
