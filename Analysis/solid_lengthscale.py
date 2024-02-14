#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Core solidification
Lengthscale comparison between concentric inward and cumulate in centre
"""
import numpy as np
import matplotlib.pyplot as plt

@np.vectorize
def r1(x,f):
    """
    outer solid shell inner radius

    Parameters
    ----------
    x : float
        fraction of mass in outer shell
    f : float
        fractional inner core radius for concentric inward

    Returns
    -------
    r1 : float
        fractional inner radius of outer solid shell

    """
    r1 = (1-(1-x)*(1-f**3))**(1/3)
    return r1

@np.vectorize
def r2(x,f):
    """
    inner solid core radius

    Parameters
    ----------
    x : float
        fraction of mass in outer shell
    f : float
        fractional inner core radius for concentric inward

    Returns
    -------
    r2 : float
        fractional inner solid core radius

    """
    r2 = (x*(1-f**3))**(1/3)
    return r2

f = np.linspace(0.7,0.99,5)
x = np.linspace(0,1,100)

lratio = np.zeros([len(f),len(x)])

for i, fval in enumerate(f):
    lratio[i,:] = (r1(x,fval)-r2(x,fval))/fval

#%% Plot results
colors = ['black','navy','cornflowerblue','skyblue','mediumturquoise']

#Range of central masses
plt.figure()
for i in range(len(f)):
    plt.plot(x,lratio[i,:],label=f'$r_i$={f[i]:.3f}',color=colors[i])
plt.xlabel('Fraction of solidified mass in core centre')
plt.ylabel('(r$_2$-r$_1$)/r$_i$')
plt.legend()
plt.ylim([0,1])

#%% All mass in centre for range of f
xvals = [0.25,0.5,0.75]
plt.figure()
for xval in xvals:
    plt.plot(f,(r1(xval,f)-r2(xval,f))/f,label=xval)
plt.xlabel('Fraction inner core radius')
plt.ylabel('(r$_2$-r$_1$)/r$_i$')
plt.legend()
plt.ylim([0,1])