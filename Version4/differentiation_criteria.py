#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Investigate affect of stagnant lid criterion on differentiation time and temperature
"""
import numpy as np
import matplotlib.pyplot as plt
from parameters import Myr

data = np.loadtxt('lid_test.csv',delimiter=',',skiprows=1)

t = data[:,0]/Myr #time in s
T = data[:,1] #temp in K at onset of differentiation
f = data[:,2] #fraction of stagnant lid thickness to total radius

melt = np.array([1480,1520,1560,1600])
colors=['black','darkblue','cornflowerblue','lightblue']

plt.figure()
plt.scatter(f,T)
for val, color in zip(melt,colors):
    plt.hlines(val,0.001,0.9,label=f'{(val-1400)/4:.0f}% melt',linestyle='dashed',color=color,linewidth=2)
plt.xlabel('$\delta_0$/r')
plt.ylabel('T at onset/K')
plt.xscale('log')
plt.legend()


plt.figure()
plt.scatter(f,t)
plt.xlabel('$\delta_0$/r')
plt.ylabel('t at onset/ Myr')
plt.xscale('log')

