#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare liquidus formalisms
"""

import numpy as np
import matplotlib.pyplot as plt
from fe_fes_liquidus import fe_fes_liquidus_bw, fe_fes_liquidus_linear
from parameters import rhoc, rhom, G, Mr_s, Mr_fe


#Create S array
Xs = np.linspace(0,40,100)
Xsd = Xs/100 #convert wt % to decimal
mrr = Mr_fe/Mr_s
x = Xsd*mrr/(1-Xsd) #mole fraction of FeS


#create pressure array
r = np.array([10,100,250,500])*1e3 #radius [m]
rc = r/2
P = 2*np.pi*G*(rc**2*rhoc+rhom**2*(r**2-rc**2))/1e9 #pressure at centre [GPa]

bw =np.zeros([len(r),100])
Teut = 1260 #Buono & Walker eutectic temp

#Calculate liquidi
linear = fe_fes_liquidus_linear(Xs)
for i, pressure in enumerate(P):
    bw[i,:] = fe_fes_liquidus_bw(Xs,pressure)

#find eutectic S
eutS = []
for i in range(len(r)):
    eutS.append(Xs[bw[i]>Teut][-1])
colors = ['navy','purple','seagreen','cornflowerblue','lightblue']

#plot
plt.figure()
plt.plot(Xs[linear>Teut],linear[linear>Teut],label='linear',color='black')
for i, rad in enumerate(r):
    plt.plot(Xs[bw[i]>Teut],bw[i,bw[i]>Teut],label=f'radius {rad/1e3:.1f} km',color=colors[i])
plt.hlines(1260,0,max(Xs),linestyle='dashed',color='black',label='eutectic temp 1bar - 6GPa')
plt.xlabel('Weight % S')
plt.ylabel('Tl /K')
plt.ylim([1000,1800])
plt.legend()
plt.title('Comparison between linear liquidus and Buono & Walker (2011)')
#plt.savefig('Plots/liquidus_comp_xs.png')

##### Plot as a function of P and x for comparison with Buono & Walker
# plt.figure()
# plt.plot(x,linear,label='linear',color='black')
# for i, rad in enumerate(P):
#     plt.plot(x[bw[i]>Teut],bw[i,bw[i]>Teut],label=f'{rad:.0g} GPa',color=colors[i])
# plt.hlines(1260,0,max(x),linestyle='dashed',color='black',label='eutectic temp 1bar - 6GPa')
# plt.xlabel('Mole fraction of FeS')
# plt.ylabel('Tl /K')
# plt.ylim([1000,1800])
# plt.xlim([0,1])
# plt.legend()
# plt.title('Comparison between linear liquidus and Buono & Walker (2011)')
#plt.savefig('Plots/liquidus_comp.png')