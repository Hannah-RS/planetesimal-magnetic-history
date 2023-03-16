#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare liquidus formalisms
"""

import numpy as np
import matplotlib.pyplot as plt
from fe_fes_liquidus import fe_fes_liquids_bw, fe_fes_liquidus
from parameters import rhoc, rhom, G, Mr_s, Mr_fe


#Create S array
Xs = np.linspace(0,80,100)
Xsd = Xs/100 #convert wt % to decimal
mole_s = Xsd*1000/Mr_s
mole_fe = (1-Xsd)*1000/Mr_fe

x = mole_s/(mole_s+mole_fe)

#create pressure array
r = np.linspace(100,1500,5)*1e3 #radius [m]
rc = r/2
P = 2*np.pi*G*(rc**2*rhoc+rhom**2*(r**2-rc**2))/1e9 #pressure at centre [GPa]

bw =np.zeros([5,100])
#Calculate liquidi
linear = fe_fes_liquidus(Xs)
for i, pressure in enumerate(P):
    bw[i,:] = fe_fes_liquids_bw(Xs,pressure)

colors = ['navy','purple','seagreen','cornflowerblue','lightblue']
#plot
plt.figure()
plt.plot(Xs,linear,label='linear',color='black')
for i, rad in enumerate(r):
    plt.plot(Xs,bw[i,:],label=f'{rad/1e3:.0f} km radius',color=colors[i])
plt.hlines(1260,0,max(Xs),linestyle='dashed',color='black',label='eutectic 1bar - 6GPa')
plt.xlabel('Xs /wt %')
plt.ylabel('Tl /K')
plt.ylim(bottom=1000)
plt.legend()
plt.title('Comparison between linear liquidus and Buono & Walker (2011)')
#plt.savefig('')