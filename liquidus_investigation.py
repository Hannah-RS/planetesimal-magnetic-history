#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare liquidus formalisms
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from fe_fes_liquidus import fe_fes_liquidus_bw, fe_fes_liquidus_linear, fe_fes_liquidus_dp
from parameters import rhoc, rhom, G, Mr_s, Mr_fe, Tml, Tms, cpc, alpha_c, Xs_eutectic


#Create S array
Xs = np.linspace(0,40,100)
Xsd = Xs/100 #convert wt % to decimal
mrr = Mr_fe/Mr_s
x = Xsd*mrr/(1-Xsd) #mole fraction of FeS


#create pressure array
r = np.array([10,100,250,500])*1e3 #radius [m]
#r = np.linspace(10,500,200)*1e3
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

#Calculate pressure derivative of liquidus
dTdP = np.zeros([len(r),100])
for i, pressure in enumerate(P):
    dTdP[i,:] = fe_fes_liquidus_dp(Xs,pressure) #divide by 1e9 as pressure derivative was for GPa
Tc = 1400 #estimate for Tc
Delta = dTdP*(rhoc*cpc)/(alpha_c*Tc)

#find lowest Xs for a given silicate melting
phi = np.linspace(0.01,0.5,200)
minS = np.zeros([len(r),len(phi)])
Tphi = Tms+phi*(Tml-Tms)
for i in range(len(r)):
    for j in range(len(phi)):
        minS[i,j]=Xs[bw[i,:]<Tphi[j]][0]

plt.figure()
plt.pcolormesh(Xs[Xs<32],r/1e3,Delta[:,Xs<32])
plt.colorbar(norm=mcolors.LogNorm(),label='$\\frac{dT_L}{dP}$')
plt.ylabel('radius /km')
plt.xlabel('$X_s$ /wt %')

#plot minimum Xs for full melting at differentiation as a function of RCMF
fig = plt.figure(tight_layout=True)
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
ax1.pcolormesh(phi,r/1e3,minS,shading='gouraud',vmin=np.min(minS), vmax=np.max(minS))
ax2.pcolormesh(Tphi,r/1e3,minS,shading='gouraud',vmin=np.min(minS), vmax=np.max(minS))
ax1.set_xlabel('Critical melt fraction')
ax1.set_ylabel('r/km')
ax2.set_xlabel('Differentiation temp /K')
#make the colorbar (took this off the internet it is a bit gross)
norm = mcolors.Normalize(vmin=np.min(minS), vmax=np.max(minS)) 
# creating ScalarMappable
sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
sm.set_array([]) 
plt.colorbar(sm,ax=ax2,label='min X$_s$')
#plt.savefig('Plots/differentiation_Xs.png')

#plot liquidus for different values of r 
colors = ['navy','purple','seagreen','cornflowerblue','lightblue']
plt.figure()
plt.plot(Xs[linear>Teut],linear[linear>Teut],label='linear',color='black')
for i, rad in enumerate(r):
    plt.plot(Xs[bw[i]>Teut],bw[i,bw[i]>Teut],label=f'radius {rad/1e3:.1f} km',color=colors[i%4])
plt.scatter(Xs_eutectic,1260,linestyle='dashed',color='black',label='eutectic temp 1bar - 6GPa',marker='X')
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