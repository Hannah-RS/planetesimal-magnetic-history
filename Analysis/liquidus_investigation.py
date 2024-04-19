#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare liquidus formalisms
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import sys
# setting path
sys.path.append('../')
from fe_fes_liquidus import fe_fes_liquidus_bw, fe_fes_liquidus_linear,\
    fe_fes_liquidus_dp, central_pressure
from parameters import rhoc, rhom, G, Mr_s, Mr_fe, Tml, Tms, cpc, alpha_c, Xs_eutectic


#Create S array
Xs = np.linspace(0,40,100)
Xsd = Xs/100 #convert wt % to decimal
mrr = Mr_fe/Mr_s
x = Xsd*mrr/(1-Xsd) #mole fraction of FeS


#create pressure array
#r = np.array([10,100,250,300,500])*1e3 #radius [m]
r = np.array([200e3]) #for one radius just use this line
#r = np.linspace(100,500,3)*1e3
rc = r/2
P = central_pressure(rhom,rhoc,r,rc) #pressure at centre [GPa]

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
#phi = np.linspace(0.2,0.5,200)
phi = np.array([0.3]) #for one melt fraction just use this line
minS = np.zeros([len(r),len(phi)])
Tphi = Tms+phi*(Tml-Tms)
for i in range(len(r)):
    for j in range(len(phi)):
        minS[i,j]=Xs[bw[i,:]<Tphi[j]][0]

if len(r) == 1: #just trying to get one value
    print(minS)
plt.figure()
plt.pcolormesh(Xs[Xs<Xs_eutectic],r/1e3,Delta[:,Xs<Xs_eutectic])
plt.colorbar(norm=mcolors.LogNorm(),label='$\\frac{dT_L}{dP}$')
plt.ylabel('radius /km')
plt.xlabel('$X_s$ /wt %')

#%%
#plot minimum Xs for full melting at differentiation as a function of RCMF
fig = plt.figure(tight_layout=True)
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
ax1.pcolormesh(phi,r/1e3,minS,shading='gouraud',vmin=np.min(minS), vmax=np.max(minS))
ax2.pcolormesh(Tphi,r/1e3,minS,shading='gouraud',vmin=np.min(minS), vmax=np.max(minS))
ax1.set_xlabel('Critical melt fraction, $\phi_C$')
ax1.set_ylabel('planetesimal radius /km')
ax2.set_xlabel('Differentiation temp /K')
#make the colorbar (took this off the internet it is a bit gross)
norm = mcolors.Normalize(vmin=np.min(minS), vmax=np.max(minS)) 
# creating ScalarMappable
sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
sm.set_array([]) 
plt.colorbar(sm,ax=ax2,label='minimum X$_{S,0}$')
#plt.savefig('../Plots/Icarus_paper/differentiation_Xs.png',dpi=500)

#%%
#plot liquidus for different values of r 
colors = ['#0292D7','#BB4DA7','#361AE5','cornflowerblue','lightblue']
plt.figure()
plt.plot(Xs[linear>Teut],linear[linear>Teut],label='linear',color='black')
for i, rad in enumerate(r):
    plt.plot(Xs[bw[i]>Teut],bw[i,bw[i]>Teut],label=f'R={rad/1e3:.0f} km',color=colors[i%4])
plt.scatter(Xs_eutectic,1260,linestyle='dashed',color='black',label='eutectic temperature 1bar - 6GPa',marker='X')
plt.xlabel('Core sulfur content /wt % S')
plt.ylabel('Liquidus temperature /K')
plt.ylim([1200,1800])
plt.legend()
#plt.savefig('../Plots/liquidus_comp_xs.png',dpi=500)


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