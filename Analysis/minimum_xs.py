#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate minimum Xs for a given planetesimal size and critical melt fraction
"""
#%% Import numpy modules
import numpy as np
import sys
# setting path
sys.path.append('../')
from fe_fes_liquidus import fe_fes_liquidus_bw, central_pressure, fe_fes_density
from parameters import rhom, Mr_s, Mr_fe, rho_exp
from solidus_calc import solidus, liquidus

#%% Choose parameters
r = np.array([120e3]) #m
rc = 0.75*r #core radius
phi = np.array([0.5]) #critical melt fraction
xwater = np.array([0.05]) #wt %

#%% Calculation

#Create S array
Xs = np.linspace(0,40,100)
Xsd = Xs/100 #convert wt % to decimal
mrr = Mr_fe/Mr_s
x = Xsd*mrr/(1-Xsd) #mole fraction of FeS
#Calculate core density
rhoc = fe_fes_density(Xs)*rho_exp
# Calculate solidus
Xsnom = 15 #chosen Xs in Tms and Tml makes 0.01 K difference
Tms = solidus(r,rc,xwater,Xsnom,rhom) 
Tml = liquidus(r,rc,xwater,Xsnom,rhom)
Tms = 1400
Tml = 1800
#pressure
rhom = 3750
P = central_pressure(rhom,rhoc,r,rc) #pressure at centre [GPa]
Tphi = Tms+phi*(Tml-Tms)

bw = fe_fes_liquidus_bw(Xs,P)
Teut = 1260 #Buono & Walker eutectic temp

#find eutectic S
eutS = Xs[bw>Teut][-1]


#find lowest Xs for a given silicate melting
minS=Xs[bw<Tphi][0]

print(f'The eutectic core sulfur content is {eutS:.1f} wt %')
print(f'The minimum initial core sulfur content is {minS:.2f} wt %')
