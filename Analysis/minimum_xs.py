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
from fe_fes_liquidus import fe_fes_liquidus_bw, central_pressure
from parameters import rhoc, rhom, Mr_s, Mr_fe, Tml, Tms

#%% Choose parameters
r = np.array([300e3]) #m
phi = np.array([0.3]) #critical melt fraction

#%% Calculation

#Create S array
Xs = np.linspace(0,40,100)
Xsd = Xs/100 #convert wt % to decimal
mrr = Mr_fe/Mr_s
x = Xsd*mrr/(1-Xsd) #mole fraction of FeS

#pressure, c
rc = r/2
P = central_pressure(rhom,rhoc,r,rc) #pressure at centre [GPa]
Tphi = Tms+phi*(Tml-Tms)

bw = fe_fes_liquidus_bw(Xs,P)
Teut = 1260 #Buono & Walker eutectic temp

#find eutectic S
eutS = Xs[bw>Teut][-1]


#find lowest Xs for a given silicate melting
minS=Xs[bw<Tphi][0]

print(f'The eutectic core sulfur content is {eutS:.1f} wt %')
print(f'The minimum initial core sulfur content is {minS:.1f} wt %')
