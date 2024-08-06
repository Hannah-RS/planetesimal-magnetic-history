#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extrapolation of reference viscosity from literature to asteroid values for a 
range of sizes and activation volumes
"""
import numpy as np
import sys
sys.path.append('../')
from solidus_calc import pmid_calc
from parameters import R, G, Tms, rhoc, rhom

def viscosity_converter(p,pref,Tref,V,E,eta0):
    """
    
    Parameters
    ----------
    p : float
        asteroid pressure [Pa]
    pref : float
        reference pressure [Pa]
    Tref : float
        reference temperature [K]
    V : float
        activation volume [m^3 mol^-1]
    E : float
        activation energy [J mol^-1]
    eta0 : float
        reference viscosity in literature [Pas]

    Returns
    -------
    eta_ast : float
        reference viscosity at mantle solidus and asteroid pressure [Pas]
    eta_fk : float
        reference viscosity at mantle solidus and asteroid pressure from Frank Kamenetskii approximation[Pas]

    """
    #Arrhenius viscosity
    eta_ast = eta0*np.exp((E+p*V)/(R*Tms)-(E+pref*V)/(R*Tref))  #Eqn 2.6 Noack & Breuer 2013
    # Frank-Kamenetskii viscosity
    eta_fk = eta0*np.exp((E+pref*V)/(R*Tref**2)*(Tref-Tms)+V/(R*Tref)*(p-pref)) #Eqn. A5 Noack & Breuer 2013
    
    return eta_ast, eta_fk


#calculate pressures for a range of asteroid sizes at mid-mantle
r = np.array([100,200,300,400,500])*1e3 #[m]
rc = r/2
rmid = 0.75*r
pmid = pmid_calc(rhom,rhoc,rc,rmid,r) # [Pa]

#activation energies from Karato and Wu 1993 and Hirth & Kohlstedt 2003
E = np.array([240,300,325,430,540,570])*1e3 # [J /mol]

#activation volumes Hirth & Kohlstedt 2003
V = np.array([2,10,20])*1e-6 #[m^3/mol]

#find length of arrays for saving values
l = len(r)
m = len(E)
n = len(V)


#%% Extrapolation case 1: Mars, Plesa et. al. 2020
eta0 = 1e21 #[Pas]
pref = 3e9  # [GPa]
Tref = 1600# [K]

#find spread in values
peta_ast = np.zeros([l,m,n])
peta_fk = np.zeros([l,m,n])

for i, p in enumerate(pmid):
    for j, Vin in enumerate(V):
            peta_ast[i,:,j], peta_fk[i,:,j] = viscosity_converter(p, pref, Tref, Vin, E, eta0)
            
# find minimum and maximum reference viscosities
print(f'From Mars: The minimum asteroid Arrhenius reference viscosity is {np.min(peta_ast):0e} Pas')
print(f'From Mars: The maximum asteroid Arrhenius reference viscosity is {np.max(peta_ast):0e} Pas')
print(f'From Mars: The minimum asteroid FK reference viscosity is {np.min(peta_fk):0e} Pas')
print(f'From Mars: The maximum asteroid FK reference viscosity is {np.max(peta_fk):0e} Pas')

#%% Extrapolation case 2: Ganymede Ruckriemen 2018 - parameters from Table 1
eta0 = np.array([1e19,1e22])
rg = 2634*1e3 #Ganymede radius [m]
rgc = 1000*1e3 #average core radius [m]
rgm = 1784*1e3 #Ganymede mantle radius [m]

rgmid = rgc + rgm/2 #mid mantle radius [m]
rhomg = 3300 #density of silicate layer [kg m^-3] from Ruckriemen 2018
rhocg = 6072 #density of core [kg m^-3] from Ruckriemen 2018
rhoi = 1200 # density of ice [kg m^-3] from Ruckriemen 2018
# Ganymede has an ice layer!
pref = 2*np.pi*G/3*(rhomg**2*(rgm**2-rgmid**2)+rhoi**2*(rg**2-rgm**2)+2*rhomg*(rhocg-rhomg)*rgc**3*((1/rgmid)-(1/rgm)))\
        -2*np.pi*G/3*(2*rhoi*(rgc**3*rhocg+(rgm**3-rgc**3)*rhomg-rgm**3*rhoi)*((1/rg)-(1/rgm)))
Tref = 1600 #[K]

#find spread in values
reta_ast = np.zeros([2,l,m,n])
reta_fk = np.zeros([2,l,m,n])

for i, p in enumerate(pmid):
    for j, Vin in enumerate(V):
        for k, Ein in enumerate(E):
            reta_ast[:,i,k,j], reta_fk[:,i,k,j] = viscosity_converter(p, pref, Tref, Vin, Ein, eta0)
            
# find minimum and maximum reference viscosities
print(f'From Ganymede: The minimum asteroid Arrhenius reference viscosity is {np.min(reta_ast):0e} Pas')
print(f'From Ganymede: The maximum asteroid Arrhenius reference viscosity is {np.max(reta_ast):0e} Pas')
print(f'From Ganymede: The minimum asteroid FK reference viscosity is {np.min(reta_fk):0e} Pas')
print(f'From Ganymede: The maximum asteroid FK reference viscosity is {np.max(reta_fk):0e} Pas')

#%% Extrapolation case 3: Scott & Kohlstedt 2006 experimental pararmeters
eta0 = 1e14 #[Pas]
Tref = 1523 #but no melt so 1400 may be more appropriate, [K]
pref = 300*1e6 #confining pressure [Pa]

#find spread in values
seta_ast = np.zeros([l,m,n])
seta_fk = np.zeros([l,m,n])

for i, p in enumerate(pmid):
    for j, Vin in enumerate(V):
            seta_ast[i,:,j], seta_fk[i,:,j] = viscosity_converter(p, pref, Tref, Vin, E, eta0)
            
# find minimum and maximum reference viscosities
print(f'From experiments: The minimum asteroid Arrhenius reference viscosity is {np.min(seta_ast):0e} Pas')
print(f'From experiments: The maximum asteroid Arrhenius reference viscosity is {np.max(seta_ast):0e} Pas')
print(f'From experiments: The minimum asteroid FK reference viscosity is {np.min(seta_fk):0e} Pas')
print(f'From experiments: The maximum asteroid FK reference viscosity is {np.max(seta_fk):0e} Pas')