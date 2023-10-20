#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extrapolation of reference viscosity from literature to asteroid values for a 
range of sizes and activation volumes
"""
import numpy as np
import sys
sys.path.append('../')
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
pmid = 2/3*np.pi*G*rhom**2*(r**2-rmid**2) #[Pa]

#activation energies from Karato and Wu 1993 and Hirth & Kohlstedt 2003
E = np.array([240,300,325,430,540,570])*1e3 # [J /mol]

#activation volumes Hirth & Kohlstedt 2003
V = np.array([2,10,20])*1e-6 #[m^3/mol]

#find length of arrays for saving values
l = len(r)
m = len(E)
n = len(V)


#%% Extrapolation case 1: Mars, Plesa et. al. 2020
eta0 = 1e21
pref = 3e9
Tref = 1600

#find spread in values
peta_ast = np.zeros([l,m,n])
peta_fk = np.zeros([l,m,n])

for i, p in enumerate(pmid):
    for j, Vin in enumerate(V):
            peta_ast[i,:,j], peta_fk[i,:,j] = viscosity_converter(p, pref, Tref, Vin, E, eta0)
            
# find minimum and maximum reference viscosities
print(f'The minimum asteroid Arrhenius reference viscosity is {np.min(peta_ast):0e} Pas')
print(f'The maximum asteroid Arrhenius reference viscosity is {np.max(peta_ast):0e} Pas')
print(f'The minimum asteroid FK reference viscosity is {np.min(peta_fk):0e} Pas')
print(f'The maximum asteroid FK reference viscosity is {np.max(peta_fk):0e} Pas')

#%% Extrapolation case 2: Ganymede Ruckriemen 2018 - parameters from Table 1
eta0 = np.array([1e19,1e22])
rg = 2634*1e3 #Ganymede radius
rgc = 1000*1e3 #average core radius [m]
rgmid = rgc + (rg-rgc)/2 #mid mantle radius [m]
rhog = 3300 #density of silicate layer [kg m^-3]
pref = 2/3*np.pi*G*rhog**2*(rg**2-rgmid**2)
Tref = 1600

#find spread in values
reta_ast = np.zeros([2,l,m,n])
reta_fk = np.zeros([2,l,m,n])

for i, p in enumerate(pmid):
    for j, Vin in enumerate(V):
        for k, Ein in enumerate(E):
            reta_ast[:,i,k,j], reta_fk[:,i,k,j] = viscosity_converter(p, pref, Tref, Vin, Ein, eta0)
            
# find minimum and maximum reference viscosities
print(f'The minimum asteroid Arrhenius reference viscosity is {np.min(reta_ast):0e} Pas')
print(f'The maximum asteroid Arrhenius reference viscosity is {np.max(reta_ast):0e} Pas')
print(f'The minimum asteroid FK reference viscosity is {np.min(reta_fk):0e} Pas')
print(f'The maximum asteroid FK reference viscosity is {np.max(reta_fk):0e} Pas')

#%% Extrapolation case 3: Scott & Kohlstedt 2006 experimental pararmeters
eta0 = 1e14
Tref = 1523 #but no melt so 1400 may be more appropriate
pref = 300*1e6 #confining pressure [Pa]

#find spread in values
seta_ast = np.zeros([l,m,n])
seta_fk = np.zeros([l,m,n])

for i, p in enumerate(pmid):
    for j, Vin in enumerate(V):
            seta_ast[i,:,j], seta_fk[i,:,j] = viscosity_converter(p, pref, Tref, Vin, E, eta0)
            
# find minimum and maximum reference viscosities
print(f'The minimum asteroid Arrhenius reference viscosity is {np.min(seta_ast):0e} Pas')
print(f'The maximum asteroid Arrhenius reference viscosity is {np.max(seta_ast):0e} Pas')
print(f'The minimum asteroid FK reference viscosity is {np.min(seta_fk):0e} Pas')
print(f'The maximum asteroid FK reference viscosity is {np.max(seta_fk):0e} Pas')