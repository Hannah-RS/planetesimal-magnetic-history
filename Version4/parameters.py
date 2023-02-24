#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parameters for coupled mantle and core evolution model. 
All numbers for the mantle from Bryson et al. (2019) Constraints on asteroid 
magnetic field evolution and the radii of meteorite parent bodies from thermal modelling.
All numbers from the core are from Bryson et al. (2019) with value from Nimmo (2009) 
"Energetics of asteroid dynamos and the role of compositional convection." given in brackets after 
or used where there was no equivalent value in Bryson 2019

If numbers from a different source it is commented next to it
"""

import numpy as np

#Constants
G=6.67e-11 # gravitational constant [N kg^-2 m^2]
year=60*60*24*365 #number of seconds in a year
Myr=1e6*year

# Size of body
r = 400e3 # radius of asteroid [m]
V = 4/3*np.pi*r**3
rc = r/2 #radius of core [m]
dr = 500 # size of cells [m]
out_interval =10*Myr #how often do you want temp profiles to be output
save_interval = 0.01*Myr # how often do you want each variable to be saved

#Surface parameters
Ts=200 # surface temperature, modelled as fixed [K]
As = 4*np.pi*r**2 # surface area [m^2]

#Undifferentiated parameters
ka = 2.16 # [W /m /K] currently same as mantle
rhoa = 3000 # kg m^-3 density of undifferentiated material
cpa = 850 # heat capacity [J /kg /K]
Tref = 1800 # viscosity reference temperature (Dodds 2021) [K] 

# Mantle parameters
cpm = 850 # heat capacity [J /kg /K] (Dodds et al.)
Lm = 400e3 #latent heat of mantle [J /kg]
Tml = 1800 # liquidus [K]
Tms = 1400 # solidus [K]
rhom = 3000 # density [kg m^-3]
alpha_m = 4e-5 # thermal expansivity of mantle [/K]
kappa = 9e-7 # thermal diffusivity of silicate [m^2 /s] - 4 options are given in Bryson (2019) have picked the one used in the figures
Rac = 1000  #critical Rayleigh number 
E = 300e3 # activation energy [J /mol]
R = 8.31 # gas constant [J /K /mol]
Tm0 = 1600 # intial mantle temp - from beginning of stage 2 in Figure 1 of Bryson (2019)
c1 = 8 # constant in boundary layer thickness Dodds (2021)
t_transition = 17*Myr #time to transition from RaH to normal Ra

#radiogenic heating
h0Al = 0.355 # [W/ kg] heating rate of Al^26 at t=0 (Dodds 2021)
Al0 = 5e-5 # 26Al/27Al ratio in accreting material (Dodds 2021)
XAl_a = 0.014 # abundance of Al in accreting material [wt % /100] (Dodds 2021)
thalf_al = 0.717*Myr # half life of Al26 [s] (Dodds 2021)
h0Fe = 0.04 # [W/ kg] heating rate of 60Fe at t=0 (Neumann 2012, converted)
Fe0 = 6e-7 # 60Fe/56FE ratio in accreting material (Cook 2021)
XFe_a = 0.22 # abundance of Fe in accreting material [wt % /100] (Neumann 2012)
thalf_fe = 2.6*Myr # half life of 60Fe [s] (Neumann 2012)

Xs_0 = 32 # initial wt % sulfur in core (needs a citation)
XFe_d = 1 - Xs_0/100 #abundance of Fe in core assuming S is only other significant phase [wt %]
Xs_to_core = XFe_a/(1-Xs_0/100)-XFe_a # wt % of S from accreted body that went into core
XAl_d = XAl_a/(1-XFe_d-Xs_to_core)
#viscosity models
default ='Bryson2' #default viscosity model
# Bryson 2019
eta0 = 1e21 # reference viscosity [Pa s] - assumed constant in this model
eta0_50 = 1e14 #viscosity of material at 50% melting [Pas]
eta_r50 = 9.358045457250838e+16 #viscosity of material at 50% melting [Pas] for Arrhenius model
alpha_n = 25 # constant in viscosity model
T0eta = 1400 # reference temperature [K]
km = kappa*cpm*rhom # thermal conductivity of silicate [W /m /K]
Tm50 = 1600 # 50% melting temperature [K]

#Sterenborg and Crowley 2012
Tcrit = Tms+(Tml-Tms)/2 #50% melting temperature

# Core parameters from Bryson (2019)
#values from Nimmo 2009 are shown after the units in square brackets if they disagree
rhoc=7800 # density of core [kg m^-3] 7019
cpc=850 # heat capacity of core [J kg^-1 K^-1] 835
kc=30 # thermal conducitivity of core material [Wm^-1 K^-1] 30
alpha_c=9.2e-5 #[K^-1] 9.2e-5
Lc = 270e3 # latent heat of core [J kg^-1] 750,000
drho=0.025 # \delta rho/rho 0.05 Nimmo (2009)
Delta=1.2 #dTm/dP/dT/dP 1.2 Nimmo (2009)
eta_c =0.01 # viscosity of core [Pa s] Dodds (2021)
bpart = 0.5 #buoyancy partitioning coefficient Nichols (2021) but based on Aubert 2009 - might want to investigate
Xs_0 = 32 # initial wt % sulfur in core (needs a citation)
Xs_eutectic = 32 # eutectic wt% S

# Initial conditions
Tc0=Tm0# intial CMB temp [K] is same as intial mantle temp from Fig 1 in Bryson (2019)
f0=0.01 #initial fractional inner core radius 

#constants for magnetic Reynolds number
lambda_mag = 1.3 # magnetic diffusivity [m^2 s^{-1}] 2
Omega = 1.75e-4 # rotational frequency (10 hr period) [s^{-1}] 6hr period in N09
beta=0.85 

# Calculated mantle parameters
cpm_p = cpm*(1+(Lm/(cpm*(Tml-Tms)))) # modified mantle heat capacity accounting for latent heat
Vm = 4/3*np.pi*(r**3-rc**3) # volume of magma ocean, assumed to be same as volume of mantle
gamma = E/(R*T0eta**2)
g = G*(Vm*rhom+4/3*np.pi*rc**3*rhoc)/r**2 # surface gravity 

#Calculated core parameters
D= np.sqrt(3*cpc/(2*np.pi*alpha_c*rhoc*G)) # scale height [m]
kappa_c = kc/(cpc*rhoc) #thermal diffusivity of core material
gc = 4/3*np.pi*rc*rhoc*G #gravitational field strength at CMB

# df/dt = D^2/(2*Tc*rc^2*f*(Delta-1)dTc/dt=B dTc/dt) Equation 7 in N09
B=D**2/(2*rc**2*(Delta-1))
Acmb=4*np.pi*rc**2