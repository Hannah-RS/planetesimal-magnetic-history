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
r = 200e3 # radius of asteroid [m]
V = 4/3*np.pi*r**3
rc = r/2 #radius of core [m]
Vc = 4/3*np.pi*rc**3
dr = 500 # size of cells [m]
out_interval =20 #how many times do you want t to be printed in the whole run
save_interval_d = 0.01*Myr # how often do you want each variable to be saved during differentiation
save_interval_t = 0.5*Myr # how often do you want each variable to be saved during thermal evolution

#Surface parameters
Ts=200 # surface temperature, modelled as fixed [K]
As = 4*np.pi*r**2 # surface area [m^2]

#Undifferentiated parameters
ka = 2.295 # [W /m /K] same as mantle (correct for post sintering - Dodds 2021)
rhoa = 3500 # kg m^-3 density of undifferentiated material (Dodds 2021 rho_b)
cpa = 850 # heat capacity [J /kg /K] (Elkins-Tanton 2011, Bryson 2019 use silicate value)
alpha_a = 4e-5 # thermal expansivity of mantle [/K]
# Mantle parameters
cpm = 850 # heat capacity [J /kg /K] (Dodds et al.)
Lm = 400e3 #latent heat of mantle [J /kg]
Tml = 1800 # mantle liquidus [K]
Tms = 1400 # mantle solidus [K]
rhom = 3000 # density [kg m^-3]
alpha_m = 4e-5 # thermal expansivity of mantle [/K]
kappa = 9e-7 # thermal diffusivity of silicate [m^2 /s] - 4 options are given in Bryson (2019) have picked the one used in the figures
Rac = 1000  #critical Rayleigh number 
E = 300e3 # activation energy [J /mol]
R = 8.31 # gas constant [J /K /mol]
Tref = 1800 # viscosity reference temperature (Dodds 2021) [K] 
Tm0 = 1600 # intial mantle temp - from beginning of stage 2 in Figure 1 of Bryson (2019)
c1 = 8 # constant in boundary layer thickness Dodds (2021)
convect_ratio = 0.4 #ratio of d0/r for onset of convection in differentiation

#radiogenic heating
h0Al = 0.355 # [W/ kg] heating rate of Al^26 at t=0 (Dodds 2021)
Al0 = 5e-5 # 26Al/27Al ratio in accreting material (Dodds 2021)
XAl_a = 0.014 # abundance of Al in accreting material [wt % /100] (Dodds 2021)
thalf_al = 0.717*Myr # half life of Al26 [s] (Dodds 2021)
h0Fe = 0.0366 # [W/ kg] heating rate of 60Fe at t=0 (Dodds thesis)
Fe0 = 0#1e-7 # 60Fe/56FE ratio in accreting material (Dodds 1e-7) (6e-7 Cook 2021)
XFe_a = 0.224 # abundance of Fe in accreting material [wt % /100] (Dodds thesis - Lodders 2021 for CV chondrites)
thalf_fe = 2.62*Myr # half life of 60Fe [s] (Dodds thesis - Ruedas 2017)

Xs_0 = 20 # initial wt % sulfur in core (needs a citation)
Mr_s = 32.07 # Pub chem [amu]
Mr_fe = 55.84 #Pub chem  [amu]
XFe_d = 1 - Xs_0/100 #abundance of Fe in core assuming S is only other significant phase [wt %]
Xs_to_core = XFe_a/(1-Xs_0/100)-XFe_a # wt % of S from accreted body that went into core
XAl_d = XAl_a/(1-XFe_d-Xs_to_core)
#viscosity models
default ='Dodds' #default viscosity model
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
Xs_eutectic = 32 # eutectic wt% S


# Initial conditions
Tc0=Tm0# intial CMB temp [K] is same as intial mantle temp from Fig 1 in Bryson (2019)
f0=0.01 #initial fractional inner core radius 

#constants for magnetic Reynolds number
lambda_mag = 1.3 # magnetic diffusivity [m^2 s^{-1}] 2
Omega = 1.75e-4 # rotational frequency (10 hr period) [s^{-1}] 6hr period in N09
beta=0.85 

# Calculated mantle parameters
Vm = 4/3*np.pi*(r**3-rc**3) # volume of magma ocean, assumed to be same as volume of mantle
gamma = E/(R*T0eta**2)
g = G*(Vm*rhom+4/3*np.pi*rc**3*rhoc)/r**2 # surface gravity 

#Calculated core parameters
D= np.sqrt(3*cpc/(2*np.pi*alpha_c*rhoc*G)) # scale height [m]
kappa_c = kc/(cpc*rhoc) #thermal diffusivity of core material
gc = 4/3*np.pi*rc*rhoc*G #gravitational field strength at CMB
Pc = 2*np.pi*G*(rc**2*rhoc+rhom**2*(r**2-rc**2))/1e9 #pressure at centre of core [GPa]

# df/dt = D^2/(2*Tc*rc^2*f*(Delta-1)dTc/dt=B dTc/dt) Equation 7 in N09
B=D**2/(2*rc**2*(Delta-1))
Acmb=4*np.pi*rc**2

#Modified specific heat capacities
from fe_fes_liquidus import fe_fes_liquidus
Tl_fe = fe_fes_liquidus(Xs_0)
Ts_fe = fe_fes_liquidus(Xs_eutectic) # FeFeS solidus [K]
#before differentiation
if Xs_0 != Xs_eutectic:
    cpa_fe = cpa + XFe_a*Lc/(Tl_fe-Ts_fe) #cp for Tfe_s < T < Tsi_s
    cpa_fesi = cpa + XFe_a*Lc/(Tl_fe-Ts_fe) + (1-XFe_a)*Lm/(Tml-Tms) #cp for Tms < T < Tfe_s
else:
    cpa_fe = None
    cpa_fesi = None
    
cpa_si = cpa + (1-XFe_a)*Lm/(Tml-Tms) #cp for Tfe_l < T < Tml
#after differentiation
cpm_p = cpm*(1+(Lm/(cpm*(Tml-Tms)))) # modified mantle heat capacity accounting for latent heat Tms < Tm <Tml

