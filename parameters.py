#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parameters for coupled mantle and core evolution model. 
The source of each number is commented next to it. If there is no comment it is from
Dodds et. al. (2021) and references within or it is a fundamental constant
"""

import numpy as np
import pandas as pd

#Constants
G = 6.67e-11 # gravitational constant [N kg^-2 m^2]
year = 60*60*24*365 #number of seconds in a year
Myr = 1e6*year #number of seconds in Myr
R = 8.31 # gas constant [J /K /mol]
mu0 = 4*np.pi*1e-7 #magnetic permeability of a vacuum [H/m]

#Run parameters
automated = True
full_save = True #do you want to save temp profiles etc or just summary stats
out_interval = 20 #how many times do you want t to be printed in the whole run
save_interval_d = 0.01*Myr # how often do you want each variable to be saved during differentiation
save_interval_t = 0.1*Myr # how often do you want each variable to be saved during thermal evolution

# Parameters that will vary
if automated == True:
    auto = pd.read_csv('auto_params.csv',skiprows=[1])
    ind = np.where((auto['status']!=1)&(auto['status']!=0))[0][0] #find index of first uncalculated run
    r = auto.loc[ind,'r']
    default = auto.loc[ind,'default']
    rcmf = auto.loc[ind,'rcmf']
    Xs_0 = auto.loc[ind,'Xs_0']
    Fe0 = auto.loc[ind,'Fe0']
    run = int(auto.loc[ind,'run'])
    t_acc_m = auto.loc[ind,'t_acc_m']
    t_end_m = auto.loc[ind,'t_end_m']
    dr = auto.loc[ind,'dr']
else: #set manually
    r = 400e3 # radius of asteroid [m]
    dr = 500 # grid size [m]
    default ='Dodds' #default viscosity model
    rcmf = 0.2 #rheologically critical melt fraction - melting required for differentiation
    Xs_0 = 28 # initial wt % sulfur in core 
    Fe0 = 1e-7 # 60Fe/56FE ratio in accreting material (Dodds 1e-7) (6e-7 Cook 2021)
    run = 1
    t_acc_m = 0.8 #accretion time [Myr]
    t_end_m = 10 # max end time [Myr]

# Size of body
rc = r/2 #radius of core [m]
Acmb = 4*np.pi*rc**2 #CMB surface area [m^2]
As = 4*np.pi*r**2 # surface area of asteroid [m^2]
V = 4/3*np.pi*r**3 #volume of asteroid [m^3]
Vc = 4/3*np.pi*rc**3 #volume of core [m^3]
Vm = V - Vc #volume of mantle [m^3]

#Surface parameters
Ts = 200 # surface temperature, modelled as fixed [K]

# Mantle parameters
cpm = 800 # heat capacity [J /kg /K] (Elkins-Tanton 2011)
Lm = 400e3 #latent heat of mantle [J /kg]
Tml = 1800 # mantle liquidus [K]
Tms = 1400 # mantle solidus [K]
rhom = 3000 # density [kg m^-3] Bryson et. al. (2019)
km = 2.16 # thermal conductivity of silicate [W /m /K] needs to be consistent with cp, kappa and rhom
kappa = 9e-7 # thermal diffusivity of silicate [m^2 /s] - 4 options are given in Dodds (2021) have picked the middle one - needs to be consistent with other choices
alpha_m = 4e-5 # thermal expansivity of mantle [/K]

Rac = 1000  #critical Rayleigh number for isoviscous convection


# Viscosity parameters
E = 300e3 # activation energy [J /mol]
Tref = 1800 # viscosity reference temperature [K] 
c1 = 8 # constant in boundary layer thickness (Sterenborg & Crowley, 2013)

#viscosity models - some parameters here are needed for viscosity comparison code
# Bryson 2019 law (all values from Bryson 2019)
eta0 = 1e21 # reference viscosity [Pa s] - assumed constant in this model
eta0_50 = 1e14 #viscosity of material at 50% melting [Pas]
alpha_n = 25 # constant in viscosity model
T0eta = 1400 # reference temperature [K]
Tm50 = 1600 # 50% melting temperature [K]
#Arrhenius model
eta_r50 = 9.358045457250838e+16 #viscosity of material at 50% melting [Pas] for Arrhenius model
#Sterenborg & Crowley (2013)
Tcrit = Tms+(Tml-Tms)/2 #50% melting temperature (Sterenborg and Crowley 2013
gamma = E/(R*T0eta**2)

#Undifferentiated parameters
#Authors tend to use same value as silicates
ka = km # [W /m /K] same as mantle (correct for post sintering - Dodds 2021)
cpa = cpm # heat capacity [J /kg /K] (Elkins-Tanton 2011 and Bryson 2019 use silicate value)
alpha_a = alpha_m # thermal expansivity of mantle [/K]
convect_ratio = 0.99 #ratio of d0/r for onset of convection in undifferentiated body
 
#radiogenic heating
h0Al = 0.355 # [W/ kg] heating rate of Al^26 at t=0 (Dodds 2021)
Al0 = 5e-5 # 26Al/27Al ratio in accreting material (Dodds 2021)
XAl_a = 0.014 # abundance of Al in accreting material [wt % /100] (Dodds 2021)
thalf_al = 0.717*Myr # half life of Al26 [s] (Dodds 2021)
h0Fe = 0.0366 # [W/ kg] heating rate of 60Fe at t=0 (Dodds thesis)
XFe_a = 0.224 # abundance of Fe in accreting material [wt % /100] (Dodds thesis - Lodders 2021 for CV chondrites)
thalf_fe = 2.62*Myr # half life of 60Fe [s] (Dodds thesis - Ruedas 2017)

#core and mantle compostion
Mr_s = 32.07 # Molar mass of sulfur Pub chem [amu]
Mr_fe = 55.84 # Molar mass of iron Pub chem  [amu]
XFe_d = 1 - Xs_0/100 #abundance of Fe in core assuming S is only other significant phase [wt %]


# Core parameters 
from fe_fes_liquidus import fe_fes_liquidus_bw, fe_fes_density
Ts_fe = 1260 # [K] Buono and Walker 2011 give 1263+-25
Xs_eutectic = 32 # Buono & Walker give 31.9-32.7 for 10-500km size bodies
cpc=850 # heat capacity of core [J kg^-1 K^-1] 
kc=30 # thermal conducitivity of core material [Wm^-1 K^-1] 
alpha_c=9.2e-5 #[K^-1] Nimmo 2009
Lc = 270e3 # latent heat of core [J kg^-1] Bryson 2015, Tarduno 2012
eta_c =0.01 # viscosity of core [Pa s] Dodds (2021)
f0=0.999 #initial fractional inner liquid core radius
#constants for magnetic Reynolds number
lambda_mag = 1.3 # magnetic diffusivity [m^2 s^{-1}] Weiss 2010
Omega = 1.75e-4 # rotational frequency (10 hr period) [s^{-1}]
Tmor = 1900 # temperature of measurements in Morard 2019 [K]
rho_exp = 1 + alpha_c*(Tmor-1600) #expansion correction as core not at 1900K
rhofe_l = fe_fes_density(0)*rho_exp # density of pure liquid iron [kg m^-3] Morard (2019)
rhofe_s = 7800 # density of pure solid iron [kg m^-3] Bryson 2015
rho_eut = fe_fes_density(Xs_eutectic)*rho_exp # density of eutectic Fe-FeS solid [kg m^-3] Morard (2019)
rhoc = fe_fes_density(Xs_0)*rho_exp # density of core [kg m^-3]
Bp_frac = 0.1 #fraction of poloidal field at CMB to total field in the core (Weiss 2010)
fohm = 1 #fraction of energy dissipated via Ohmic dissipation in the dynamo (Weiss 2010)
cu = 1.65 #  Aubert 2009
cb = 1.31 # Aubert 2009

#Calculated parameters
rhoa = 1/(XFe_a/rhofe_s +(1-XFe_a)/rhom) # kg m^-3 density of undifferentiated material (Sturtz 2022b eqn. 1)
XAl_d = (rhoa/rhom*(r**3/(r**3-rc**3)))*XAl_a
kappa_c = kc/(cpc*rhoc) #thermal diffusivity of core material
g = G*(Vm*rhom+4/3*np.pi*rc**3*rhoc)/r**2 # surface gravity [m/s^2]
gc = 4/3*np.pi*rc*rhoc*G #gravitational field strength at CMB [m/s^2]
Pc = 2*np.pi*G*(rc**2*rhoc+rhom**2*(r**2-rc**2))/1e9 #pressure at centre of core [GPa]
Tl_fe = fe_fes_liquidus_bw(Xs_0,Pc)
t_cond_core = dr**2/kappa_c #conductive timestep for core
t_cond_mantle = dr**2/kappa #conductive timestep for mantle

if automated == True:
    step_m = auto.loc[ind,'dt']*t_cond_core
else:
    step_m = 0.1*t_cond_core
    
#Modified specific heat capacities
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

