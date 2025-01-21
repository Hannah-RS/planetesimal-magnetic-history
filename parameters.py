#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parameters for thermal evolution and dynamo generation model.
The reference for each value is commented next to it. If there is no comment it is from
Dodds et. al. (2021) and references within or it is a fundamental constant.
For non-automated runs, change parameters in this file before running solver.py
"""

import numpy as np
import pandas as pd
import sys
from solidus_calc import solidus, liquidus
#Constants
G = 6.67e-11 # gravitational constant [N kg^-2 m^2]
year = 60*60*24*365 #number of seconds in a year
Myr = 1e6*year #number of seconds in Myr
R = 8.31 # gas constant [J /K /mol]
mu0 = 4*np.pi*1e-7 #magnetic permeability of a vacuum [H/m]

#Run parameters
automated = False #do you want to run a series of automated runs
full_save = True #do you want to save temp profiles etc or just summary stats
B_save = False #do you want to save field strengths and Rem
rhoa_var = False #is the undifferentiated density calculated based on inputs - false for Psyche runs
out_interval = 20 #how many times do you want t to be printed in the whole run
save_interval_d = 0.01*Myr # how often do you want each variable to be saved during differentiation
save_interval_t = 0.1*Myr # how often do you want each variable to be saved during thermal evolution
save_interval_mag = 10*Myr #minimum size for detectable gap in dynamo generation

# Parameters that will vary
if automated == True:
    folder = sys.argv[1]
    auto = pd.read_csv(f'{folder}auto_params.csv',skiprows=[1])
    ind = np.where((auto['status']!=1)&(auto['status']!=0)&(auto['status']!=-1))[0][0] #find index of first uncalculated run
    r = auto.loc[ind,'r']
    rcr = auto.loc[ind,'rcr']
    default = auto.loc[ind,'default']
    rcmf = auto.loc[ind,'rcmf']
    eta0 = auto.loc[ind,'eta0']
    beta = auto.loc[ind,'beta']
    w = auto.loc[ind,'w']
    etal = auto.loc[ind,'etal']
    alpha_n = auto.loc[ind,'alpha_n']
    Xs_0 = auto.loc[ind,'Xs_0']
    Fe0 = auto.loc[ind,'Fe0']
    run = int(auto.loc[ind,'run'])
    t_start_m = auto.loc[ind,'t_start_m']
    t_end_m = auto.loc[ind,'t_end_m']
    dr = auto.loc[ind,'dr']
    icfrac = auto.loc[ind,'icfrac']
    xwater = auto.loc[ind,'xwater']
    accrete = auto.loc[ind,'accrete']
else: #set manually
    r = 100e3 # radius of asteroid [m]
    rcr = 0.5 #core radius as a fraction of asteroid radius
    dr = 500 # grid size [m]
    default ='vary' #default viscosity model
    rcmf = 0.3 #rheologically critical melt fraction - melting required for differentiation
    eta0 = 1e19 #reference viscosity at Tms [Pas]
    beta =0.0225 #E/RTref^2
    w = 5 #width of log linear region [K]
    etal = 10 #liquid viscsoity [Pas]
    alpha_n = 30 #melt weakening (diffusion creep)
    Xs_0 = 30# initial wt % sulfur in core 
    Fe0 = 1e-8 # 60Fe/56FE ratio in accreting material (Dodds 1e-7) (6e-7 Cook 2021)
    run = 23
    t_start_m = 2 #start time - accretion time if accrete = True, differentiation time if accrete = False [Myr]
    t_end_m = 1000 # max end time [Myr]
    icfrac = 0 #fraction of solidified material that forms a passive inner core during solidification
    xwater = 0.05 #water content of mantle [wt %]
    accrete = False #do you want to include accretion to differentiation

# Size of body
#rc = rcr*r #radius of core [m]
n_cells = int(r/dr) +1 #number of cells needed to span the body including one at the centre
nccells = round((n_cells-3)*(rcr))+2 #number of cells needed to span core (inc. centre and CMB)
nmcells = n_cells - nccells +1 #number of cells needed to span mantle plus one extra for CMB
rc = (nccells-1)*dr #radius of core [m], subtract one for centre
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
rhom = 3000 # density [kg m^-3] Bryson et. al. (2019)
#Tml = 1800 # mantle liquidus [K]
#Tms = 1400 # mantle solidus [K]
Tms = solidus(r,rc,xwater,Xs_0,rhom) # mantle solidus [K] Katz et al. 2003
Tml = liquidus(r,rc,xwater,Xs_0,rhom) # mantle liquidus [K] Katz et al. 2003
km = 2.16 # thermal conductivity of silicate [W /m /K] needs to be consistent with cp, kappa and rhom
kappa = 9e-7 # thermal diffusivity of silicate [m^2 /s] - 4 options are given in Dodds (2021) have picked the middle one - needs to be consistent with other choices
alpha_m = 4e-5 # thermal expansivity of mantle [/K]
conv_tol = 0.9 #convective tolerance Ra/Ra_crit<conv_tol = no convection
Rac = 1000  #critical Rayleigh number for isoviscous convection
frht = beta #viscous temperature scale - consistent with Deschamps & Villela (2021)

# Viscosity parameters

#old viscosity models - some parameters here are needed for viscosity comparison code
# Bryson 2019 law (all values from Bryson 2019)
eta0_50 = 1e14 #viscosity of material at 50% melting [Pas]
T0eta = 1400 # reference temperature [K]
Tm50 = 1600 # 50% melting temperature [K]
#Arrhenius model
eta_r50 = 9.358045457250838e+16 #viscosity of material at 50% melting [Pas] for Arrhenius model
#Sterenborg & Crowley (2013)
Tcrit = Tms+(Tml-Tms)/2 #50% melting temperature (Sterenborg and Crowley 2013
E = 300e3 # activation energy [J /mol]
Tref = 1800 # viscosity reference temperature [K] 
c1 = 8 # constant in boundary layer thickness (Sterenborg & Crowley, 2013)
gamma = E/(R*T0eta**2)
Dw = 1e-3 #C^sol_H/C^melt_H partition coefficient of water between solid and melt
mmr = 140 #mineral Mr used for converting ppm to H/1e6 Si
#my model - also uses alpha_n from above
Trcmf = rcmf*(Tml-Tms)+Tms #temperature at critical melt fraction
if frht == 'old':
    frht = (E/(R*Tref**2))
    
#Undifferentiated parameters
#Authors tend to use same value as silicates
ka = km # [W /m /K] same as mantle (correct for post sintering - Dodds 2021)
cpa = cpm # heat capacity [J /kg /K] (Elkins-Tanton 2011 and Bryson 2019 use silicate value)
alpha_a = alpha_m # thermal expansivity of mantle [/K]
convect_ratio = 0.99 #ratio of d0/r for onset of convection in undifferentiated body
tdiff_max = 5*Myr #maximum time for differentiation to occur

#radiogenic heating
h0Al = 0.355 # [W/ kg] heating rate of Al^26 at t=0 (Dodds 2021)
Al0 = 5e-5 # 26Al/27Al ratio in accreting material (Dodds 2021)
XAl_a = 0.014 # abundance of Al in accreting material [wt % /100] (Dodds 2021)
thalf_al = 0.717*Myr # half life of Al26 [s] (Dodds 2021)
h0Fe = 0.0366 # [W/ kg] heating rate of 60Fe at t=0 (Dodds thesis)
thalf_fe = 2.62*Myr # half life of 60Fe [s] (Dodds thesis - Ruedas 2017)

#core and mantle compostion
Mr_s = 32.07 # Molar mass of sulfur Pub chem [amu]
Mr_fe = 55.84 # Molar mass of iron Pub chem  [amu]
XFe_d = 1 - Xs_0/100 #abundance of Fe in core assuming S is only other significant phase [wt %]


# Core parameters 
from fe_fes_liquidus import fe_fes_liquidus_bw, fe_fes_liquidus_bw_min,\
    fe_fes_density, central_pressure
from scipy.optimize import root_scalar
Ts_fe = 1260 # [K] Buono and Walker 2011 give 1263+-25
# Use root finding on FeFeS liquidus for 300km body pressure, 30.05 wt% core
Xs_eutectic = 33 # np.round(root_scalar(fe_fes_liquidus_bw_min,bracket=(25,34),args=(0.154)).root,1) 
cpc=850 # heat capacity of core [J kg^-1 K^-1] 
kc=30 # thermal conducitivity of core material [Wm^-1 K^-1] 
alpha_c=9.2e-5 #[K^-1] Nimmo 2009
Lc = 270e3 # latent heat of core [J kg^-1] Bryson 2015, Tarduno 2012
eta_c =0.01 # viscosity of core [Pa s] Dodds (2021)
f0=0.999999 #initial fractional inner liquid core radius
#constants for magnetic Reynolds number
lambda_mag = 1.3 # magnetic diffusivity [m^2 s^{-1}] Weiss 2010
Omega = 1.75e-4 # rotational frequency (10 hr period) [s^{-1}]
Tmor = 1900 # temperature of measurements in Morard 2019 [K]
rho_exp = 1 + alpha_c*(Tmor-1600) #expansion correction as core not at 1900K
rhofe_l = fe_fes_density(0)*rho_exp # density of pure liquid iron [kg m^-3] Morard (2019)
rhofe_s = 7800 # density of pure solid iron [kg m^-3] Bryson 2015
rho_eut = fe_fes_density(Xs_eutectic)*rho_exp # density of eutectic Fe-FeS solid [kg m^-3] Morard (2019)
rhoc = fe_fes_density(Xs_0)*rho_exp # density of core [kg m^-3]
fohm = 1 #fraction of energy dissipated via Ohmic dissipation in the dynamo (Weiss 2010)
cu = 1.31 #  Aubert 2009
cb = 0.23 # Davies et. al. 2022 median value of c for Bdip,cmb
temp_tol = 1e-8 #minimum usable temp difference

#Calculated parameters
kappa_c = kc/(cpc*rhoc) #thermal diffusivity of core material
g = G*(Vm*rhom+4/3*np.pi*rc**3*rhoc)/r**2 # surface gravity [m/s^2]
gc = 4/3*np.pi*rc*rhoc*G #gravitational field strength at CMB [m/s^2]
Pc = central_pressure(rhom,rhoc,r,rc) #central pressure [GPa]
Tl_fe = fe_fes_liquidus_bw(Xs_0,Pc)
t_cond_core = dr**2/kappa_c #conductive timestep for core
t_cond_mantle = dr**2/kappa #conductive timestep for mantle

if rhoa_var == True: #fixed iron abundance, compatible rc, rhoa responds
    XFe_a = 0.224 # abundance of Fe in accreting material [wt % /100] (Dodds thesis - Lodders 2021 for CV chondrites)
    rhoa = (rhoc*rc**3 + rhom*(r**3-rc**3))/r**3# kg m^-3 density of undifferentiated material (i.e bulk density)    
else: #fix bulk density, compatible rc, accreted iron abundance responds
    rhoa = 4000 # undifferentiated density [kg m^-3]
    XFe_a = rcr**3*(rhoc/rhoa)*XFe_d
    
XAl_d = (rhoa/rhom*(r**3/(r**3-rc**3)))*XAl_a

if automated == True:
    step_m = auto.loc[ind,'dt']*t_cond_core
else:
    step_m = 0.075*t_cond_core
    
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

