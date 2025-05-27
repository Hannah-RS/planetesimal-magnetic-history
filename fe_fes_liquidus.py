#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions required to calculate Fe-FeS liquidus and Fe-FeS density

"""
import numpy as np

def central_pressure(rhom,rhoc,r,rc):
    """
    Pressure at the centre of a planetesimal assuming constant core and mantle density
    Eqn. 31 in Sanderson et. al. (2024)
    
    Parameters
    ----------
    rhom : float
        mantle density [kg m^-3]
    rhoc : float
        core density [kg m^-3]
    r : float
        planetesimal radius [m]
    rc : float
        core radius [m]

    Returns
    -------
    Pc : float
        central pressure [GPa]

    """
    G = 6.67e-11 # gravitational constant [N kg^-2 m^2]
    Pcpa = 2*np.pi*G/3*(rc**2*(rhoc**2-rhom**2)+r**2*rhom**2+2*rhom*(rhoc-rhom)*rc**2*(1-(rc/r))) #[Pa]
    Pc = Pcpa/1e9 #convert to GPa
    return Pc

def weight_perc_to_mole_frac(Xs):
    """
    Convert from wt % sulfur to mole fraction of sulfur.

    Parameters
    ----------
    s : float
        sulfur content [wt %]

    Returns
    -------
    x : float
        mole fraction of S in Fe-FeS

    """
    Mr_s = 32.07 # Pub chem [amu]
    Mr_fe = 55.84 #Pub chem  [amu]
    Xsd = Xs/100 #convert wt % to decimal

    mrr = Mr_fe/Mr_s
    x = Xsd*mrr/(1-Xsd) #mole fraction of FeS in Fe,FeS
    return x

def weight_perc_to_at_frac(Xs):
    """
    Convert from wt % sulfur to atom fraction of sulfur.

    Parameters
    ----------
    Xs : float
        sulfur content [wt %]

    Returns
    -------
    at : float
        atom fraction of S in Fe-FeS

    """
    x = weight_perc_to_mole_frac(Xs)
    at = x/(1+x)
    return at

def fe_fes_density(Xs):
    """
    Fe-FeS density
    Eqn. 37 in Sanderson et. al. (2024)

    Parameters
    ----------
    Xs : float
        sulfur content [wt %]

    Returns
    -------
    rho : float
        density [kg m^-3]
    """
    at = weight_perc_to_at_frac(Xs)
    rho = (-3108*at**2)+(-5176*at)+6950
    return rho

def fe_fes_density_err(Xs,err_Xs):
    """
    Uncertainty on Fe-FeS density based on Eqn. 37 in Sanderson et. al. (2024)

    Parameters
    ----------
    Xs : float
        sulfur content [wt %]
    err_Xs : float
        uncertainty in sulfur content [wt %]

    Returns
    -------
    rho_err : float
        uncertainty in density [kg m^-3]
    """
    xs_at = weight_perc_to_at_frac(Xs)
    err_at = weight_perc_to_at_frac(err_Xs)
    drhodxs = (-3108*2*xs_at-5176)
    rho_err = np.abs(drhodxs*err_at)
    return rho_err

def fe_fes_liquidus_linear(Xs):
    """
    Fe-FeS linear liquidus approximation
    Parameters
    ----------
    Xs : float
        sulfur content [wt %]

    Returns
    -------
    Liquidus temp [K]

    """
    return 1810-18*Xs

def fe_fes_liquidus_bw(Xs,P):
    """
    Fe-FeS liquidus from Buono & Walker (2011) - eqn 29. Valid for sub-eutectic S contents up to 10GPa. 
    Eqn. 29 in Sanderson et. al. (2024)
    
    Parameters
    ----------
    Xs : float
        sulfur content [wt %]
    P : float
        pressure [GPa]       

    Returns
    -------
    Liquidus temp [K]

    """
    x = weight_perc_to_mole_frac(Xs)
    
    A = (-2.4724*P**4) + (28.025*P**3) + (9.1404*P**2) + (581.71*P) + 3394.8
    B = (1.7978*P**4) + (-6.7881*P**3) + (-197.69*P**2) + (-271.69*P) + (-8219.5)
    C = (-0.1702*P**4) + (-9.3959*P**3) + (163.53*P**2) + (-319.35*P) + 5698.6
    D = (-0.2308*P**4) + (7.1*P**3) + (-64.118*P**2) + 105.98*P + (-1621.9)
    E = 0.2302*P**4 + (-5.3688*P**3) + 38.124*P**2 + (-46.681*P) + 1813.8
    
    Tl = A*x**4 + B*x**3 + C*x**2 + D*x**1 + E
    
    return Tl

def fe_fes_liquidus_bw_min(Xs,P):
    """
    Fe-FeS liquidus from Buono & Walker (2011) - eqn 29 
    with eutectic temp subtracted off for minimisation
    Valid for sub-eutectic S contents up to 10GPa. 
    Eutectic temp is 1263+-25K from Buono & Walker

    Parameters
    ----------
    Xs : float
        sulfur content [wt %]
    P : float
        pressure [GPa]       

    Returns
    -------
    Liquidus temp - Eutectic temp [K]

    """
    x = weight_perc_to_mole_frac(Xs)
    
    A = (-2.4724*P**4) + (28.025*P**3) + (9.1404*P**2) + (581.71*P) + 3394.8
    B = (1.7978*P**4) + (-6.7881*P**3) + (-197.69*P**2) + (-271.69*P) + (-8219.5)
    C = (-0.1702*P**4) + (-9.3959*P**3) + (163.53*P**2) + (-319.35*P) + 5698.6
    D = (-0.2308*P**4) + (7.1*P**3) + (-64.118*P**2) + 105.98*P + (-1621.9)
    E = 0.2302*P**4 + (-5.3688*P**3) + 38.124*P**2 + (-46.681*P) + 1813.8
    
    Tl = A*x**4 + B*x**3 + C*x**2 + D*x**1 + E -1260
    
    return Tl

def fe_fes_liquidus_dp(Xs,P):
    """
    dTl/dP using Fe-FeS liquidus from Buono & Walker (2011) - eqn 29. 
    Valid for sub-eutectic S contents up to 10GPa. 

    Parameters
    ----------
    Xs : float
        sulfur content [wt %]
    P : float
        pressure [GPa]       

    Returns
    -------
    dTl/dP 

    """
    x = weight_perc_to_mole_frac(Xs)
       
    F = 4*(-2.4724*P**3)+3*(28.025*P**2)+2*(9.1404*P)+(581.71)
    G = 4*(1.7978*P**3)+3*(-6.7881*P**2)+2*(-197.69*P)+(-271.69)
    H = 4*(-0.1702*P**3)+3*(-9.3959*P**2)+2*(163.53*P)+(-319.35)
    I = (-4*0.2308*P**3) + 3*7.1*P**2 + (-2*64.118*P) + 105.98
    J = 4*0.2302*P**3 + (-3*5.3688*P**2) + 2*38.124*P + (-46.681)
    
    dTl_dPg = F*x**4 + G*x**3 + H*x**2 + I*x**1 + J #in K/GPa
    dTl_dP = dTl_dPg/1e9 #in K/ Pa
    
    return dTl_dP
