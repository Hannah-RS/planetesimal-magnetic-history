#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate  magnetic Reynolds numbers
"""
def Rem_therm(Fdrive):
    """
    thermally drived magnetic Reynolds number (equation 16 from Bryson et al. 2019)

    Parameters
    ----------
    Fdrive: float
        FCMB-Fad , heat flux available for driving convection

    Returns
    -------
    Magnetic reynolds number for a thermally driven dynamo

    """
    import numpy as np
    
    #calculate utherm (eqn 17 in Bryson 2019)
    from parameters import G, alpha_c, rc, cpc, Omega, lambda_mag
    
    utherm = ((2*np.pi*G*alpha_c*rc*Fdrive)/(cpc*Omega))**0.5
    
    return utherm*rc/lambda_mag

def Rem_comp(t,f):
    """
    compositionally drived magnetic Reynolds number 
    Rem= ucomp*l/lambda
    (ucomp is equation 1 from Nimmo 2009 which is equivalent to equation 16 from Bryson et al. 2019 but am using Fb = delta rho dri/dt as not
     solidifying at the eutectic)

    Parameters
    ----------
    t: float
        time
    f: float
        fractional inner core radius

    Returns
    -------
    Magnetic reynolds number for a compositionally driven dynamo

    """
    import numpy as np
    from parameters import G,  rc, drho, rhoc, Omega, lambda_mag, f0 #drho here is delta_rho/rho
    
    #calculate dfdt
    "two options here - use analytical expressions in already written functions or numerically take an array of t and f and approximate"
    #calculate numerically
    dfdt = np.gradient(f,t) #this might do weird things when f is constant so use a mask to replace those values
    fmask = np.ma.masked_values(f== f0,f,True)
    dfdt[fmask]=0 #set values with initial inner core size = 0 as core is not growing
    
    #calculate ucomp (eqn 1 in Nimmo 2009 Fb = delta rho dri/dt)
    
    
    ucomp = 0.85*Omega**(-1/5)*rc*(1-f)**(1/5)*(4/3*np.pi*G*drho*rhoc*dfdt)**(2/5)
    
    
    Re_c = ucomp*rc/lambda_mag
    
    
    return Re_c