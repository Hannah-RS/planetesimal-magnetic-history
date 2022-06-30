#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Expression for dTc/dt using Equation 70 from Vol 8, Treatise on Geophysics (Nimmo 2007) and Table 1 from Nimmo, F. (2009). Qnergetics of asteroid dynamos and the role of compositional convection. 
Can toggle on and off effects of gravitational potential energy release, radiogenic heating, secular cooling, latent heat. 

dTc/dt=(Qcmb-Qr)/(Qst+Qgt+Qlt) where Qst, Qgt, QLt are the expressions for Qs, Qg, Ql divided by dTc/dt
Written as a function of t, Tm, Tc, f to comply with order expected by vectorised argument in solve_ivp

"""
def dTcdt_calc(t,Tm,Tc,f,d0,Qr=True,Qs=True,Qg=True,Ql=True):
    """

    Parameters
    ----------
    t : float
        time, s
    Tm: float
        base of mantle temperature
    Tc : float
        core temperature
    f : float
        fractional inner core radius 
    d0: float
        stagnant lid thickness
    Qr : boolean, optional
        power contribution from radiogenic heating. The default is True.
    Qs : boolean, optional
        power contribution from secular cooling The default is True.
    Qg : boolean, optional
        power contribution from compositional buoyancy. The default is True.
    Ql : boolean, optional
        power contribution from latent heat. The default is True.

    Returns
    -------
    dTcdt : float
            rate of change of core temperature (if negative core is cooling)

    """

    #check at least one term containing dTc/dt is non zero
    if Qs == False and Qg == False and Ql == False:
        raise ValueError('At least one source of power production due to core cooling must be non-zero')
    else: 
        den = 0 # create a variable which will be the value of the denominator, will be changed later as the criterion has passed
        
    #import constants and parameters    
    from parameters import gamma, km, Ts, Acmb
        
    #calculate CMB heat flux
    Fcmb = (Tc-Tm)*gamma*km*(Tm-Ts)/d0
    Qcmb = Acmb*Fcmb
    
    num = Qcmb #create a variable which has the value of the numerator, assign Qcmb


    if Qr == True:
        from qr_func import Qr
        num = num - Qr(t,Tc) #add radiogenic contribution
    else: pass
    
    
    if Qs == True:
        from qs_func import Qst
        den = den + Qst(Tc) #add secular heating contribution
    else: pass
        
    if Qg == True:
        from qg_func import Qgt
        den = den + Qgt(Tc,f) #add compositional buoyancy contribution
    else: pass   
    
    if Ql == True:
        from ql_func import Qlt
        den = den + Qlt(Tc,f) #add latent heat contribution
    else: pass
    
    #now have added all terms calculate dTc/dt  
        
    return num/den