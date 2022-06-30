#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Expression for dTc/dt using Equation 6 and Table 1 from Nimmo, F. (2009). Energetics of asteroid dynamos and the role of compositional convection. 
Can toggle on and off effects of conductivity, radiogenic heating, secular cooling, latent heat. 

dTc/dt=(Ephi-Ek-Er)/(Est+Egt+Elt) where Est, Egt, ELt are the expressions for Es, Eg, El divided by dTc/dt
Written as a function of t, Tc to comply with order expected by solve_ivp

"""
def dTcdt_calc(t,Tc,Ek=True,Er=True,Es=True,Eg=True,El=True):
    """

    Parameters
    ----------
    t : float
        time, s
    Tc : float
        core temperature
    Ek : boolean, optional
        entropy contribution from core conductivity. The default is True. If True the contribution is included
    Er : boolean, optional
        entropy contribution from radiogenic heating. The default is True.
    Es : boolean, optional
        entropy contribution from secular cooling The default is True.
    Eg : boolean, optional
        entropy contribution from compositional buoyancy. The default is True.
    El : boolean, optional
        entropy contribution from latent heat. The default is True.

    Returns
    -------
    dTcdt : float
            rate of change of core temperature (if negative core is cooling)

    """
    #check at least one term containing dTc/dt is non zero
    if Es == False and Eg == False and El == False:
        raise ValueError('At least one source of entropy production due to core cooling must be non-zero')
    else: 
        den = 0 # create a variable which will be the value of the denominator, will be changed later as the criterion has passed

    #import constants and parameters
    from parameters import Tc0, f0, B
    
    #calculate f for current value of Tc
    from f_def import f_calc
    f = f_calc(Tc,Tc0,f0,B)                                                                                                                                         
    
    #run through optional functions importing those required - go through the statements and add to one function
    
    from ephi_func import Ephi #always need this function
    
    num = Ephi(f) #create a variable which has the value of the numerator, assign Ephi
    
    if Ek == True:
        from ek_func import Ek
        num = num + Ek() #add conductive contribution
    else: pass

    if Er == True:
        from er_func import Er
        num = num + Er(Tc) #add radiogenic contribution
    else: pass
    
    
    if Es == True:
        from es_func import Est
        den = den + Est(Tc) #add secular heating contribution
    else: pass
        
    if Eg == True:
        from eg_func import Egt
        den = den + Egt(Tc,f) #add compositional buoyancy contribution
    else: pass   
    
    if El == True:
        from el_func import Elt
        den = den + Elt(Tc,f) #add latent heat contribution
    else: pass

    #now have added all terms calculate dTc/dt
    
    return num/den