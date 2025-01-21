#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from parameters import eta0, gamma, alpha_n, Tms, Tml, eta0_50, Tm50, Tcrit, T0eta, default, xwater
from variable_viscosity import eta_calc

def viscosity(Tm,model = default,Tms=Tms,Tml=Tml,xwater=xwater):
    """
    Overall viscosity function which can call different viscosity models 
    
    Parameters
    ----------
    Tm : float
        mantle temperature [K]
    model: str, default specified in run parameters
        viscosity model to use, options are 'vary' (Sanderson et. al. 2025 model), 'Bryson', 'Dodds'
    Tms : float
        solidus temperature [K], default in parameters file
    Tml : float
        liquidus temperature [K], default in parameters file
    xwater : float
        concentration of water in wt %, default in parameters file
    Returns
    -------
    eta: float
        mantle viscosity [Pas]

    """
        
    if model == 'vary':
        # Sanderson et. al. 2024 model
        eta = eta_calc(Tm,Tms,Tml,xwater)
        
    elif model == 'Dodds':   
    # Viscosity model from Dodds et al (2021)
        log10_eta = 64 - Tm/29 - 5*np.tanh((Tm-1625)/15)
        eta = 10**log10_eta
    elif model == 'Bryson':
        # Viscosity model from Bryson et al. (2019)
        if type(Tm)==np.ndarray: # if an array will need to loop
            n = len(Tm)
            eta = np.zeros([n],dtype='float64')
            for i in range(n):
                if Tm[i] > 1650:
                    eta[i] = 100
                elif Tm[i] <= 1650 and Tm[i] > 1600:
                    eta[i] = eta0_50*np.exp(-gamma*(Tm[i]-Tm50))*np.exp(-alpha_n*(Tm[i]-Tm50)/50)
                elif Tm[i] <= 1600:
                    eta[i] = eta0*np.exp(-gamma*(Tm[i]-T0eta))*np.exp(-alpha_n*(Tm[i]-Tms)/(Tml-Tms))
        else: #just an integer perform once
           if Tm > 1650:
               eta = 100
           elif Tm <= 1650 and Tm > 1600:
               eta = eta0_50*np.exp(-gamma*(Tm-Tm50))*np.exp(-alpha_n*(Tm-Tm50)/50)
           elif Tm <= 1600:
               eta = eta0*np.exp(-gamma*(Tm-T0eta))*np.exp(-alpha_n*(Tm-Tms)/(Tml-Tms))
    
    elif model == 'Sterenborg': 
        #model from Sterenborg & Crowley (2012) 
        # use rounded solidus and liquidus temperatures from Bryson (2019)
        if type(Tm)==np.ndarray: # if an array will need to loop
            n = len(Tm)
            eta = np.zeros([n],dtype='float64')
            for i in range(n):
                if Tm[i] <= Tcrit:
                    eta[i] = eta0*np.exp(-gamma*(Tm[i]-T0eta))*np.exp(-alpha_n*(Tm[i]-Tms)/(Tml-Tms)) 
                else:
                    eta[i] = 1
        else: #if an integer just do once
            if Tm <= Tcrit:
                eta = eta0*np.exp(-gamma*(Tm-T0eta))*np.exp(-alpha_n*(Tm-Tms)/(Tml-Tms)) 
            else:
                eta = 1
    
    else: 
         raise ValueError('Invalid viscosity model chosen')
     
    return eta