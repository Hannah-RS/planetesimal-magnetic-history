#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from parameters import eta0, gamma, alpha_n, Tms, Tml, eta0_50, eta_r50, Tm50, E, R, Tcrit, T0eta, default
def viscosity(Tm, model = default):
    """
    Different viscosity models
    
    Parameters
    ----------
    Tm : float
        mantle temperature
    model: str, default Dodds
        viscosity model to use, options are Bryson, Dodds

    Returns
    -------
    Silicate viscosity

    """
    if model == 'Dodds':   
    # Viscosity model from Dodds et al (2021)
        log10_eta = 64 - Tm/29 - 5*np.tanh((Tm-1625)/15)
        eta = 10**log10_eta
        
    elif model == 'Bryson':
        # Viscosity model from Bryson et al. (2019)
        if type(Tm)==np.ndarray: # if an array will need to loop
            n = len(Tm)
            eta = np.zeros([n])
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
    elif model == 'Bryson2':
        # Viscosity model from Bryson et al. (2019)
        if type(Tm)==np.ndarray: # if an array will need to loop
            n = len(Tm)
            eta = np.zeros([n])
            for i in range(n):
                if Tm[i] > 1550:
                    eta[i] = 100
                elif Tm[i] <= 1550 and Tm[i] > 1500:
                    eta[i] = 3e17*np.exp(-gamma*(Tm[i]-1500))*np.exp(-alpha_n*(Tm[i]-1500)/38)
                elif Tm[i] <= 1500:
                    eta[i] = eta0*np.exp(-gamma*(Tm[i]-T0eta))*np.exp(-alpha_n*(Tm[i]-Tms)/(Tml-Tms))

        else: #just an integer perform once
           if Tm > 1550:
               eta = 100
           elif Tm <= 1550 and Tm > 1500:
               eta = 3e17*np.exp(-gamma*(Tm-1500))*np.exp(-alpha_n*(Tm-1500)/38)
           elif Tm <= 1500:
               eta = eta0*np.exp(-gamma*(Tm-T0eta))*np.exp(-alpha_n*(Tm-Tms)/(Tml-Tms))
    elif model == 'Bryson3':
        # Viscosity model from Bryson et al. (2019)
        if type(Tm)==np.ndarray: # if an array will need to loop
            n = len(Tm)
            eta = np.zeros([n])
            for i in range(n):
                if Tm[i] > 1650:
                    eta[i] = 100
                elif Tm[i] <= 1650 and Tm[i] > 1600:
                    eta[i] = 9e6*np.exp(-gamma*(Tm[i]-1600))*np.exp(-alpha_n*(Tm[i]-1600)/120)
                elif Tm[i] <= 1600:
                    eta[i] = 1e14*np.exp(-gamma*(Tm[i]-T0eta))*np.exp(-alpha_n*(Tm[i]-Tms)/(Tml-Tms))
        else: #just an integer perform once
           if Tm > 1650:
               eta = 100
           elif Tm <= 1650 and Tm > 1600:
               eta = 9e6*np.exp(-gamma*(Tm-1600))*np.exp(-alpha_n*(Tm-1600)/120)
           elif Tm <= 1600:
               eta = 1e14*np.exp(-gamma*(Tm-T0eta))*np.exp(-alpha_n*(Tm-Tms)/(Tml-Tms))
    elif model =='Robuchon-Bryson':
        # arrhenius functional form from Robuchon & Nimmo (2011) with constants and three piece from Bryson (2019)

        #use a linear profile for log10(eta) for T>1650 to get same end point as Bryson
        if type(Tm) == np.ndarray: #check if an array before applying the condition
            n = len(Tm)
            eta = np.zeros([n])
            for i in range(n):
                if Tm[i] > 1650:
                    eta[i] = 100
                elif Tm[i] <= 1650 and Tm[i] > 1600:
                    eta[i] = 10**((2-np.log10(eta_r50))*(Tm[i]-1600)/50+np.log10(eta_r50))
                elif Tm[i] <= 1600:
                    eta[i] = eta0*np.exp(-E/R*(1/Tm[i]-1/T0eta))*np.exp(-alpha_n*(Tm[i]-Tms)/(Tml-Tms))
        
        else: 
            if Tm > 1650:
                eta = 100
            elif Tm <= 1650 and Tm > 1600:
                eta = 10**((2-np.log10(eta_r50))*(Tm-1600)/50+np.log10(eta_r50))
            elif Tm <= 1600:
                eta = eta0*np.exp(-E/R*(1/Tm-1/T0eta))*np.exp(-alpha_n*(Tm-Tms)/(Tml-Tms))
                  
        
    
    elif model == 'Sterenborg': #slightly different as use rounded solidus and liquidus temperatures from Bryson (2019)
        if type(Tm)==np.ndarray: # if an array will need to loop
            n = len(Tm)
            eta = np.zeros([n])
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
    
    elif model == 'Arrhenius':   
        # Arrhenius model for viscosity
        eta = eta0*np.exp(E/R*(1/Tm-1/1600))*np.exp(-0.5*45) #semi-based on Robuchon and Nimmo 2011, eqn 8, assuming 50% melt fraction
    
    elif model == 'new':   
    # adapted Viscosity model from Dodds et al (2021)
        log10_eta = 38 - Tm/55 - 6*np.tanh((Tm-1625)/15)
        eta = 10**log10_eta

    else: 
         raise ValueError('Invalid viscosity model chosen')
     
    return eta