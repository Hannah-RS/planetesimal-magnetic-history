#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for the IIE workflow
"""
import numpy as np
import pandas as pd

def dynamo_check(tdatalow,tdataup,t,Rem,Remc,tsolid_start):
    """
    Check if the dynamo is on within the range of each Ar-Ar radiometric date
    and whether the core is solidifying at each time
    Parameters
    ----------
    tdatalow : array
        Radiometric dates (lower bound) [Myr]
    tdataup : array
        Radiometric dates (upper bound) [Myr]
    t : array
        Model times [Myr]
    Rem : array
        Rolling-averaged Magnetic Reynolds number
    Remc : float
        Critical magnetic Reynolds number
    tsolid_start : float
        Onset of core solidification [Myr]
    Returns
    -------
    f : bool
        True if the dynamo is on at every radiometric date
    c : list
        List of booleans for each radiometric date,
        True if the core is solidifying at that time
    """
    c = []
    fcheck = [] #counter for dynamo on/off
    for tlow, tup in zip(tdatalow,tdataup):
        if tup < tsolid_start: #check if core is solidifying at each time
            c.append(False)
        else:
            c.append(True)
        if np.any(t>=tup)==True: #check if time is before core fully solidified
            if np.any(Rem[(t>=tlow)&(t<=tup)]> Remc): #check if dynamo is on in time period
                fcheck.append(True)
            else:
                fcheck.append(False)
        else:
            fcheck.append(False)
    if all(fcheck) == True:
        f = True
    else:
        f = False
    return f, c

def paleo_check(tdatalow,tdataup,t,B,Brel,Brel_err):
    """
    Check if the relative paleointensity of the model data is within the range
    of the experimental data
    Parameters
    ----------
    tdatalow : array
        Radiometric dates (lower bounds) [Myr]
    tdataup : array
        Radiometric dates (upper bounds) [Myr]
    t : array
        Model times [Myr]
    B : array
        Rolling-averaged magnetic field strength [T]
    Brel : array
        Relative paleointensities of the experimental data
    Brel_err : array
        Uncertainties in the relative paleointensities of the experimental data
    Returns
    -------
    f : bool
        True if the relative paleointensity of the model data is within the 
        range of the experimental data
    """
    #find model B for each radiometric date
    Bmodel = np.zeros([2,3])
    i = 0
    for tlow, tup in zip(tdatalow,tdataup):
        Bmodel[0,i] = B[t>=tlow][0]
        Bmodel[1,i] = B[t>=tup][0]
        i += 1
    #normalise using max average B
    Brelmodel = Bmodel/np.max((Bmodel[:,0]+Bmodel[:,1])/2)
    fcheck = [] #counter for relative paleointensity test
    i = 0
    for Brellow, Brelup in zip(Brelmodel[:,0],Brelmodel[:,1]):
        #one of values must lie in range
        if ((Brellow >= Brel[i]-Brel_err[i]) & (Brellow <= Brel[i]+Brel_err[i])) \
            | ((Brelup >= Brel[i]-Brel_err[i]) & (Brelup <= Brel[i]+Brel_err[i])):
            fcheck.append(True)
        else:
            fcheck.append(False)
    if all(fcheck) == True:
        f = True
    else:
        f = False
    return f