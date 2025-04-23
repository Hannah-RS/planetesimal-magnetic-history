#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Dynamo functions for the pallasite workflow
"""
import numpy as np
import pandas as pd

def dynamo_check(Bdat,rind,t,temp,Rem,Remc,B,tsolid_start,rplot):
    """
    Check if the dynamo is on when each position cools through 593K
    and whether the core is solidifying at each time
    Parameters
    ----------
    Bdat : array
        Magnetic field strength from data [muT]
    rind : array
        Indices within the mantle for each meteorite
    t : array
        Model output time [Myr]
    temp : array
        Model mantle temperature profile [K]
    Rem : array
        Rolling-averaged Magnetic Reynolds number
    Remc : float
        Critical magnetic Reynolds number
    B : array
        Rolling-averaged magnetic field strength [muT]
    tsolid_start : float
        Onset of core solidification [Myr]
    rplot : array
        Radial grid in mantle [km]
    Returns
    -------
    f : bool
        True if the dynamo is on for the correct meteorites
    depth : array
        Depths for each meteorite (upper and lower bounds) [km]
    tdata : array
        Date for each magnetic remanence (upper and lower bounds) [Myr]
    Bmodel : array
        Upper and lower bounds on rolling-averaged magnetic field strength for each meteorite [T]
    c : list
        List of booleans for each radiometric date,
        True if the core is solidifying at that time
    """
    c = [] #core solidification check
    fcheck = np.zeros([5]) #counter for dynamo on/off
    depth = np.zeros([5,2]) #depths for each meteorite (upper and lower bounds) [km]
    tdata = np.zeros([5,2]) #dates for each magnetic remanence (upper and lower bounds) [Myr]
    Bmodel = np.zeros([5,2]) #upper and lower bounds on rolling-averaged magnetic field strength for each meteorite [T]
    for i, rval in enumerate(rind[:,0]):
        frcheck = np.zeros([2]) #counter for either time for remanance working
        if np.any(temp[:,rval]<=593): #check if cools below 593K
            tval = np.where(temp[:,rval]<=593)[0][0] #find first time cool below 593K
            #save depths, times, field strength
            depth[i,0] = rplot[-1] - rplot[rind[i,0]] 
            tdata[i,0] = t[np.where(temp[:,rind[i,0]]<=593)[0][0]]
            Bmodel[i,0] = B[tval] #save model B
            if ((Rem[tval] >= Remc) & (Bdat[i] > 0)) | ((Rem[tval] < Remc) & (Bdat[i] <=0)): 
                #dynamo is on/off and should be
                frcheck[0] = 1
                
        #do for upper bound
        if np.any(temp[:,rind[i,1]]<=593): #check if cools below 593K
            tval = np.where(temp[:,rind[i,1]]<=593)[0][0] #find first time cool below 593K
            #save depths, times, field strength
            depth[i,1] = rplot[-1] - rplot[rind[i,1]]
            tdata[i,1] = t[np.where(temp[:,rind[i,1]]<=593)[0][0]]
            Bmodel[i,1] = B[tval] #save model B
            if ((Rem[tval] >= Remc) & (Bdat[i] > 0)) | ((Rem[tval] < Remc) & (Bdat[i] <=0)): 
                #dynamo is on/off and should be
                frcheck[1] = 1    
        
        if np.any(frcheck == 1): #if either time works
            fcheck[i] = 1
            if frcheck[0] == 1: #use first value for c, B
                tval = np.where(temp[:,rval]<=593)[0][0] #find first time cool below 593K
            else: #use second value for c, B
                tval = np.where(temp[:,rind[i,1]]<=593)[0][0]
            if t[tval] >= tsolid_start: #check if core is solidifying
                c.append(True)
            else:
                c.append(False)
        else:
            c.append(False)

    if np.all(fcheck==1): #if all meteorites work
        f = True
    else:
        f = False
        #overwrite any saved data
        depth = np.zeros([5,2])
        tdata = np.zeros([5,2])
        Bmodel = np.zeros([5,2])
        
    return f, depth, tdata, Bmodel, c

def paleo_check(Bmodel,Brel,Brel_err):
    """
    Check if the relative paleointensity of the model data is within the range
    of the experimental data
    Parameters
    ----------
    Bmodel : array
        Predicted rolling-averaged magnetic field strength for each meteorite [muT]
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
