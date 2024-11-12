#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function for averaging B and Rem, accoutns for a lot of edge cases so quite complicated
"""
import numpy as np
import pandas as pd

def average_B_rem(B,Rem,t,xs,xs_eut,tsolid_start,Rem_c=10):
    """
    Average B and Rem values post onset of solidification

    Parameters
    ----------
    B : np.array
        Magnetic field strength [muT]
    Rem : np.array
        Magnetic Reynolds number 
    t : np.array
        Time [Myr]
    xs : np.array
        Core sulfur content [wt%]
    xs_eut : float
        Core sulfur content at eutectic [wt%]
    tsolid_start : float
        Time of solidification onset [Myr]
    Rem_c : float, optional
        Critical Rem value. Default is 10.
    
    Returns
    -------
    Bplot : np.array
        B values with averaged B values post onset of solidification [muT]
    Remplot : np.array
        Rem values with averaged Rem values post onset of solidification
    """
    #find when core reaches eutectic composition
    t_eut = t[xs>=xs_eut][0]
    #process B values
    Bplot = B[t<tsolid_start] #before solidification
    Bdf = pd.Series(B[(t>=tsolid_start)&(t<t_eut)]) # convert to pandas series
    Blate = B[t>=t_eut] #after eutectic
    tsolids = t[(t>=tsolid_start)&(t<t_eut)]
    #split into 3 series
    t1 = 10 #first threshold [Myr]
    t2 = 100 # second threshold [Myr]
    wn1 = 10 #window for rolling average - 1Ma
    wn2 = 100 #window for rolling average - 10Ma
    if t_eut < 250:
        wn3 = 250 #window for rolling average - 25Ma
    else:
        wn3 = 500 #window for rolling average - 50Ma
    Bdf_short = Bdf[tsolids < t1]
    Bdf_med = Bdf[(tsolids >=t1)&(tsolids<t2)]
    Bdf_long = Bdf[tsolids>=t2]
    if len(Bdf_short)>1:
        Bshort_av = Bdf_short.rolling(window=wn1,center=False).mean() #calculate rolling average
        if len(Bdf_short) > wn1:
            Badd = np.ones([wn1])*np.average(Bdf_short[:wn1]) #average initial window width
            Bplot = np.concatenate([Bplot,Badd,Bshort_av.values[wn1:]])
        else:
            Badd = np.ones([len(Bdf_short)])*np.average(Bdf_short) #average initial window width
            Bplot = np.concatenate([Bplot,Badd])           
    if len(Bdf_med)>1:
        if len(Bdf_med) < wn2: #if less than window width
            wn2 = wn1 #use previous window width
        if len(Bdf_short)>wn2: #concatenate part of previous array to cover blank space of rolling average
            Bdf_med = pd.concat([Bdf_short[-wn2:],Bdf_med])
            Bmed_av = Bdf_med.rolling(window=wn2,center=False).mean() #calculate rolling average
            Bplot = np.concatenate([Bplot,Bmed_av.values[wn2:]])
        else: # previous array too small to concatenate
            Bmed_av = Bdf_med.rolling(window=wn2,center=False).mean()
            Badd = np.ones([wn2])*np.average(Bdf_med[:wn2]) #average initial window width
            Bplot = np.concatenate([Bplot,Badd,Bmed_av.values[wn2:]])
    if len(Bdf_long)>1:
        if len(Bdf_long) < wn2: #if less than medium window width
            Badd = np.ones([len(Bdf_long)])*np.average(Bdf_long) #average initial window width
            Bplot = np.concatenate([Bplot,Badd])
        else: #enough to calculate rolling average
            if len(Bdf_long) < wn3: #if less than window width
                    wn3 = wn2
            if len(Bdf_med)> wn3: #concatenate part of previous array to cover blank space of rolling average
                Bdf_long = pd.concat([Bdf_med[-wn3:],Bdf_long])
                Blong_av = Bdf_long.rolling(window=wn3,center=False).mean() #calculate rolling average
                Bplot = np.concatenate([Bplot,Blong_av.values[wn3:]])
            else: # previous array too short to concatenate    
                Blong_av = Bdf_long.rolling(window=wn3,center=False).mean()
                Badd = np.ones([wn3])*np.average(Bdf_long[:wn3]) #average initial window width
                Bplot = np.concatenate([Bplot,Badd,Blong_av.values[wn3:]])  
    #process Rem values
    Remplot = Rem[t<tsolid_start] #before solidification
    Remdf = pd.Series(Rem[(t>=tsolid_start)&(t<t_eut)]) # convert to pandas series
    Remlate = Rem[t>=t_eut] #after eutectic
    Remdf_short = Remdf[tsolids < t1]
    Remdf_med = Remdf[(tsolids >=t1)&(tsolids<t2)]
    Remdf_long = Remdf[tsolids>=t2]
    if len(Remdf_short)>1:
        Remshort_av = Remdf_short.rolling(window=wn1,center=False).mean() #calculate rolling average
        if len(Remdf_short)>wn1:
            Remadd = np.ones([wn1])*np.average(Remdf_short[:wn1]) #average initial window width
            Remplot = np.concatenate([Remplot,Remadd,Remshort_av.values[wn1:]])
        else: 
            Remadd = np.ones([len(Remdf_short)])*np.average(Remdf_short) #average initial window width
            Remplot = np.concatenate([Remplot,Remadd])
    if len(Remdf_med)>1:
        if len(Remdf_med) < wn2: #if less than window width
            wn2 = wn1
        if len(Remdf_short)>wn2: #concatenate part of previous array to cover blank space of rolling average
            Remdf_med = pd.concat([Remdf_short[-wn2:],Remdf_med])
            Remmed_av = Remdf_med.rolling(window=wn2,center=False).mean()
            Remplot = np.concatenate([Remplot,Remmed_av.values[wn2:]])
        else: # previous array too small to concatenate
            Remmed_av = Remdf_med.rolling(window=wn2,center=False).mean()
            Remadd = np.ones([wn2])*np.average(Remdf_med[:wn2])
            Remplot = np.concatenate([Remplot,Remadd,Remmed_av.values[wn2:]])
    if len(Remdf_long)>1:
        if len(Remdf_long) < wn2: #too small for rolling average
            Rem_add = np.ones([len(Remdf_long)])*np.average(Remdf_long)
            Remplot = np.concatenate([Remplot,Rem_add])
        else:
            if len(Remdf_long) < wn3:
                wn3 = wn2
            if len(Remdf_med)>wn3: #concatenate part of previous array to cover blank space of rolling average
                Remdf_long = pd.concat([Remdf_med[-wn3:],Remdf_long])
                Remlong_av = Remdf_long.rolling(window=wn3,center=False).mean()
                Remplot = np.concatenate([Remplot,Remlong_av.values[wn3:]])
            else:   # previous array too small to concatenate
                Remlong_av = Remdf_long.rolling(window=wn3,center=False).mean()
                Remadd = np.ones([wn3])*np.average(Remdf_long[:wn3])
                Remplot = np.concatenate([Remplot,Remadd,Remlong_av.values[wn3:]])    
    #add late values
    Bplot = np.concatenate([Bplot,Blate])
    Remplot = np.concatenate([Remplot,Remlate])
    #filter B by Rem
    Bplot[Remplot<Rem_c]=0

    return Bplot, Remplot

def average_B_Psyche(B,Rem,t,xs,xs_eut,tsolid_start,Rem_c=10):
    """
    Average B and Rem values post onset of solidification

    Parameters
    ----------
    B : np.array
        Magnetic field strength [muT]
    Rem : np.array
        Magnetic Reynolds number 
    t : np.array
        Time [Myr]
    xs : np.array
        Core sulfur content [wt%]
    xs_eut : float
        Core sulfur content at eutectic [wt%]
    tsolid_start : float
        Time of solidification onset [Myr]
    Rem_c : float, optional
        Critical Rem value. Default is 10.
    
    Returns
    -------
    Bplot : np.array
        B values with averaged B values post onset of solidification [muT]
    Remplot : np.array
        Rem values with averaged Rem values post onset of solidification
    """
    #find when core reaches eutectic composition
    t_eut = t[xs>=xs_eut][0]
    #process B values
    Bplot = B[t<tsolid_start] #before solidification
    Bdf = pd.Series(B[(t>=tsolid_start)&(t<t_eut)]) # convert to pandas series
    Blate = B[t>=t_eut] #after eutectic
    tsolids = t[(t>=tsolid_start)&(t<t_eut)]
    #split into 2 series
    t1 = 5 #first threshold [Myr]
    wn1 = 10 #window for rolling average - 1Ma
    wn2 = 50 #window for rolling average - 5Ma
    Bdf_short = Bdf[tsolids < t1]
    Bdf_med = Bdf[(tsolids >=t1)]
    if len(Bdf_short)>1:
        Bshort_av = Bdf_short.rolling(window=wn1,center=False).mean() #calculate rolling average
        if len(Bdf_short) > wn1:
            Badd = np.ones([wn1])*np.average(Bdf_short[:wn1]) #average initial window width
            Bplot = np.concatenate([Bplot,Badd,Bshort_av.values[wn1:]])
        else:
            Badd = np.ones([len(Bdf_short)])*np.average(Bdf_short) #average initial window width
            Bplot = np.concatenate([Bplot,Badd])           
    if len(Bdf_med)>1:
        if len(Bdf_med) < wn2: #if less than window width
            wn2 = wn1 #use previous window width
        if len(Bdf_short)>wn2: #concatenate part of previous array to cover blank space of rolling average
            Bdf_med = pd.concat([Bdf_short[-wn2:],Bdf_med])
            Bmed_av = Bdf_med.rolling(window=wn2,center=False).mean() #calculate rolling average
            Bplot = np.concatenate([Bplot,Bmed_av.values[wn2:]])
        else: # previous array too small to concatenate
            Bmed_av = Bdf_med.rolling(window=wn2,center=False).mean()
            Badd = np.ones([wn2])*np.average(Bdf_med[:wn2]) #average initial window width
            Bplot = np.concatenate([Bplot,Badd,Bmed_av.values[wn2:]])
    
    #process Rem values
    Remplot = Rem[t<tsolid_start] #before solidification
    Remdf = pd.Series(Rem[(t>=tsolid_start)&(t<t_eut)]) # convert to pandas series
    Remlate = Rem[t>=t_eut] #after eutectic
    Remdf_short = Remdf[tsolids < t1]
    Remdf_med = Remdf[(tsolids >=t1)]

    if len(Remdf_short)>1:
        Remshort_av = Remdf_short.rolling(window=wn1,center=False).mean() #calculate rolling average
        if len(Remdf_short)>wn1:
            Remadd = np.ones([wn1])*np.average(Remdf_short[:wn1]) #average initial window width
            Remplot = np.concatenate([Remplot,Remadd,Remshort_av.values[wn1:]])
        else: 
            Remadd = np.ones([len(Remdf_short)])*np.average(Remdf_short) #average initial window width
            Remplot = np.concatenate([Remplot,Remadd])
    if len(Remdf_med)>1:
        if len(Remdf_med) < wn2: #if less than window width
            wn2 = wn1
        if len(Remdf_short)>wn2: #concatenate part of previous array to cover blank space of rolling average
            Remdf_med = pd.concat([Remdf_short[-wn2:],Remdf_med])
            Remmed_av = Remdf_med.rolling(window=wn2,center=False).mean()
            Remplot = np.concatenate([Remplot,Remmed_av.values[wn2:]])
        else: # previous array too small to concatenate
            Remmed_av = Remdf_med.rolling(window=wn2,center=False).mean()
            Remadd = np.ones([wn2])*np.average(Remdf_med[:wn2])
            Remplot = np.concatenate([Remplot,Remadd,Remmed_av.values[wn2:]])
       
    #add late values
    Bplot = np.concatenate([Bplot,Blate])
    Remplot = np.concatenate([Remplot,Remlate])
    #filter B by Rem
    Bplot[Remplot<Rem_c]=0

    return Bplot, Remplot