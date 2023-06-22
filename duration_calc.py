#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

def on_off_test(tarray,out_array,threshold,save_interval):
    """
    Function for finding start and end periods e.g. for a dynamo

    Parameters
    ----------
    tarray : float
        time series
    out_array : float
        output to test condition on - must be same length as time series
    threshold : float
        threshold value to compare out_array
    save_interval : float
        spacing between time saves

    Returns
    -------
    tstart: float
        start times when out_array > threshold (same units as tarray), uncertainty = save_interval
    tend: float
        end times for when out_array > threshold (same units as tarray), uncertainty = save_interval
    duration : float
        length of time when out_array > threshold, uncertainty = 2*save_interval 
    """
    #check arrays are same length
    assert(len(tarray)==len(out_array))
    t = tarray[out_array>threshold]
    #find difference between sucessive tvalues
    tdiff = np.ediff1d(t)
    tend_ind = np.where(tdiff>save_interval)[0]
    tend = t[tend_ind]
    tstart_ind = np.where(tdiff>save_interval)[0]+1
    tstart = t[tstart_ind]
    #append first and last array values
    tend = np.append(tend,t[-1])
    tstart = np.insert(tstart, 0, t[0])
    duration = tend-tstart
    
    return tstart, tend, duration

#test function
tarray = np.linspace(1,20,20)
Rem = np.array([1,20,20,20,4,4,20,20,20,3,8,11,11,11,11,11,0,0,20,20])
threshold =10
save_interval = 1

print(on_off_test(tarray,Rem,threshold,save_interval))


