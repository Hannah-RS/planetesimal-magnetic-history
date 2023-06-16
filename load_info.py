#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
from parameters import Myr

def load_run_info(run,file,conduction=False):
    """
    Load run info for plotting

    Parameters
    ----------
    run : int
        run number
    file : string
        path to file where run data is stored
    conduction : bool
        whether the mantle switched to conduction

    Returns
    -------
    r : float
        radius of asteroid [m]
    dr : float
        grid spacing [m]
    tstart : float
        start time of simulation [Myr after CAIs]
    tend : float
        max possible time of simulation [Myr]
    tstep : float
        step size in integration [Myr]
    tsolid : float
        time at which solidifcation finishes [Myr]
    dt : float
        variable output frequency [Myr]
    viscosity : string
        viscosity profile

    """
    run_info = pd.read_csv(file,delimiter=',',skiprows=[1])
    r = run_info[run_info['run']==run]['r'].to_numpy()[0] #radius [m]
    dr = run_info[run_info['run']==run]['dr'].to_numpy()[0] #grid spacing [m]
    tstart = run_info[run_info['run']==run]['t_acc'].to_numpy()[0] #accretion time [Myr after CAIs]
    tend = run_info[run_info['run']==run]['t_end_max'].to_numpy()[0] # max possible time of simulation [Myr]
    tstep = run_info[run_info['run']==run]['step_size'].to_numpy()[0] #step size in integration [Myr]
    tsolid = run_info[run_info['run']==run]['t_end_actual'].to_numpy()[0] #time at which solidifcation finishes [Myr] if tsolid == tend then core may not have finished solidifying
    cond_t = run_info[run_info['run']==run]['cond_t'].to_numpy()[0]
    if conduction == True:
        cond_t = cond_t/Myr #time at which mantle switched to conduction
    dt = run_info[run_info['run']==run]['save_interval_t'].to_numpy()[0] #variable output frequency
    viscosity = run_info[run_info['run']==run]['default'].iloc[0] #viscosity profile
    
    return r, dr, tstart, tend, tstep, tsolid, cond_t, dt, viscosity