#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
from plotting_constants import Myr

def load_run_info(run,file):
    """
    Load run info for plotting

    Parameters
    ----------
    run : int
        run number
    file : string
        path to file where run data is stored

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
    dt : float
        variable output frequency [Myr]
    viscosity : string
        viscosity profile

    """
    run_info = pd.read_csv(file,delimiter=',',skiprows=[1])  
    r = run_info[run_info['run']==run]['r'].to_numpy()[0] #radius [m]
    tstep = run_info[run_info['run']==run]['step_size'].to_numpy()[0] #step size in integration [Myr]
    dr = run_info[run_info['run']==run]['dr'].to_numpy()[0] #grid spacing [m]
    tstart = run_info[run_info['run']==run]['t_acc'].to_numpy()[0] #accretion time [Myr after CAIs]
    dt = run_info[run_info['run']==run]['save_interval_t'].to_numpy()[0] #variable output frequency
    viscosity = run_info[run_info['run']==run]['default'].iloc[0] #viscosity profile
    
    return r, dr, tstart, tend, tstep, dt, viscosity

def load_run_results(run,file):
    """
    Load run results for plotting

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
    run_results : pandas dataframe
        all calculated results for a given run

    """
    run_info = pd.read_csv(file,delimiter=',',skiprows=[1])  
    run_results = run_info[run_info['run']==run]
    
    return run_results
    return r, dr, tstart, tend, tstep, tsolid, cond_t, dt, viscosity