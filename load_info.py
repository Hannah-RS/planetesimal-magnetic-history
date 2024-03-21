#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd

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
    rcr : float
        ratio of core radius to asteroid radius
    dr : float
        grid spacing [m]
    tstart : float
        start time of simulation [Myr after CAIs]
    tend : float
        max possible time of simulation [Myr]
    tstep : float
        step size in integration [Myr]
    viscosity : string
        viscosity profile

    """
    run_info = pd.read_csv(file,delimiter=',',skiprows=[1])  
    r = run_info[run_info['run']==run]['r'].to_numpy()[0] #radius [m]
    rcr = run_info[run_info['run']==run]['rcr'].to_numpy()[0] #core radius/asteroid radius
    tstep = run_info[run_info['run']==run]['dt'].to_numpy()[0] #step size in integration [Myr]
    dr = run_info[run_info['run']==run]['dr'].to_numpy()[0] #grid spacing [m]
    tstart = run_info[run_info['run']==run]['t_acc_m'].to_numpy()[0] #accretion time [Myr after CAIs]
    viscosity = run_info[run_info['run']==run]['default'].iloc[0] #viscosity profile
    icfrac = run_info[run_info['run']==run]['icfrac'].to_numpy()[0] #fraction of core mass in centre
    return r, rcr, dr, tstart, tstep, viscosity, icfrac

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
    run_results.reset_index(inplace=True,drop=True) #reset indices of sub data frame
    return run_results

def combine_info(folder,params,results,save=False):
    """
    Combine results from runs with their parameters

    Parameters
    ----------
    folder : str
        relative location of output files
    params : str
        name of parameters file
    results : str
        name of results file
    save : bool
        should the dataframe be saved to file, default False

    Returns
    -------
    data : dataframe
        run parameters and results

    """
    indata = pd.read_csv(folder+params,delimiter=',',skiprows=[1],header=0) #import run parameters
    outdata = pd.read_csv(folder+results,delimiter=',',skiprows=[1],header=0) #import run results
    data = pd.merge(indata, outdata, on="run") #join together on matching runs
    if save == True:
        data.to_csv(folder+'all_results.csv',index=False)
    return data
