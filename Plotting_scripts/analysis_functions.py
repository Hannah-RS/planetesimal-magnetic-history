#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for loading data from runs and combining run info and results
"""
import pandas as pd

def combine_info(folder,params,results,Bfile):
    """
    Combine results from runs with magfield duration data - only works for one on-off period at the moment

    Parameters
    ----------
    folder : str
        relative location of output files
    params : str
        name of parameters file
    results : str
        name of results file
    Bfile : str
        names of magnetic field files

    Returns
    -------
    data : dataframe
        run parameters and results

    """
    indata = pd.read_csv(folder+params,delimiter=',',skiprows=[1]) #import run parameters
    outdata = pd.read_csv(folder+results,delimiter=',',skiprows=[1]) #import run results
    data = pd.merge(indata, outdata, on="run") #join together on matching runs
    
    for file in Bfile:
        magdata = pd.read_csv(folder+file,delimiter=',',skiprows=[1]) #import B field data
    #rename columns
        newname = magdata.columns[2:]+'_'+magdata.loc[0,'label']
        magdata.rename(columns={"start":newname[0],"end":newname[1],"duration":newname[2]},inplace=True)
        #print()
        magdata.drop('label',axis=1,inplace = True) #drop label column
        data = pd.merge(data,magdata,on='run')
    return data

