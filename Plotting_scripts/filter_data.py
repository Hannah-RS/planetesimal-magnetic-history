#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function for filtering the data except for two variable parameters
"""
import pandas as pd

def filter_two(data,var1,var2,r,Xs_0,Fe0,rcmf,frht,etal,eta0,alpha_n):
    """

    Parameters
    ----------
    data : pandas dataframe
        unfiltered data
    var1 : str
        first variable
    var2 : str
        second variable
    r : float
        radius [km]
    Xs_0 : float
        sulfur conten [wt %]
    Fe0 : float
        ^{60}Fe/^{56}Fe
    rcmf : float
        rheologically critical melt fraction
    frht : float
        Arrhenius slope
    etal : float
        liquid viscosity [Pas]
    eta0 : float
        reference viscosity [Pas]
    alpha_n : float
        melt weakening exponent

    Returns
    -------
    data_out : pandas DataFram
        filtered data

    """
    #apply sucessive filters and skip chosen variables
    if (var1 != 'r') & (var2 !='r'):
        data = data[data['r']==r]
    if (var1 != 'Xs_0') & (var2 !='Xs_0'):
        data = data[data['Xs_0']==Xs_0]
    if (var1 != 'Fe0') & (var2 !='Fe0'):
        data = data[data['Fe0']==Fe0]
    if (var1 != 'rcmf') & (var2 !='rcmf'):
        data = data[data['rcmf']==rcmf]
    if (var1 != 'frht') & (var2 !='frht'):
        data = data[data['frht']==frht]
    if (var1 != 'eta0') & (var2 !='eta0'):
        data = data[data['eta0']==eta0]
    if (var1 != 'etal') & (var2 !='etal'):
        data = data[data['etal']==etal]
    if (var1 != 'alpha_n') & (var2 !='alpha_n'):
        data = data[data['alpha_n']==alpha_n]
        
    if len(data)==0:
        raise ValueError('Invalid parameter choice - no data remains')
    return data