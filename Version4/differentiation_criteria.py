#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Investigate affect of stagnant lid criterion on differentiation time and temperature
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from parameters import Myr

data = pd.read_csv('lid_test.csv',delimiter=',',header=[0],skiprows=[1])

data.loc[data['core_conv']==np.inf,'core_conv'] = 50 #if core didn't start convecting 50 is a lower limit

labels = ['diff_time','Tdiff','peakT','tmax','strat_end','core_conv']

######################## Absolute plots ######################################
yaxis_lab = ['differentiation time/Myr','Temp at differentiation/K','Peak temperature /K','Time of peak temperature /Myr','End of core thermal stratification /Myr','Onset of core convection/Myr']
for val, ylab in zip(labels[:-1],yaxis_lab[:-1]):
    plt.figure(tight_layout=True)
    plt.scatter(data['convect_ratio'],data[val])
    plt.xlabel('$\delta_0$/r')
    plt.ylabel(f'{ylab}')
    plt.savefig(f'Plots/Lid_test/{val}_abs.png')
    
#plot core convection separately
plt.figure()
plt.scatter(data.loc[data['core_conv']!=50,'convect_ratio'],data.loc[data['core_conv']!=50,'core_conv'])
plt.scatter(data.loc[data['core_conv']==50,'convect_ratio'],data.loc[data['core_conv']==50,'core_conv'],label='lower limit')
plt.xlabel('$\delta_0$/r')
plt.ylabel('Onset of core convection/Myr')
plt.legend()
plt.savefig(f'Plots/Lid_test/core_conv_abs.png')
##################### Percentage plots #######################################
#Don't make one for convective core as not all cores started convecting
for val, ylab in zip(labels[:-1],yaxis_lab[:-1]):
    
    plt.figure(tight_layout=True)
    plt.scatter(data['convect_ratio'],(np.sqrt((data[val]-data[val].mean())**2)/data[val].mean()))
    plt.xlabel('$\delta_0$/r')
    plt.ylabel(f'Fractional RMS deviation in \n {ylab}')
    plt.savefig(f'Plots/Lid_test/{val}_rel.png')

######## Differentiation temp as a fraction of peak temp ##################
plt.figure(tight_layout=True)
plt.scatter(data['convect_ratio'],data['Tdiff']/data['peakT'])
plt.xlabel('$\delta_0$/r')
plt.ylabel(f'Differentiation temp/peak temp')
plt.savefig(f'Plots/Lid_test/temp_compare.png')