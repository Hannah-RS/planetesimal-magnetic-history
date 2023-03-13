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
data.loc[data['strat_end']==np.inf,'strat_end'] = 50 #if core didn't start convecting 50 is a lower limit
data['melt frac'] = (data['Tdiff']-1400)/400

################## Melt frac and differentiation temp ######################
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()
ax1.scatter(data['convect_ratio'],data['Tdiff'],color='cornflowerblue')
ax2.scatter(data['convect_ratio'],data['melt frac'],color='cornflowerblue')
ax1.set_ylabel('Differentiation temp /K')
ax1.set_xlabel('$\delta_0/r$')
ax2.set_ylabel('silicate $\phi$ at differentiation')
#plt.savefig('Plots/Lid_test/phi_proxy.png')

######################## Absolute plots ######################################
yaxis_lab = ['differentiation time/Myr','Peak temperature /K','Time of peak temperature /Myr','Onset of erosion of \n core thermal stratification /Myr','End of erosion of \n core thermal stratification /Myr','$F_{CMB} > F_{ad}$','Temp at differentiation/K']
for val, ylab in zip(data.columns[1:-2],yaxis_lab[:-1]):
    xaxis = 'melt frac' #options are 'melt frac' or 'convect_ratio'
    plt.figure(tight_layout=True)
    if val =='core_conv':
       #plot Fcmb > Fad
       plt.scatter(data.loc[data['core_conv']!=50,xaxis],data.loc[data['core_conv']!=50,'core_conv'])
       plt.scatter(data.loc[data['core_conv']==50,xaxis],data.loc[data['core_conv']==50,'core_conv'],label='lower limit')
       plt.legend() 
    elif val == 'strat_end':
        plt.scatter(data.loc[data['strat_end']!=50,xaxis],data.loc[data['strat_end']!=50,'strat_end'])
        plt.scatter(data.loc[data['strat_end']==50,xaxis],data.loc[data['strat_end']==50,'strat_end'],label='lower limit')
        plt.legend()
    else:
        plt.scatter(data[xaxis],data[val])
    plt.ylabel(f'{ylab}')
        
    if xaxis == 'melt frac':
        plt.xlabel('$\phi$ at differentiation')
        plt.savefig(f'Plots/Lid_test/{val}_abs_phi.png')
    else:
        plt.xlabel('$\delta_0$/r')
        #plt.savefig(f'Plots/Lid_test/{val}_abs_dr.png')

    

#plt.savefig(f'Plots/Lid_test/core_conv_abs.png')
##################### Percentage plots #######################################
#Don't make one for convective core as not all cores started convecting
for val, ylab in zip(data.columns[1:-3],yaxis_lab[:-2]):
    
    plt.figure(tight_layout=True)
    plt.scatter(data['convect_ratio'],(np.sqrt((data[val]-data[val].mean())**2)/data[val].mean()))
    plt.xlabel('$\delta_0$/r')
    plt.ylabel(f'Fractional RMS deviation in \n {ylab}')
    #plt.savefig(f'Plots/Lid_test/{val}_rel.png')

######## Differentiation temp as a fraction of peak temp ##################
plt.figure(tight_layout=True)
plt.scatter(data['convect_ratio'],data['Tdiff']/data['peakT'])
plt.xlabel('$\delta_0$/r')
plt.ylabel(f'Differentiation temp/peak temp')
#plt.savefig(f'Plots/Lid_test/temp_compare.png')

################ Put Fcmb > Fad and min_unstable ==0 on same axes
plt.figure(tight_layout=True)
plt.scatter(data.loc[data['core_conv']!=50,'convect_ratio'],data.loc[data['core_conv']!=50,'core_conv'],marker='o',color='black',label='Super adiabatic heat flux')
plt.scatter(data.loc[data['core_conv']==50,'convect_ratio'],data.loc[data['core_conv']==50,'core_conv'],marker='^',color='black',label='lower limit')
plt.scatter(data.loc[data['strat_end']!=50,'convect_ratio'],data.loc[data['strat_end']!=50,'strat_end'],marker='o',color='blue',label='Erosion of stratifciation')
plt.scatter(data.loc[data['strat_end']==50,'convect_ratio'],data.loc[data['strat_end']==50,'strat_end'],marker='^',color='blue',label='lower limit')
plt.legend(loc='upper left')
plt.xlabel('$\delta_0$/r')
plt.ylabel('Time/Myr')
#plt.savefig('Plots/Lid_test/Onset_of_convection_dr.png')

plt.figure(tight_layout=True)
plt.scatter(data.loc[data['core_conv']!=50,'melt frac'],data.loc[data['core_conv']!=50,'core_conv'],marker='o',color='black',label='Super adiabatic heat flux')
plt.scatter(data.loc[data['core_conv']==50,'melt frac'],data.loc[data['core_conv']==50,'core_conv'],marker='^',color='black',label='lower limit')
plt.scatter(data.loc[data['strat_end']!=50,'melt frac'],data.loc[data['strat_end']!=50,'strat_end'],marker='o',color='blue',label='Erosion of stratifciation')
plt.scatter(data.loc[data['strat_end']==50,'melt frac'],data.loc[data['strat_end']==50,'strat_end'],marker='^',color='blue',label='lower limit')
plt.legend(loc='upper right')
plt.ylabel('Time /Myr')
plt.xlabel('Melt fraction at differentiation')
#plt.savefig('Plots/Lid_test/Onset_of_convection_phi.png')
