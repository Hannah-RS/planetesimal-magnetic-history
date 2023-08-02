#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Timing comparison plot
"""
from scatter_function import make_scatter, make_sub_scatter
import sys
# setting path
sys.path.append('../')

import matplotlib.pyplot as plt
import pandas as pd

folder = 'Fullrun2'
data = pd.read_csv(f'../Results_combined/{folder}/all_sucess_info.csv',delimiter=',',skiprows=[1],header=0,index_col=False)
data['r']=data['r']/1e3 #rescale to km

#filter for the variable you are interested in
var1 = 'r' #variable of interest for these plots
varlab1 = 'radius /km' 
logval1 = False
save = False #do you want to save the figures
#fixed values of other quantities
eta0=1e21
r = 100
rcmf = 0.2
frht=0.011
Xs_0 = 30
Fe0 = 1e-7
alpha_n = 25
etal = 10

#apply sucessive filters and skip chosen variables
if (var1 != 'r'):
    data = data[data['r']==r]
if (var1 != 'Xs_0'):
    data = data[data['Xs_0']==Xs_0]
if (var1 != 'Fe0'):
    data = data[data['Fe0']==Fe0]
if (var1 != 'rcmf'):
    data = data[data['rcmf']==rcmf]
if (var1 != 'frht'):
    data = data[data['frht']==frht]
if (var1 != 'eta0'):
    data = data[data['eta0']==eta0]
if (var1 != 'etal'):
    data = data[data['etal']==etal]
if (var1 != 'alpha_n'):
    data = data[data['alpha_n']==alpha_n]
#################### Thermal history timings ###############################
#only one variable on y axis, keep all other variables fixed and filter data 
#for final version could try James suggestion and only display a few y axis values in horizontal bars
plt.figure(figsize=[15,5])
if var1 =='frht':
    plt.title(f'Thermal history and dynamo timings \n for r={r}, X$_{{s,0}}$={Xs_0}, rcmf={rcmf}, $\\eta_0$={eta0}, $^{{60}}$Fe/$^{{56}}$Fe ={Fe0} ')
elif var1 == 'rcmf':
    plt.title(f'Thermal history and dynamo timings \n for r={r}, X$_{{s,0}}$={Xs_0}, frht={frht}, $\\eta_0$={eta0}, $^{{60}}$Fe/$^{{56}}$Fe ={Fe0}')
elif var1 =='eta0':
    plt.title(f'Thermal history and dynamo timings \n for r={r}, X$_{{s,0}}$={Xs_0}, frht{frht}, rcmf={rcmf}')
elif var1 =='Fe0':
    plt.title(f'Thermal history and dynamo timings \n for r={r}, X$_{{s,0}}$={Xs_0}, frht{frht}, rcmf={rcmf},$\\eta_0$={eta0},')
plt.scatter(data['diff_time'],data[var1],label='differentiation',marker='s',color='darkred')
plt.scatter(data['tmax'],data[var1],label='peak mantle temp',marker='+',color='firebrick')
plt.scatter(data['tcoremax'],data[var1],label='peak core temp',marker='x',color='rosybrown')
plt.scatter(data['tstrat_remove'],data[var1],label='erosion of core stratification',marker='o',color='indianred')
plt.scatter(data['terode'],data[var1],label='core stratification removed',marker='p',color='lightcoral')
plt.scatter(data['mac_on'],data[var1],label='MAC dynamo start',marker='<',color='darkblue')
plt.scatter(data['mac_off'],data[var1],label='MAC dynamo end',marker='>',color='darkblue')
plt.scatter(data['cia_on'],data[var1],label='CIA dynamo start',marker='3',color='royalblue')
plt.scatter(data['cia_off'],data[var1],label='CIA dynamo end',marker='4',color='royalblue')
plt.scatter(data['fcond_t'],data[var1],label='end of mantle convection',marker='*',color='palevioletred')
plt.scatter(data['comp_on'],data[var1],label='compositional dynamo start',marker='^',color='skyblue')
plt.scatter(data['comp_off'],data[var1],label='compositional dynamo end',marker='v',color='skyblue')
plt.scatter(data['tsolid'],data[var1],label='core solidified',marker='d',color='pink')
plt.xlabel('Time/Myr')
plt.xscale('log')
if logval1 == True:
    plt.yscale('log')
plt.ylabel(varlab1)
plt.legend(ncols=2,bbox_to_anchor=(1,0.9))
if save == True:
    plt.savefig(f'../Plots/{folder}all_timings_{var1}.png',dpi=450,bbox_inches='tight') 
