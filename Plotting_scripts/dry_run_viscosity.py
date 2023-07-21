#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make scatter plots for viscosity test run
"""
from scatter_function import make_scatter, make_sub_scatter
import sys
# setting path
sys.path.append('../')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

folder = 'w_test'
data = pd.read_csv(f'../Results_combined/{folder}/all_sucess_info.csv',delimiter=',',skiprows=[1],header=0,index_col=False)
data['r']=data['r']/1e3 #rescale to km

#filter for the variable you are interested in
var1 = 'r' #variable of interest for these plots
varlab1 = 'radius /km' 
logval1 = False #should the yaxis of this variable be scaled logarithmically
var2 = 'w'
varlab2 = 'w /K'
logval2 = False 
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
###################### Effect on dynamo timing ################################
#### var2 as x axis, var1 as colour
plt.figure(tight_layout=True,figsize=[15,15])
plt.suptitle(f'Effect of {varlab2} on dynamo timing - 1st generation period')
#MAC dynamo
name = 'MAC'
if np.all(np.isnan(data[f'start1_{name}'])==True):
    pass #don't plot this dynamo
else:
    make_sub_scatter(data[var2],data[f'start1_{name}'],varlab2,f'{name} dynamo start /Myr',3,3,1,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
    make_sub_scatter(data[var2],data[f'end1_{name}'],varlab2,f'{name} dynamo end /Myr',3,3,2,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
    make_sub_scatter(data[var2],data[f'duration1_{name}'],varlab2,f'{name}dynamo duration /Myr',3,3,3,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
#CIA dynamo
name = 'CIA'
make_sub_scatter(data[var2],data[f'start1_{name}'],varlab2,f'{name} dynamo start /Myr',3,3,4,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
make_sub_scatter(data[var2],data[f'end1_{name}'],varlab2,f'{name} dynamo end /Myr',3,3,5,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
make_sub_scatter(data[var2],data[f'duration1_{name}'],varlab2,f'{name}dynamo duration /Myr',3,3,6,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
#compositional dynamo
name = 'comp'
make_sub_scatter(data[var2],data[f'start1_{name}'],varlab2,f'{name} dynamo start /Myr',3,3,7,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
make_sub_scatter(data[var2],data[f'end1_{name}'],varlab2,f'{name} dynamo end /Myr',3,3,8,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
make_sub_scatter(data[var2],data[f'duration1_{name}'],varlab2,f'{name}dynamo duration /Myr',3,3,9,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
#if save == True:
plt.savefig(f'../Plots/{folder}dynamo_duration1_{var2}{var1}_col.png')

plt.figure(tight_layout=True,figsize=[15,15])
plt.suptitle(f'Effect of {varlab2} on dynamo timing - second generation period')
#MAC dynamo
name = 'MAC'
if np.all(np.isnan(data[f'start2_{name}'])==True):
    pass #don't plot this dynamo
else:
    make_sub_scatter(data[var2],data[f'start2_{name}'],varlab2,f'{name} dynamo start /Myr',3,3,1,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
    make_sub_scatter(data[var2],data[f'end2_{name}'],varlab2,f'{name} dynamo end /Myr',3,3,2,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
    make_sub_scatter(data[var2],data[f'duration2_{name}'],varlab2,f'{name}dynamo duration /Myr',3,3,3,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
#CIA dynamo
name = 'CIA'
if np.all(np.isnan(data[f'start2_{name}'])==True):
    pass #don't plot this dynamo
else:
    make_sub_scatter(data[var2],data[f'start2_{name}'],varlab2,f'{name} dynamo start /Myr',3,3,4,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
    make_sub_scatter(data[var2],data[f'end2_{name}'],varlab2,f'{name} dynamo end /Myr',3,3,5,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
    make_sub_scatter(data[var2],data[f'duration2_{name}'],varlab2,f'{name}dynamo duration /Myr',3,3,6,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
#compositional dynamo
name = 'comp'
if np.all(np.isnan(data[f'start2_{name}'])==True):
    pass #don't plot this dynamo
else:
    make_sub_scatter(data[var2],data[f'start2_{name}'],varlab2,f'{name} dynamo start /Myr',3,3,7,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
    make_sub_scatter(data[var2],data[f'end2_{name}'],varlab2,f'{name} dynamo end /Myr',3,3,8,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
    make_sub_scatter(data[var2],data[f'duration2_{name}'],varlab2,f'{name}dynamo duration /Myr',3,3,9,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
#if save == True:
plt.savefig(f'../Plots/{folder}dynamo_duration2_{var2}{var1}_col.png')
#### var2 as x axis, y axis as var1, colour as dynamo value
# plt.figure(tight_layout=True,figsize=[15,15])
# plt.suptitle(f'Effect of {varlab2} and {varlab1} on dynamo timing')
# #MAC dynamo
# name = 'MAC'
# if np.all(np.isnan(data[f'start1_{name}'])==True):
#     pass #don't plot this dynamo
# else:
#     make_sub_scatter(data[var2],data[var1],varlab2,varlab1,3,3,1,colour=data[f'start1_{name}'],log=[logval2,logval1,False],colourlabel=f'{name} dynamo start /Myr',ss=5)
#     make_sub_scatter(data[var2],data[var1],varlab2,varlab1,3,3,2,colour=data[f'end1_{name}'],log=[logval2,logval1,False],colourlabel=f'{name} dynamo end /Myr',ss=5)
#     make_sub_scatter(data[var2],data[var1],varlab2,varlab1,3,3,3,colour=data[f'duration1_{name}'],log=[logval2,logval1,False],colourlabel=f'{name}dynamo duration /Myr',ss=5)
# #CIA dynamo
# name = 'CIA'
# make_sub_scatter(data[var2],data[var1],varlab2,varlab1,3,3,4,colour=data[f'start1_{name}'],log=[logval2,logval1,False],colourlabel=f'{name} dynamo start /Myr',ss=5)
# make_sub_scatter(data[var2],data[var1],varlab2,varlab1,3,3,5,colour=data[f'end1_{name}'],log=[logval2,logval1,False],colourlabel=f'{name} dynamo end /Myr',ss=5)
# make_sub_scatter(data[var2],data[var1],varlab2,varlab1,3,3,6,colour=data[f'duration1_{name}'],log=[logval2,logval1,False],colourlabel=f'{name}dynamo duration /Myr',ss=5)
# #compositional dynamo
# name = 'comp'
# make_sub_scatter(data[var2],data[var1],varlab2,varlab1,3,3,7,colour=data[f'start1_{name}'],log=[logval2,logval1,False],colourlabel=f'{name} dynamo start /Myr',ss=5)
# make_sub_scatter(data[var2],data[var1],varlab2,varlab1,3,3,8,colour=data[f'end1_{name}'],log=[logval2,logval1,False],colourlabel=f'{name} dynamo end /Myr',ss=5)
# make_sub_scatter(data[var2],data[var1],varlab2,varlab1,3,3,9,colour=data[f'duration1_{name}'],log=[logval2,logval1,False],colourlabel=f'{name}dynamo duration /Myr',ss=5)
# if save == True:
#     plt.savefig(f'../Plots/{folder}dynamo_duration1_{var2}{var1}.png')
    
#### var1 as x axis, radius as size
# plt.figure(tight_layout=True,figsize=[15,15])
# plt.suptitle(f'Effect of {varlab1} on dynamo timing')
# #MAC dynamo
# name = 'MAC'
# if np.all(np.isnan(data[f'start1_{name}'])==True):
#     pass #don't plot this dynamo
# else:
#     make_sub_scatter(data[var1],data[f'start1_{name}'],varlab1,f'{name} dynamo start /Myr',3,3,1,log=[logval1,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
#     make_sub_scatter(data[var1],data[f'end1_{name}'],varlab1,f'{name} dynamo end /Myr',3,3,2,log=[logval1,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
#     make_sub_scatter(data[var1],data[f'duration1_{name}'],varlab1,f'{name}dynamo duration/Myr',3,3,3,log=[logval1,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
# #CIA dynamo
# name = 'CIA'
# make_sub_scatter(data[var1],data[f'start1_{name}'],varlab1,f'{name} dynamo start /Myr',3,3,4,log=[logval1,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
# make_sub_scatter(data[var1],data[f'end1_{name}'],varlab1,f'{name} dynamo end /Myr',3,3,5,log=[logval1,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
# make_sub_scatter(data[var1],data[f'duration1_{name}'],varlab1,f'{name}dynamo duration /Myr',3,3,6,log=[logval1,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
# #compositional dynamo
# name = 'comp'
# make_sub_scatter(data[var1],data[f'start1_{name}'],varlab1,f'{name} dynamo start /Myr',3,3,7,log=[logval1,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
# make_sub_scatter(data[var1],data[f'end1_{name}'],varlab1,f'{name} dynamo end /Myr',3,3,8,log=[logval1,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
# make_sub_scatter(data[var1],data[f'duration1_{name}'],varlab1,f'{name}dynamo duration /Myr',3,3,9,log=[logval1,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
# if save == True:
#     plt.savefig(f'../Plots/{folder}dynamo_duration1_{var1}.png')
    
##################### Effect on dynamo strength ###############################
#### var2 as x axis, var1 as colour
plt.figure(tight_layout=True,figsize=[15,15])
plt.suptitle('Maximum field strength')
names = ['ml','mac','cia','comp']
label = ['flux balance','MAC','CIA','compositional']
for i, name,  in enumerate(names):
    make_sub_scatter(data[var2],data[f'maxB_{name}']/1e-6,varlab2,f'max B {label[i]} /$\mu$T',4,1,i+1,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
if save == True:
    plt.savefig(f'../Plots/{folder}dynamo_strength_{var2}{var1}.png',dpi=450)

##################### Peak temperature ###############################
#### second variable as x axis, first variable as colour
make_scatter(data[var2], data['peakT'], varlab2, 'peak mantle temperature /K',log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=2.5)
if save == True:
    plt.savefig(f'../Plots/{folder}peakT_axis_{var2}{var1}.png',dpi=450) 

#### var2 as x axis, var1 as y axis time of peak as colour
make_scatter(data[var2], data['peakT'], varlab2, 'peak mantle temperature /K',colour=data['tmax'],colourlabel='time of peak/Myr',ss=2.5)
if save == True:
   plt.savefig(f'../Plots/{folder}peakT_time.png',dpi=450)
   
   
#### peak core temperature
#### var2 as xaxis, var1 as y, peak T as colour
# make_scatter(data[var2],data[var1],varlab2,varlab1,colour=data['peak_coreT'],log=[logval2,logval1,False],colourlabel='peak core temperature /K',ss=2.5)
# if save == True:
#    plt.savefig(f'../Plots/{folder}peakcoreT_colour_{var1}.png',dpi=450)

#### var2 as x axis, var1 as colour
make_scatter(data[var2], data['peak_coreT'], varlab2, 'peak core temperature /K',log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=2.5)
if save == True:
    plt.savefig(f'../Plots/{folder}peakcoreT_axis_{var1}.png',dpi=450) 

#### var2 as x axis, time of peak as colour
make_scatter(data[var2], data['peak_coreT'], varlab2, 'peak core temperature /K',colour=data['tcoremax'],colourlabel='time of peak/Myr',ss=2.5)
if save == True:
   plt.savefig(f'../Plots/{folder}peakcoreT_time.png',dpi=450)   
   
#################### Differentiation time ####################################
#second variable as x axis, first variable as colour
make_scatter(data[var2], data['diff_time'], varlab2, 'differentiation time /Myr',log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=2.5)
if save == True:
   plt.savefig(f'../Plots/{folder}difft_axis_{var1}.png',dpi=450) 

#################### Core solidification time ###############################
#### second variable as x axis, first variable as colour
make_scatter(data[var2], data['tsolid'], varlab2, 'core solidification time /Myr',log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=2.5)
if save == True:
   plt.savefig(f'../Plots/{folder}tsolid_axis_{var1}.png',dpi=450) 

#################### Thermal history timings ###############################
#only one variable on y axis, keep all other variables fixed and filter data 
#for final version could try James suggestion and only display a few y axis values in horizontal bars
#this plot needs customising each time
#filter for fixed value of var2
var1 ='w'
data_fil = data#[data[var2]==var2val]
plt.figure(figsize=[15,5])
if var1 =='frht':
    plt.title(f'Thermal history and dynamo timings \n for r={r}, X$_{{s,0}}$={Xs_0}, rcmf={rcmf}, $\\eta_0$={eta0}, $^{{60}}$Fe/$^{{56}}$Fe ={Fe0} ')
elif var1 == 'rcmf':
    plt.title(f'Thermal history and dynamo timings \n for r={r}, X$_{{s,0}}$={Xs_0}, frht={frht}, $\\eta_0$={eta0}, $^{{60}}$Fe/$^{{56}}$Fe ={Fe0}')
elif var1 =='eta0':
    plt.title(f'Thermal history and dynamo timings \n for r={r}, X$_{{s,0}}$={Xs_0}, frht{frht}, rcmf={rcmf}')
elif var1 =='Fe0':
    plt.title(f'Thermal history and dynamo timings \n for r={r}, X$_{{s,0}}$={Xs_0}, frht{frht}, rcmf={rcmf},$\\eta_0$={eta0},')
plt.scatter(data_fil['diff_time'],data_fil[var1],label='differentiation',marker='s',color='darkred')
plt.scatter(data_fil['tmax'],data_fil[var1],label='peak mantle temp',marker='+',color='firebrick')
plt.scatter(data_fil['tcoremax'],data_fil[var1],label='peak core temp',marker='x',color='rosybrown')
plt.scatter(data_fil['tstrat_remove'],data_fil[var1],label='erosion of core stratification',marker='o',color='indianred')
plt.scatter(data_fil['terode'],data_fil[var1],label='core stratification removed',marker='p',color='lightcoral')
plt.scatter(data_fil['start1_MAC'],data_fil[var1],label='MAC dynamo start',marker='<',color='darkblue')
plt.scatter(data_fil['end1_MAC'],data_fil[var1],label='MAC dynamo end',marker='>',color='darkblue')
plt.scatter(data_fil['start1_CIA'],data_fil[var1],label='CIA dynamo start',marker='3',color='royalblue')
plt.scatter(data_fil['end1_CIA'],data_fil[var1],label='CIA dynamo end',marker='4',color='royalblue')
plt.scatter(data_fil['fcond_t'],data_fil[var1],label='end of mantle convection',marker='*',color='palevioletred')
plt.scatter(data_fil['start1_comp'],data_fil[var1],label='compositional dynamo start',marker='^',color='skyblue')
plt.scatter(data_fil['end1_comp'],data_fil[var1],label='compositional dynamo end',marker='v',color='skyblue')
plt.scatter(data_fil['tsolid'],data_fil[var1],label='core solidified',marker='d',color='pink')
plt.xlabel('Time/Myr')
plt.xscale('log')
if logval1 == True:
    plt.yscale('log')
plt.ylabel(varlab1)
plt.legend(ncols=2,bbox_to_anchor=(1,0.9))
if save == True:
    plt.savefig(f'../Plots/{folder}all_timings_{var1}.png',dpi=450,bbox_inches='tight') 

