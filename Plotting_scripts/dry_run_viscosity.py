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

folder = 'Fullrun1/'
data_in = pd.read_csv(f'../Results_combined/{folder}/all_sucess_info.csv',delimiter=',',skiprows=[1],header=0,index_col=False)
data_in['r']=data_in['r']/1e3 #rescale to km

#filter for the variable you are interested in
var = 'eta0' #variable of interest for these plots
varlab = '$\\eta_0$' 
logval = True #should the yaxis of this variable be scaled logarithmically
#fixed values of other quantities
#eta0=1e21
rcmf = 0.5
frht=0.02375
Xs_0 = 30.25
Fe0 = 1e-7
if var =='frht':
    data = data_in[(data_in['eta0']==eta0)&(data_in['rcmf']==rcmf)&(data_in['Xs_0']==Xs_0)&(data_in['Fe0']==Fe0)]#runs where other two variables are reference values
elif var == 'rcmf':
    data = data_in[(data_in['eta0']==eta0)&(data_in['frht']==frht)&(data_in['Xs_0']==Xs_0)&(data_in['Fe0']==Fe0)]#runs where other two variables are reference values
elif var == 'eta0':
    data = data_in[(data_in['rcmf']==rcmf)&(data_in['frht']==frht)&(data_in['Xs_0']==Xs_0)&(data_in['Fe0']==Fe0)]
save = False #do you want to save the figures
###################### Effect on dynamo timing ################################
#r as x axis, rcmf as colour
x='r'
xlab = 'radius/km'
xlog = False
col = var
collab = varlab
plt.figure(tight_layout=True,figsize=[15,15])
plt.suptitle('Effect of radius on dynamo timing')
#MAC dynamo
name = 'MAC'
if np.all(np.isnan(data[f'start_{name}'])==True):
    pass #don't plot this dynamo
else:
    make_sub_scatter(data[x],data[f'start_{name}'],xlab,f'{name} dynamo start /Myr',3,3,1,log=[xlog,False,logval],colour=data[col],colourlabel=collab,ss=5)
    make_sub_scatter(data[x],data[f'end_{name}'],xlab,f'{name} dynamo end /Myr',3,3,2,log=[xlog,False,logval],colour=data[col],colourlabel=collab,ss=5)
    make_sub_scatter(data[x],data[f'duration_{name}'],xlab,f'{name}dynamo duration /Myr',3,3,3,log=[xlog,False,logval],colour=data[col],colourlabel=collab,ss=5)
#CIA dynamo
name = 'CIA'
make_sub_scatter(data[x],data[f'start_{name}'],xlab,f'{name} dynamo start /Myr',3,3,4,log=[xlog,False,logval],colour=data[col],colourlabel=collab,ss=5)
make_sub_scatter(data[x],data[f'end_{name}'],xlab,f'{name} dynamo end /Myr',3,3,5,log=[xlog,False,logval],colour=data[col],colourlabel=collab,ss=5)
make_sub_scatter(data[x],data[f'duration_{name}'],xlab,f'{name}dynamo duration /Myr',3,3,6,log=[xlog,False,logval],colour=data[col],colourlabel=collab,ss=5)
#compositional dynamo
name = 'comp'
make_sub_scatter(data[x],data[f'start_{name}'],xlab,f'{name} dynamo start /Myr',3,3,7,log=[xlog,False,logval],colour=data[col],colourlabel=collab,ss=5)
make_sub_scatter(data[x],data[f'end_{name}'],xlab,f'{name} dynamo end /Myr',3,3,8,log=[xlog,False,logval],colour=data[col],colourlabel=collab,ss=5)
make_sub_scatter(data[x],data[f'duration_{name}'],xlab,f'{name}dynamo duration /Myr',3,3,9,log=[xlog,False,logval],colour=data[col],colourlabel=collab,ss=5)
if save == True:
    plt.savefig(f'../Plots/{folder}dynamo_duration_r{var}_col.png')

#r as x axis, y axis as rcmf, colour as dynamo value
x='r'
xlab = 'radius/km'
xlog = False #shpuld x be plotted logarithmically
y = var
ylab = varlab
plt.figure(tight_layout=True,figsize=[15,15])
plt.suptitle('Effect of radius and critical melt fraction on dynamo timing')
#MAC dynamo
name = 'MAC'
if np.all(np.isnan(data[f'start_{name}'])==True):
    pass #don't plot this dynamo
else:
    make_sub_scatter(data[x],data[y],xlab,ylab,3,3,1,colour=data[f'start_{name}'],log=[xlog,logval,False],colourlabel=f'{name} dynamo start /Myr',ss=5)
    make_sub_scatter(data[x],data[y],xlab,ylab,3,3,2,colour=data[f'end_{name}'],log=[xlog,logval,False],colourlabel=f'{name} dynamo end /Myr',ss=5)
    make_sub_scatter(data[x],data[y],xlab,ylab,3,3,3,colour=data[f'duration_{name}'],log=[xlog,logval,False],colourlabel=f'{name}dynamo duration /Myr',ss=5)
#CIA dynamo
name = 'CIA'
make_sub_scatter(data[x],data[y],xlab,ylab,3,3,4,colour=data[f'start_{name}'],log=[xlog,logval,False],colourlabel=f'{name} dynamo start /Myr',ss=5)
make_sub_scatter(data[x],data[y],xlab,ylab,3,3,5,colour=data[f'end_{name}'],log=[xlog,logval,False],colourlabel=f'{name} dynamo end /Myr',ss=5)
make_sub_scatter(data[x],data[y],xlab,ylab,3,3,6,colour=data[f'duration_{name}'],log=[xlog,logval,False],colourlabel=f'{name}dynamo duration /Myr',ss=5)
#compositional dynamo
name = 'comp'
make_sub_scatter(data[x],data[y],xlab,ylab,3,3,7,colour=data[f'start_{name}'],log=[xlog,logval,False],colourlabel=f'{name} dynamo start /Myr',ss=5)
make_sub_scatter(data[x],data[y],xlab,ylab,3,3,8,colour=data[f'end_{name}'],log=[xlog,logval,False],colourlabel=f'{name} dynamo end /Myr',ss=5)
make_sub_scatter(data[x],data[y],xlab,ylab,3,3,9,colour=data[f'duration_{name}'],log=[xlog,logval,False],colourlabel=f'{name}dynamo duration /Myr',ss=5)
if save == True:
    plt.savefig(f'../Plots/{folder}dynamo_duration_r{var}.png')
    
#rcmf as x axis, radius as size
x = var
xlab = varlab
plt.figure(tight_layout=True,figsize=[15,15])
plt.suptitle(f'Effect of {varlab} on dynamo timing')
#MAC dynamo
name = 'MAC'
if np.all(np.isnan(data[f'start_{name}'])==True):
    pass #don't plot this dynamo
else:
    make_sub_scatter(data[x],data[f'start_{name}'],xlab,f'{name} dynamo start /Myr',3,3,1,log=[logval,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
    make_sub_scatter(data[x],data[f'end_{name}'],xlab,f'{name} dynamo end /Myr',3,3,2,log=[logval,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
    make_sub_scatter(data[x],data[f'duration_{name}'],xlab,f'{name}dynamo duration/Myr',3,3,3,log=[logval,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
#CIA dynamo
name = 'CIA'
make_sub_scatter(data[x],data[f'start_{name}'],xlab,f'{name} dynamo start /Myr',3,3,4,log=[logval,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
make_sub_scatter(data[x],data[f'end_{name}'],xlab,f'{name} dynamo end /Myr',3,3,5,log=[logval,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
make_sub_scatter(data[x],data[f'duration_{name}'],xlab,f'{name}dynamo duration /Myr',3,3,6,log=[logval,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
#compositional dynamo
name = 'comp'
make_sub_scatter(data[x],data[f'start_{name}'],xlab,f'{name} dynamo start /Myr',3,3,7,log=[logval,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
make_sub_scatter(data[x],data[f'end_{name}'],xlab,f'{name} dynamo end /Myr',3,3,8,log=[logval,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
make_sub_scatter(data[x],data[f'duration_{name}'],xlab,f'{name}dynamo duration /Myr',3,3,9,log=[logval,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
if save == True:
    plt.savefig(f'../Plots/{folder}dynamo_duration_{var}.png')
    
##################### Effect on dynamo strength ###############################
#r as x axis, sulfur as colour
plt.figure(tight_layout=True,figsize=[15,15])
plt.suptitle('Maximum field strength')
names = ['ml','mac','cia','comp']
label = ['flux balance','MAC','CIA','compositional']
for i, name,  in enumerate(names):
    make_sub_scatter(data['r'],data[f'maxB_{name}']/1e-6,'radius/km',f'max B {label[i]} /$\mu$T',4,1,i+1,log=[False,False,logval],colour=data[var],colourlabel=varlab,ss=5)
if save == True:
    plt.savefig(f'../Plots/{folder}dynamo_strength_{var}.png',dpi=450)

##################### Peak temperature ###############################
#r as x axis, rcmf as colour
make_scatter(data['r'], data['peakT'], 'radius /km', 'peak mantle temperature /K',log=[False,False,logval],colour=data[var],colourlabel=varlab,ss=2.5)
if save == True:
   plt.savefig(f'../Plots/{folder}peakT_axis_{var}.png',dpi=450) 

#r as x axis, time of peak as colour
make_scatter(data['r'], data['peakT'], 'radius /km', 'peak mantle temperature /K',colour=data['tmax'],colourlabel='time of peak/Myr',ss=2.5)
if save == True:
   plt.savefig(f'../Plots/{folder}peakT_time.png',dpi=450)
   
   
# peak core temperature
# r as xaxis, rcmf as y, peak T as colour
make_scatter(data['r'],data[var],'radius /km',varlab,colour=data['peak_coreT'],log=[False,logval,False],colourlabel='peak core temperature /K',ss=2.5)
if save == True:
   plt.savefig(f'../Plots/{folder}peakcoreT_colour_{var}.png',dpi=450)

#r as x axis, rcmf as colour
make_scatter(data['r'], data['peak_coreT'], 'radius /km', 'peak core temperature /K',log=[False,False,logval],colour=data[var],colourlabel=varlab,ss=2.5)
if save == True:
   plt.savefig(f'../Plots/{folder}peakcoreT_axis_{var}.png',dpi=450) 

#r as x axis, time of peak as colour
make_scatter(data['r'], data['peak_coreT'], 'radius /km', 'peak core temperature /K',colour=data['tcoremax'],colourlabel='time of peak/Myr',ss=2.5)
if save == True:
   plt.savefig(f'../Plots/{folder}peakcoreT_time.png',dpi=450)   
#################### Thermal history timings ###############################
#only one variable on y axis, keep all other variables fixed and filter data 
#for final version could try James suggestion and only display a few y axis values in horizontal bars
r=200
data_fil = data[(data['Xs_0']==Xs_0)&(data['r']==r)]
plt.figure(figsize=[15,5])
if var =='frht':
    plt.title(f'Thermal history and dynamo timings \n for r={r}, X$_{{s,0}}$={Xs_0}, rcmf={rcmf}, $\\eta_0$={eta0}, $^{{60}}$Fe/$^{{56}}$Fe ={Fe0} ')
elif var == 'rcmf':
    plt.title(f'Thermal history and dynamo timings \n for r={r}, X$_{{s,0}}$={Xs_0}, frht={frht}, $\\eta_0$={eta0}, $^{{60}}$Fe/$^{{56}}$Fe ={Fe0}')
elif var =='eta0':
    plt.title(f'Thermal history and dynamo timings \n for r={r}, X$_{{s,0}}$={Xs_0}, frht{frht}, rcmf={rcmf}')
plt.scatter(data_fil['diff_time'],data_fil[var],label='differentiation',marker='s',color='darkred')
plt.scatter(data_fil['tmax'],data_fil[var],label='peak mantle temp',marker='+',color='firebrick')
plt.scatter(data_fil['tcoremax'],data_fil[var],label='peak core temp',marker='x',color='rosybrown')
plt.scatter(data_fil['tstrat_remove'],data_fil[var],label='erosion of core stratification',marker='o',color='indianred')
plt.scatter(data_fil['terode'],data_fil[var],label='core stratification removed',marker='p',color='lightcoral')
plt.scatter(data_fil['start_MAC'],data_fil[var],label='MAC dynamo start',marker='<',color='darkblue')
plt.scatter(data_fil['end_MAC'],data_fil[var],label='MAC dynamo end',marker='>',color='darkblue')
plt.scatter(data_fil['start_CIA'],data_fil[var],label='CIA dynamo start',marker='3',color='royalblue')
plt.scatter(data_fil['end_CIA'],data_fil[var],label='CIA dynamo end',marker='4',color='royalblue')
plt.scatter(data_fil['fcond_t'],data_fil[var],label='end of mantle convection',marker='*',color='palevioletred')
plt.scatter(data_fil['start_comp'],data_fil[var],label='compositional dynamo start',marker='^',color='skyblue')
plt.scatter(data_fil['end_comp'],data_fil[var],label='compositional dynamo end',marker='v',color='skyblue')
plt.scatter(data_fil['tsolid'],data_fil[var],label='core solidified',marker='d',color='pink')
plt.xlabel('Time/Myr')
plt.xscale('log')
if logval == True:
    plt.yscale('log')
plt.ylabel(varlab)
plt.legend(ncols=1,bbox_to_anchor=(1,0.9))
if save == True:
   plt.savefig(f'../Plots/{folder}all_timings_{var}.png',dpi=450,bbox_inches='tight') 

