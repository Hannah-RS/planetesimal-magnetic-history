#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make scatter plots for Xs-r test run
"""
from scatter_function import make_scatter, make_sub_scatter

import matplotlib.pyplot as plt
import sys
# setting path
sys.path.append('../')
from load_info import combine_info

data = combine_info('../Results_combined/Xs_r_tests/','auto_params.csv','run_results.csv',['MAC_onoff.csv','CIA_onoff.csv','comp_onoff.csv','coreconv_onoff.csv'])
data['r']=data['r']/1e3 #rescale to km

save = False #do you want to save the figures
###################### Effect on dynamo timing ################################
#r as x axis, sulfur as colour
plt.figure(tight_layout=True,figsize=[15,15])
plt.suptitle('Effect of radius on dynamo timing')
#MAC dynamo
name = 'MAC'
make_sub_scatter(data['r'],data[f'start_{name}'],'radius/km',f'{name} dynamo start /Myr',3,3,1,colour=data['Xs_0'],colourlabel='X$_{S,0}$',ss=5)
make_sub_scatter(data['r'],data[f'end_{name}'],'radius/km',f'{name} dynamo end /Myr',3,3,2,colour=data['Xs_0'],colourlabel='X$_{S,0}$',ss=5)
make_sub_scatter(data['r'],data[f'duration_{name}'],'radius/km',f'{name}dynamo duration /Myr',3,3,3,colour=data['Xs_0'],colourlabel='X$_{S,0}$',ss=5)
#CIA dynamo
name = 'CIA'
make_sub_scatter(data['r'],data[f'start_{name}'],'radius/km',f'{name} dynamo start /Myr',3,3,4,colour=data['Xs_0'],colourlabel='X$_{S,0}$',ss=5)
make_sub_scatter(data['r'],data[f'end_{name}'],'radius/km',f'{name} dynamo end /Myr',3,3,5,colour=data['Xs_0'],colourlabel='X$_{S,0}$',ss=5)
make_sub_scatter(data['r'],data[f'duration_{name}'],'radius/km',f'{name}dynamo duration /Myr',3,3,6,colour=data['Xs_0'],colourlabel='X$_{S,0}$',ss=5)
#compositional dynamo
name = 'comp'
make_sub_scatter(data['r'],data[f'start_{name}'],'radius/km',f'{name} dynamo start /Myr',3,3,7,colour=data['Xs_0'],colourlabel='X$_{S,0}$',ss=5)
make_sub_scatter(data['r'],data[f'end_{name}'],'radius/km',f'{name} dynamo end /Myr',3,3,8,colour=data['Xs_0'],colourlabel='X$_{S,0}$',ss=5)
make_sub_scatter(data['r'],data[f'duration_{name}'],'radius/km',f'{name}dynamo duration /Myr',3,3,9,colour=data['Xs_0'],colourlabel='X$_{S,0}$',ss=5)
if save == True:
    plt.savefig('../Plots/Xs_r_tests/dynamo_duration_r.png')

#r as x axis, y axis as sulfur, colour as dynamo value
plt.figure(tight_layout=True,figsize=[15,15])
plt.suptitle('Effect of radius and sulfur content on dynamo timing')
#MAC dynamo
name = 'MAC'
make_sub_scatter(data['r'],data['Xs_0'],'radius/km','X$_{S,0}$',3,3,1,colour=data[f'start_{name}'],colourlabel=f'{name} dynamo start /Myr',ss=5)
make_sub_scatter(data['r'],data['Xs_0'],'radius/km','X$_{S,0}$',3,3,2,colour=data[f'end_{name}'],colourlabel=f'{name} dynamo end /Myr',ss=5)
make_sub_scatter(data['r'],data['Xs_0'],'radius/km','X$_{S,0}$',3,3,3,colour=data[f'duration_{name}'],colourlabel=f'{name}dynamo duration /Myr',ss=5)
#CIA dynamo
name = 'CIA'
make_sub_scatter(data['r'],data['Xs_0'],'radius/km','X$_{S,0}$',3,3,4,colour=data[f'start_{name}'],colourlabel=f'{name} dynamo start /Myr',ss=5)
make_sub_scatter(data['r'],data['Xs_0'],'radius/km','X$_{S,0}$',3,3,5,colour=data[f'end_{name}'],colourlabel=f'{name} dynamo end /Myr',ss=5)
make_sub_scatter(data['r'],data['Xs_0'],'radius/km','X$_{S,0}$',3,3,6,colour=data[f'duration_{name}'],colourlabel=f'{name}dynamo duration /Myr',ss=5)
#compositional dynamo
name = 'comp'
make_sub_scatter(data['r'],data['Xs_0'],'radius/km','X$_{S,0}$',3,3,7,colour=data[f'start_{name}'],colourlabel=f'{name} dynamo start /Myr',ss=5)
make_sub_scatter(data['r'],data['Xs_0'],'radius/km','X$_{S,0}$',3,3,8,colour=data[f'end_{name}'],colourlabel=f'{name} dynamo end /Myr',ss=5)
make_sub_scatter(data['r'],data['Xs_0'],'radius/km','X$_{S,0}$',3,3,9,colour=data[f'duration_{name}'],colourlabel=f'{name}dynamo duration /Myr',ss=5)
if save == True:
    plt.savefig('../Plots/Xs_r_tests/dynamo_duration_rxs.png')
    
#sulfur as x axis, radius as size
plt.figure(tight_layout=True,figsize=[15,15])
plt.suptitle('Effect of initial sulfur content on dynamo timing')
#MAC dynamo
name = 'MAC'
make_sub_scatter(data['Xs_0'],data[f'start_{name}'],'X$_{s,0}$/ wt%',f'{name} dynamo start /Myr',3,3,1,size=data['r'],sizelabel='radius /km',ss=1.5)
make_sub_scatter(data['Xs_0'],data[f'end_{name}'],'X$_{s,0}$/ wt%',f'{name} dynamo end /Myr',3,3,2,size=data['r'],sizelabel='radius /km',ss=1.5)
make_sub_scatter(data['Xs_0'],data[f'duration_{name}'],'X$_{s,0}$/ wt%',f'{name}dynamo duration/Myr',3,3,3,size=data['r'],sizelabel='radius /km',ss=1.5)
#CIA dynamo
name = 'CIA'
make_sub_scatter(data['Xs_0'],data[f'start_{name}'],'X$_{s,0}$/ wt%',f'{name} dynamo start /Myr',3,3,4,size=data['r'],sizelabel='radius /km',ss=1.5)
make_sub_scatter(data['Xs_0'],data[f'end_{name}'],'X$_{s,0}$/ wt%',f'{name} dynamo end /Myr',3,3,5,size=data['r'],sizelabel='radius /km',ss=1.5)
make_sub_scatter(data['Xs_0'],data[f'duration_{name}'],'X$_{s,0}$/ wt%',f'{name}dynamo duration /Myr',3,3,6,size=data['r'],sizelabel='radius /km',ss=1.5)
#compositional dynamo
name = 'comp'
make_sub_scatter(data['Xs_0'],data[f'start_{name}'],'X$_{s,0}$/ wt%',f'{name} dynamo start /Myr',3,3,7,size=data['r'],sizelabel='radius /km',ss=1.5)
make_sub_scatter(data['Xs_0'],data[f'end_{name}'],'X$_{s,0}$/ wt%',f'{name} dynamo end /Myr',3,3,8,size=data['r'],sizelabel='radius /km',ss=1.5)
make_sub_scatter(data['Xs_0'],data[f'duration_{name}'],'X$_{s,0}$/ wt%',f'{name}dynamo duration /Myr',3,3,9,size=data['r'],sizelabel='radius /km',ss=1.5)
if save == True:
    plt.savefig('../Plots/Xs_r_tests/dynamo_duration_Xs.png')
    
##################### Effect on dynamo strength ###############################
#r as x axis, sulfur as colour
plt.figure(tight_layout=True,figsize=[15,15])
plt.suptitle('Maximum field strength')
names = ['ml','mac','cia','comp']
label = ['flux balance','MAC','CIA','compositional']
for i, name,  in enumerate(names):
    make_sub_scatter(data['r'],data[f'maxB_{name}']/1e-6,'radius/km',f'max B {label[i]} /$\mu$T',4,1,i+1,colour=data['Xs_0'],colourlabel='X$_{S,0}$',ss=5)
if save == True:
    plt.savefig('../Plots/Xs_r_tests/dynamo_strength.png',dpi=450)

##################### Peak temperature ###############################
#r as x axis, sulfur as colour
make_scatter(data['r'], data['peakT'], 'radius /km', 'peak mantle temperature /K',colour=data['Xs_0'],colourlabel='X$_{S,0}$',ss=2.5)
if save == True:
   plt.savefig('../Plots/Xs_r_tests/peakT_axis.png',dpi=450) 

#r as x axis, time of peak as colour
make_scatter(data['r'], data['peakT'], 'radius /km', 'peak mantle temperature /K',colour=data['tmax'],colourlabel='time of peak/Myr',ss=2.5)
if save == True:
   plt.savefig('../Plots/Xs_r_tests/peakT_time.png',dpi=450)
   
# r as xaxis, sulfur as y, peak T as colour
make_scatter(data['r'],data['Xs_0'],'radius /km','X$_{S,0}$',colour=data['peakT'],colourlabel='peak mantle temperature /K',ss=2.5)
if save == True:
   plt.savefig('../Plots/Xs_r_tests/peakT_colour.png',dpi=450) 
   
# peak core temperature
# r as xaxis, sulfur as y, peak T as colour
make_scatter(data['r'],data['Xs_0'],'radius /km','X$_{S,0}$',colour=data['peak_coreT'],colourlabel='peak core temperature /K',ss=2.5)
if save == True:
   plt.savefig('../Plots/Xs_r_tests/peakcoreT_colour.png',dpi=450)

#r as x axis, sulfur as colour
make_scatter(data['r'], data['peak_coreT'], 'radius /km', 'peak core temperature /K',colour=data['Xs_0'],colourlabel='X$_{S,0}$',ss=2.5)
if save == True:
   plt.savefig('../Plots/Xs_r_tests/peakcoreT_axis.png',dpi=450) 

#r as x axis, time of peak as colour
make_scatter(data['r'], data['peak_coreT'], 'radius /km', 'peak core temperature /K',colour=data['tcoremax'],colourlabel='time of peak/Myr',ss=2.5)
if save == True:
   plt.savefig('../Plots/Xs_r_tests/peakcoreT_time.png',dpi=450)   
#################### Thermal history timings ###############################
#y axis is r, x axis is time
plt.figure(figsize=[15,5])
plt.title('Thermal history timing')
plt.scatter(data['diff_time'],data['r'],label='differentiation',marker='x',c=data['Xs_0'])
plt.scatter(data['tmax'],data['r'],label='peak mantle temp',marker='+',c=data['Xs_0'])
plt.scatter(data['tstrat_remove'],data['r'],label='erosion of core stratification',marker='o',c=data['Xs_0'])
plt.scatter(data['terode'],data['r'],label='core stratification removed',marker='v',c=data['Xs_0'])
plt.scatter(data['fcond_t'],data['r'],label='end of mantle convection',marker='*',c=data['Xs_0'])
plt.scatter(data['tsolid'],data['r'],label='core solidified',marker='^',c=data['Xs_0'])
plt.xlabel('Time/Myr')
plt.xscale('log')
plt.ylabel('radius/km')
plt.colorbar(label='x$_{S,0}$')
plt.legend(ncols=1,bbox_to_anchor=(1.45,0.9))
if save == True:
    plt.savefig('../Plots/Xs_r_tests/thermal_timings.png',dpi=450,bbox_inches='tight')

#dynamo timings only
plt.figure(figsize=[15,5])
plt.title('Dynamo timings')
plt.scatter(data['start_MAC'],data['r'],label='MAC dynamo start',marker='<',c=data['Xs_0'])
plt.scatter(data['end_MAC'],data['r'],label='MAC dynamo end',marker='>',c=data['Xs_0'])
plt.scatter(data['start_CIA'],data['r'],label='CIA dynamo start',marker='3',c=data['Xs_0'])
plt.scatter(data['end_CIA'],data['r'],label='CIA dynamo end',marker='4',c=data['Xs_0'])
plt.scatter(data['start_comp'],data['r'],label='compositional dynamo start',marker='^',c=data['Xs_0'])
plt.scatter(data['end_comp'],data['r'],label='compositional dynamo end',marker='v',c=data['Xs_0'])
plt.xlabel('Time/Myr')
plt.xscale('log')
plt.ylabel('radius/km')
plt.colorbar(label='x$_{S,0}$')
plt.legend(ncols=1,bbox_to_anchor=(1.45,0.9))
ax = plt.gca()
leg = ax.get_legend()
for i, lab in enumerate(leg.legendHandles):
    lab.set_edgecolor('black')
    lab.set_facecolor([[0,0,0,1]])
    #lab.set_facecolor('black')
if save == True:
   plt.savefig('../Plots/Xs_r_tests/dynamo_timings.png',dpi=450,bbox_inches='tight') 

#dynamo and thermal timings
#y axis is r, x axis is time
plt.figure(figsize=[15,5])
plt.title('Thermal history and dynamo timings')
plt.scatter(data['diff_time'],data['r'],label='differentiation',marker='x',c=data['Xs_0'])
plt.scatter(data['tmax'],data['r'],label='peak mantle temp',marker='+',c=data['Xs_0'])
plt.scatter(data['tstrat_remove'],data['r'],label='erosion of core stratification',marker='o',c=data['Xs_0'])
plt.scatter(data['terode'],data['r'],label='core stratification removed',marker='p',c=data['Xs_0'])
plt.scatter(data['start_MAC'],data['r'],label='MAC dynamo start',marker='<',c=data['Xs_0'])
plt.scatter(data['end_MAC'],data['r'],label='MAC dynamo end',marker='>',c=data['Xs_0'])
plt.scatter(data['start_CIA'],data['r'],label='CIA dynamo start',marker='3',c=data['Xs_0'])
plt.scatter(data['end_CIA'],data['r'],label='CIA dynamo end',marker='4',c=data['Xs_0'])
plt.scatter(data['fcond_t'],data['r'],label='end of mantle convection',marker='*',c=data['Xs_0'])
plt.scatter(data['start_comp'],data['r'],label='compositional dynamo start',marker='^',c=data['Xs_0'])
plt.scatter(data['end_comp'],data['r'],label='compositional dynamo end',marker='v',c=data['Xs_0'])
plt.scatter(data['tsolid'],data['r'],label='core solidified',marker='d',c=data['Xs_0'])
plt.xlabel('Time/Myr')
plt.xscale('log')
plt.ylabel('radius/km')
plt.colorbar(label='x$_{S,0}$')
plt.legend(ncols=1,bbox_to_anchor=(1.45,0.9))
ax = plt.gca()
leg = ax.get_legend()
for i, lab in enumerate(leg.legendHandles):
    lab.set_edgecolor('black')
    lab.set_facecolor([[0,0,0,1]])
    #lab.set_facecolor('black')
if save == True:
   plt.savefig('../Plots/Xs_r_tests/all_timings.png',dpi=450,bbox_inches='tight') 

#one sulfur content
#choose fixed sulfur values and filter data
Xs_val = 30 
data_fil = data[data['Xs_0']==Xs_val]
plt.figure(figsize=[15,5])
plt.title(f'Thermal history and dynamo timings \n core sulfur content = {Xs_val} wt%')
plt.scatter(data_fil['diff_time'],data_fil['r'],label='differentiation',marker='s',color='darkred')
plt.scatter(data_fil['tmax'],data_fil['r'],label='peak mantle temp',marker='+',color='firebrick')
plt.scatter(data_fil['tcoremax'],data_fil['r'],label='peak core temp',marker='x',color='rosybrown')
plt.scatter(data_fil['tstrat_remove'],data_fil['r'],label='erosion of core stratification',marker='o',color='indianred')
plt.scatter(data_fil['terode'],data_fil['r'],label='core stratification removed',marker='p',color='lightcoral')
plt.scatter(data_fil['start_MAC'],data_fil['r'],label='MAC dynamo start',marker='<',color='darkblue')
plt.scatter(data_fil['end_MAC'],data_fil['r'],label='MAC dynamo end',marker='>',color='darkblue')
plt.scatter(data_fil['start_CIA'],data_fil['r'],label='CIA dynamo start',marker='3',color='royalblue')
plt.scatter(data_fil['end_CIA'],data_fil['r'],label='CIA dynamo end',marker='4',color='royalblue')
plt.scatter(data_fil['fcond_t'],data_fil['r'],label='end of mantle convection',marker='*',color='palevioletred')
plt.scatter(data_fil['start_comp'],data_fil['r'],label='compositional dynamo start',marker='^',color='skyblue')
plt.scatter(data_fil['end_comp'],data_fil['r'],label='compositional dynamo end',marker='v',color='skyblue')
plt.scatter(data_fil['tsolid'],data_fil['r'],label='core solidified',marker='d',color='pink')
plt.xlabel('Time/Myr')
plt.xscale('log')
plt.ylabel('radius/km')
plt.legend(ncols=1,bbox_to_anchor=(1,0.9))
if save == True:
   plt.savefig('../Plots/Xs_r_tests/all_timings_fixS.png',dpi=450,bbox_inches='tight') 
   
#early times
plt.figure(figsize=[15,5])
plt.scatter(data['diff_time'],data['r'],label='differentiation',marker='x',c=data['Xs_0'])
plt.scatter(data['tmax'],data['r'],label='peak mantle temp',marker='+',c=data['Xs_0'])
plt.scatter(data['tstrat_remove'],data['r'],label='erosion of core stratification',marker='o',c=data['Xs_0'])
plt.scatter(data['terode'],data['r'],label='core stratification removed',marker='p',c=data['Xs_0'])
plt.scatter(data['start_MAC'],data['r'],label='MAC dynamo start',marker='<',c=data['Xs_0'])
plt.scatter(data['end_MAC'],data['r'],label='MAC dynamo end',marker='>',c=data['Xs_0'])
plt.scatter(data['start_CIA'],data['r'],label='CIA dynamo start',marker='3',c=data['Xs_0'])
plt.scatter(data['end_CIA'],data['r'],label='CIA dynamo end',marker='4',c=data['Xs_0'])
plt.scatter(data['fcond_t'],data['r'],label='end of mantle convection',marker='*',c=data['Xs_0'])
plt.xlabel('Time/Myr')
plt.xscale('log')
plt.xlim(right=10)
plt.ylabel('radius/km')
plt.colorbar(label='x$_{S,0}$')
plt.title('Thermal history and dynamo timings  \n 0.8-10Myr')
plt.legend(ncols=1,bbox_to_anchor=(1.45,0.8))
ax = plt.gca()
leg = ax.get_legend()
for i, lab in enumerate(leg.legendHandles):
    lab.set_edgecolor('black')
    lab.set_facecolor([[0,0,0,1]])
    #lab.set_facecolor('black')
if save == True:
   plt.savefig('../Plots/Xs_r_tests/all_timings_early.png',dpi=450,bbox_inches='tight') 
   
#late times
plt.figure(figsize=[15,3.5])
plt.scatter(data['end_CIA'],data['r'],label='CIA dynamo end',marker='4',c=data['Xs_0'])
plt.scatter(data['start_comp'],data['r'],label='compositional dynamo start',marker='^',c=data['Xs_0'])
plt.scatter(data['end_comp'],data['r'],label='compositional dynamo end',marker='v',c=data['Xs_0'])
plt.scatter(data['tsolid'],data['r'],label='core solidified',marker='d',c=data['Xs_0'])
plt.xlabel('Time/Myr')
plt.xscale('log')
plt.xlim(left=10)
plt.ylabel('radius/km')
plt.colorbar(label='x$_{S,0}$')
plt.title('Thermal history and dynamo timings  \n >10Myr')
plt.legend(ncols=1)
ax = plt.gca()
leg = ax.get_legend()
for i, lab in enumerate(leg.legendHandles):
    lab.set_edgecolor('black')
    lab.set_facecolor([[0,0,0,1]])
    #lab.set_facecolor('black')
if save == True:
   plt.savefig('../Plots/Xs_r_tests/all_timings_late.png',dpi=450,bbox_inches='tight')
