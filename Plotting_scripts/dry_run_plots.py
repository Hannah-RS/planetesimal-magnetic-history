#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make scatter plots for Xs-r test run
"""
from analysis_functions import combine_info
from scatter_function import make_scatter, make_sub_scatter

import matplotlib.pyplot as plt
import sys
# setting path
sys.path.append('../')
from plotting_constants import Myr

data = combine_info('../Results_combined/Xs_r_tests/','auto_params.csv','run_results.csv',['MAC_onoff.csv','CIA_onoff.csv','comp_onoff.csv','coreconv_onoff.csv'])
data['r']=data['r']/1e3 #rescale to km
data['cond_t'] = data['cond_t']/Myr #rescale by Myr

save = False #do you want to save the figures
###################### Effect on dynamo timing ################################
#r as x axis, sulfur as colour
plt.figure(tight_layout=True,figsize=[15,15])
#MAC dynamo
name = 'MAC'
make_sub_scatter(data['r'],data[f'start_{name}'],'radius/km',f'{name} dynamo start /Myr',3,3,1,colour=data['Xs_0'],colourlabel='X$_{S,0}$',ss=5)
make_sub_scatter(data['r'],data[f'end_{name}'],'radius/km',f'{name} dynamo end /Myr',3,3,2,colour=data['Xs_0'],colourlabel='X$_{S,0}$',ss=5)
make_sub_scatter(data['r'],data[f'duration_{name}'],'radius/km',f'{name}dynamo start /Myr',3,3,3,colour=data['Xs_0'],colourlabel='X$_{S,0}$',ss=5)
#CIA dynamo
name = 'CIA'
make_sub_scatter(data['r'],data[f'start_{name}'],'radius/km',f'{name} dynamo start /Myr',3,3,4,colour=data['Xs_0'],colourlabel='X$_{S,0}$',ss=5)
make_sub_scatter(data['r'],data[f'end_{name}'],'radius/km',f'{name} dynamo end /Myr',3,3,5,colour=data['Xs_0'],colourlabel='X$_{S,0}$',ss=5)
make_sub_scatter(data['r'],data[f'duration_{name}'],'radius/km',f'{name}dynamo start /Myr',3,3,6,colour=data['Xs_0'],colourlabel='X$_{S,0}$',ss=5)
#compositional dynamo
name = 'comp'
make_sub_scatter(data['r'],data[f'start_{name}'],'radius/km',f'{name} dynamo start /Myr',3,3,7,colour=data['Xs_0'],colourlabel='X$_{S,0}$',ss=5)
make_sub_scatter(data['r'],data[f'end_{name}'],'radius/km',f'{name} dynamo end /Myr',3,3,8,colour=data['Xs_0'],colourlabel='X$_{S,0}$',ss=5)
make_sub_scatter(data['r'],data[f'duration_{name}'],'radius/km',f'{name}dynamo start /Myr',3,3,9,colour=data['Xs_0'],colourlabel='X$_{S,0}$',ss=5)
if save == True:
    plt.savefig('../Plots/Xs_r_tests/dynamo_duration_r.png')

#sulfur as x axis, radius as size
plt.figure(tight_layout=True,figsize=[15,15])
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
names = ['ml','mac','cia','comp']
for i, name in enumerate(names):
    make_sub_scatter(data['r'],data[f'maxB_{name}'],'radius/km',f'max B {name}/$\mu$T',4,1,i+1,colour=data['Xs_0'],colourlabel='X$_{S,0}$',ss=5)
if save == True:
    plt.savefig('../Plots/Xs_r_tests/dynamo_strength.png')
