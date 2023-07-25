#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make scatter plots to compare effects of two variables on the following parameters:
Dynamo start
Dynamo stop
Duration
Ng l10
Ngl100
Peak B
End of convection
fconv_T
Diff t
Core solid
Peak mantle T
Peak core T
"""
from scatter_function import make_scatter, make_sub_scatter
import sys
# setting path
sys.path.append('../')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

folder = 'Fullrun2'
data = pd.read_csv(f'../Results_combined/{folder}/all_sucess_info.csv',delimiter=',',skiprows=[1],header=0,index_col=False)
data['r']=data['r']/1e3 #rescale to km

#filter for the variable you are interested in
var1 = 'r' #variable of interest for these plots
varlab1 = 'radius /km' 
logval1 = False #should the yaxis of this variable be scaled logarithmically
var2 = 'rcmf'
varlab2 = 'rcmf'
logval2 = False 
save = False #do you want to save the figures
#fixed values of other quantities
eta0=1e14
r = 100
rcmf = 0.2
frht=0.005
Xs_0 = 28.5
Fe0 = 0
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
    
if len(data)==0:
    raise ValueError('Invalid parameter choice - no data remains')
###################### Effect on dynamo timing ################################
#### var2 as x axis, var1 as colour
plt.figure(tight_layout=True,figsize=[15,15])
plt.suptitle(f'Effect of {varlab2} on dynamo timing')
#MAC dynamo
name = 'mac'
if np.all(np.isnan(data[f'{name}_on'])==True):
    pass #don't plot this dynamo
else:
    make_sub_scatter(data[var2],data[f'{name}_on'],varlab2,f'{name} dynamo start /Myr',3,3,1,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
    make_sub_scatter(data[var2],data[f'{name}_off'],varlab2,f'{name} dynamo end /Myr',3,3,2,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
    make_sub_scatter(data[var2],data[f'{name}_dur'],varlab2,f'{name}dynamo duration /Myr',3,3,3,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
#CIA dynamo
name = 'cia'
make_sub_scatter(data[var2],data[f'{name}_on'],varlab2,f'{name} dynamo start /Myr',3,3,4,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
make_sub_scatter(data[var2],data[f'{name}_off'],varlab2,f'{name} dynamo end /Myr',3,3,5,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
make_sub_scatter(data[var2],data[f'{name}_dur'],varlab2,f'{name}dynamo duration /Myr',3,3,6,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
#compositional dynamo
name = 'comp'
make_sub_scatter(data[var2],data[f'{name}_on'],varlab2,f'{name} dynamo start /Myr',3,3,7,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
make_sub_scatter(data[var2],data[f'{name}_off'],varlab2,f'{name} dynamo end /Myr',3,3,8,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
make_sub_scatter(data[var2],data[f'{name}_dur'],varlab2,f'{name}dynamo duration /Myr',3,3,9,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
if save == True:
    plt.savefig(f'../Plots/{folder}dynamo_duration1_{var2}{var1}_col.png')

#### var2 as x axis, y axis as var1, colour as dynamo value
# plt.figure(tight_layout=True,figsize=[15,15])
# plt.suptitle(f'Effect of {varlab2} and {varlab1} on dynamo timing')
# #MAC dynamo
# name = 'mac'
# if np.all(np.isnan(data[f'{name}_on'])==True):
#     pass #don't plot this dynamo
# else:
#     make_sub_scatter(data[var2],data[var1],varlab2,varlab1,3,3,1,colour=data[f'{name}_on'],log=[logval2,logval1,False],colourlabel=f'{name} dynamo start /Myr',ss=5)
#     make_sub_scatter(data[var2],data[var1],varlab2,varlab1,3,3,2,colour=data[f'{name}_off'],log=[logval2,logval1,False],colourlabel=f'{name} dynamo end /Myr',ss=5)
#     make_sub_scatter(data[var2],data[var1],varlab2,varlab1,3,3,3,colour=data[f'{name}_dur'],log=[logval2,logval1,False],colourlabel=f'{name}dynamo duration /Myr',ss=5)
# #CIA dynamo
# name = 'cia'
# make_sub_scatter(data[var2],data[var1],varlab2,varlab1,3,3,4,colour=data[f'{name}_on'],log=[logval2,logval1,False],colourlabel=f'{name} dynamo start /Myr',ss=5)
# make_sub_scatter(data[var2],data[var1],varlab2,varlab1,3,3,5,colour=data[f'{name}_off'],log=[logval2,logval1,False],colourlabel=f'{name} dynamo end /Myr',ss=5)
# make_sub_scatter(data[var2],data[var1],varlab2,varlab1,3,3,6,colour=data[f'{name}_dur'],log=[logval2,logval1,False],colourlabel=f'{name}dynamo duration /Myr',ss=5)
# #compositional dynamo
# name = 'comp'
# make_sub_scatter(data[var2],data[var1],varlab2,varlab1,3,3,7,colour=data[f'{name}_on'],log=[logval2,logval1,False],colourlabel=f'{name} dynamo start /Myr',ss=5)
# make_sub_scatter(data[var2],data[var1],varlab2,varlab1,3,3,8,colour=data[f'{name}_off'],log=[logval2,logval1,False],colourlabel=f'{name} dynamo end /Myr',ss=5)
# make_sub_scatter(data[var2],data[var1],varlab2,varlab1,3,3,9,colour=data[f'{name}_dur'],log=[logval2,logval1,False],colourlabel=f'{name}dynamo duration /Myr',ss=5)
# if save == True:
#     plt.savefig(f'../Plots/{folder}dynamo_duration1_{var2}{var1}.png')
    
#### var1 as x axis, radius as size
# plt.figure(tight_layout=True,figsize=[15,15])
# plt.suptitle(f'Effect of {varlab1} on dynamo timing')
# #MAC dynamo
# name = 'mac'
# if np.all(np.isnan(data[f'{name}_on'])==True):
#     pass #don't plot this dynamo
# else:
#     make_sub_scatter(data[var1],data[f'{name}_on'],varlab1,f'{name} dynamo start /Myr',3,3,1,log=[logval1,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
#     make_sub_scatter(data[var1],data[f'{name}_off'],varlab1,f'{name} dynamo end /Myr',3,3,2,log=[logval1,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
#     make_sub_scatter(data[var1],data[f'{name}_dur'],varlab1,f'{name}dynamo duration/Myr',3,3,3,log=[logval1,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
# #CIA dynamo
# name = 'cia'
# make_sub_scatter(data[var1],data[f'{name}_on'],varlab1,f'{name} dynamo start /Myr',3,3,4,log=[logval1,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
# make_sub_scatter(data[var1],data[f'{name}_off'],varlab1,f'{name} dynamo end /Myr',3,3,5,log=[logval1,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
# make_sub_scatter(data[var1],data[f'{name}_dur'],varlab1,f'{name}dynamo duration /Myr',3,3,6,log=[logval1,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
# #compositional dynamo
# name = 'comp'
# make_sub_scatter(data[var1],data[f'{name}_on'],varlab1,f'{name} dynamo start /Myr',3,3,7,log=[logval1,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
# make_sub_scatter(data[var1],data[f'{name}_off'],varlab1,f'{name} dynamo end /Myr',3,3,8,log=[logval1,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
# make_sub_scatter(data[var1],data[f'{name}_dur'],varlab1,f'{name}dynamo duration /Myr',3,3,9,log=[logval1,False,False],size=data['r'],sizelabel='radius /km',ss=1.5)
# if save == True:
#     plt.savefig(f'../Plots/{folder}dynamo_duration1_{var1}.png')

################## Intermittence of dynamo ###################################
#### var2 as x axis, var1 as colour
names = ['mac','cia','comp']
plt.figure(tight_layout=True,figsize=[15,15])
plt.suptitle(f'Effect of {varlab2} on dynamo intermittence')
for i, name in enumerate(names):
    if np.all(data[f'{name}_n']==0):
        pass #don't plot this dynamo - no intermittence
        print(f'No {name} intermittence')
    else:
        make_sub_scatter(data[var2],data[f'{name}_ngl10'],varlab2,f'{name} intermittence <10Myr spacing',3,3,3*i+1,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
        make_sub_scatter(data[var2],data[f'{name}_ngl100'],varlab2,f'{name} intermittence 10Myr < x < 100Myr spacing',3,3,3*i+2,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
        make_sub_scatter(data[var2],data[f'{name}_n']-(data[f'{name}_ngl10']+data[f'{name}_ngl100']),varlab2,f'{name} intermittence > 100Myr spacing',3,3,3*i+3,log=[logval2,False,logval1],colour=data[var1],colourlabel=varlab1,ss=5)
    i = i+1
if save == True:
    plt.savefig(f'../Plots/{folder}dynamo_intermit_{var2}{var1}_col.png')
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
    
###################### End of convection ########################################
plt.figure(tight_layout=True,figsize=[15,5])
plt.suptitle('Cessation of convection')
make_sub_scatter(data[var2], data['lconv_t'], varlab2,'Beginning of buffering /Myr', 1, 4, 1,colour=data[var1],colourlabel=varlab1,log=[logval2,False,logval1])
make_sub_scatter(data[var2], data['fcond_t'], varlab2, 'End of buffering /Myr', 1, 4, 2,colour=data[var1],colourlabel=varlab1,log=[logval2,False,logval1])
make_sub_scatter(data[var2], data['fcond_t']-data['lconv_t'], varlab2, 'Duration of buffering /Myr', 1, 4, 3,colour=data[var1],colourlabel=varlab1,log=[logval2,False,logval1])
make_sub_scatter(data[var2], data['lconv_T'], varlab2, 'temperature at end of convection /K', 1, 4, 4,colour=data[var1],colourlabel=varlab1,log=[logval2,False,logval1])
if save == True:
    plt.savefig(f'../Plots/{folder}conv_end_{var2}{var1}.png',dpi=450) 
    
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