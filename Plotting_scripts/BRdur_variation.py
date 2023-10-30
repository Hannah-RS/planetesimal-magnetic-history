#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plots for variation across a maxB, maxRem and duration for the varied parameters
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

#%% Set up folders etc.
folder = 'Paper_run2/'
subfolders = {'rcmf':1,'eta0':2,'beta':3,'etal':4,'Xs_0':5,'Fe0':6,'alpha_n':7,'r':8}
labels = {'rcmf':'$\\phi_{{RCMF}}$','eta0':'$\\eta_0$','beta':'$\\beta$','etal':'$\\eta_l$ ','Xs_0':'$X_{{s,0}}$','Fe0':'$^{{60}}Fe/^{{56}}Fe$','alpha_n':'$\\alpha_n$','r':'radius'}
units = {'rcmf':'','eta0':'Pas','beta':'$K^{-1}$','etal':'Pas','Xs_0':'wt %','Fe0':'','alpha_n':'','r':'km'}
logs =[False,True,False,True,False,True,False,False]
variables = ['rcmf','eta0','beta','etal','Xs_0','Fe0','alpha_n','r']
Myr = 365*24*3600*1e6 #number of s in Myr

#%% Calculate max and min values for all variables
var_results = pd.read_csv(f'../Results_combined/{folder}/all_sucess_info.csv',skiprows=[1])

#first row for min, second for max
B = np.zeros([2,len(variables)]) 
Rem = np.zeros([2,len(variables)])
dur = np.zeros([8,len(variables)]) #dynamo one, dynamo 2, gap, total duration


for i, var in enumerate(variables):

    path = '../Results_combined/'+folder+f"params_{subfolders[var]}/"
    
    #find run numbers
    var_data = pd.read_csv(path+'auto_params.csv',skiprows=[1])
    minrun = min(var_data['run'])
    maxrun = max(var_data['run'])
    nrun = len(var_data)
    #get relevant section of results
    data = var_results[(var_results['run']>=minrun)&(var_results['run']<=maxrun)].copy(deep=True)
    data.reset_index(inplace=True,drop=True)

    #get max and min values
    B[0,i]=data['max_B'].min()/1e-6 
    B[1,i]=data['max_B'].max()/1e-6
    
    Rem[0,i]=data['max_R'].min()
    Rem[1,i]=data['max_R'].max()
    
    #duration 
    dur1 = data['magoff_1']-data['magon_1']
    dur2 = data['magoff_2']-data['magon_2']
    totdur = dur1+dur2
    gap = (data['magon_2']-data['magoff_1'])[data['magon_2']>0]
    #gap[gap<0] = 0 #fill in no gap as 0
    dur[0,i]= np.min(dur1)
    dur[1,i]=np.max(dur1)
    dur[2,i]= np.min(dur2)
    dur[3,i]=np.max(dur2)
    dur[4,i]= np.min(gap)
    dur[5,i]=np.max(gap)
    dur[6,i]= np.min(totdur)
    dur[7,i]=np.max(totdur)
    
    
#%% Make bar plots
ticks = []
for var in variables:
    ticks.append(labels[var])
save = True   

#%% Rem and B
width = 0.3
fig = plt.figure()
ax = plt.axes()
ax2 = ax.twinx()
ln1 = ax.bar(variables,B[1,:],width=width,bottom=B[0,:],label='Max field strength',color='#4477AA')
ln2 = ax2.bar(variables,Rem[1,:],width=width/2,bottom=Rem[0,:],label='Max Reynolds number',color='#AA3377',align='edge')
ax.set_ylabel('Range of maximum field strengths /$\\mu$T')
ax2.set_ylabel('Range of maximum Rem')
ax.tick_params(axis='y',colors='#4477AA')
ax2.tick_params(axis='y',colors='#AA3377')
ax.set_xlabel('Variable')
ax.set_xticks(variables,ticks)
if save == True:
    plt.savefig(f'../Plots/{folder}/RemB_summary.png',bbox_inches='tight')
    
#%% Duration - for each dynamo and gap, maybe too confusing?
width=0.75
x = 2*np.arange(len(variables))
offset = width/3

fig = plt.figure()
ax = plt.axes()
ax2 = ax.twinx()
ln1 = ax.bar(x-offset,dur[1,:],width=width/3,bottom=dur[0,:],label='First dynamo',color='#4477AA')
ln2 = ax.bar(x,dur[3,:],width=width/3,bottom=dur[2,:],label='Second dynamo',color='#66CCEE')
ln3 = ax2.bar(x+offset,dur[5,:],width=width/3,bottom=dur[4,:],label='Gap in dynamo generation',color='#AA3377')
ax.set_ylabel('Range of dynamo durations/Myr')
ax2.set_ylabel('Range of gaps in generation /Myr')
ax2.set_ylim(bottom=0)
ax.tick_params(axis='y',colors='#4477AA')
ax2.tick_params(axis='y',colors='#AA3377')
ax.set_xlabel('Variable')
ax.set_xticks(x,ticks)
ax.legend(framealpha=1,bbox_to_anchor=[0.7,-0.2])
if save == True:
    plt.savefig(f'../Plots/{folder}/duration_summary.png',bbox_inches='tight')
    
#%% Simpler duration
width=0.75
x = 2*np.arange(len(variables))
offset = width/2

fig = plt.figure()
ax = plt.axes()
ax2 = ax.twinx()
ln1 = ax.bar(x-offset,dur[7,:],width=width/2,bottom=dur[6,:],label='First dynamo',color='#4477AA')
ln3 = ax2.bar(x+offset,dur[5,:],width=width/2,bottom=dur[4,:],label='Gap in dynamo generation',color='#AA3377')
ax.set_ylabel('Range of dynamo durations/Myr')
ax2.set_ylabel('Range of gaps in generation /Myr')
ax.tick_params(axis='y',colors='#4477AA')
ax2.tick_params(axis='y',colors='#AA3377')
ax2.set_ylim(bottom=0)
ax.set_xlabel('Variable')
ax.set_xticks(x,ticks)
ax2.legend([ln1,ln3],['Total duration','Gap between dynamos'],framealpha=1,bbox_to_anchor=[0.7,-0.2])
if save == True:
    plt.savefig(f'../Plots/{folder}/duration_summary_simple.png',bbox_inches='tight')