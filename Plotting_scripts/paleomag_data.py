#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot the paleomagnetic record
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%% Import data
paleo = pd.read_csv("../meteorite_paleomagnetism.csv",skiprows=[1])

#replace nan values in timing info
values = {"rel_age_lower": 0, "rel_age_upper": 0}
paleo2 = paleo.fillna(value=values)
#remove data with no timing info
paleo2 = paleo2[(paleo2['rel_age_lower']!=0)|(paleo2['rel_age_upper']!=0)]
#calculate midpoints
paleo2.loc[:,'midpoints'] = (paleo2.loc[:,'rel_age_lower']+paleo2.loc[:,'rel_age_upper'])/2
#calculate error
paleo2.loc[:,'error'] = paleo2.loc[:,'midpoints']-paleo2.loc[:,'rel_age_lower']
#%% Process data
mclass = paleo2['Classification'].unique()
yplot = np.arange(2*len(mclass)+1) #value to plot a class against
midpoints = (yplot[1:]+yplot[:-1])/2
fig, axes = plt.subplots(ncols=2,sharey='row',gridspec_kw={'width_ratios': [3,1]},figsize=[15,5],tight_layout=True)
xlim = [300,12] #x axis limits
title = ['Full record',f'First {xlim[1]} Myr']
for ax, xlim, title in zip(axes,xlim,title):
    #plot all classes
    for i, met in enumerate(mclass):
        mdata = paleo2.loc[paleo2['Classification']==met,:] #filter by class
        mdata.reset_index(inplace=True)
     
        for j in range(len(mdata)): #plot each value in class - probably want this to be fill between 
            if mdata.loc[j,'magnetised']==True:
                if mdata.loc[j,'radiometric']==True:
                    fcol = '#5b93eb'
                    ecol = '#0000af'
                    hatch = ''
                else: 
                    fcol = '#b5b5b1'
                    ecol = '#2a2728'
                    hatch = ''
            else:
                if mdata.loc[j,'radiometric']==True:
                    fcol = 'white'
                    ecol = '#0000af'
                    hatch = '/'
                else:
                    fcol = 'white'
                    ecol = '#2a2728'
                    hatch = '/'
            if (mdata.loc[j,'rel_age_lower']!=0)&(mdata.loc[j,'rel_age_upper']!=0): #if range of values fill between
                ax.fill_betweenx(yplot[2*i:(2*i+2)],mdata.loc[j,'rel_age_lower'],mdata.loc[j,'rel_age_upper'],alpha=0.5,facecolor=fcol,edgecolor=ecol,hatch=hatch)
            
            elif mdata.loc[j,'rel_age_lower']!=0: #plot single point for now
                ax.scatter(mdata.loc[j,'rel_age_lower'],yplot[2*i],marker='|',color=fcol)
                ax.arrow(mdata.loc[j,'rel_age_lower'],yplot[2*i],mdata.loc[j,'rel_age_lower']*0.3,0,head_width=0.5,head_length=mdata.loc[j,'rel_age_lower']*0.05,capstyle='butt',edgecolor=ecol,facecolor=fcol)
                if xlim > 100: #for big plot
                    ax.annotate(' > 870 Myr',(xlim-10,midpoints[2*i]),(xlim-50,midpoints[2*i]),arrowprops=dict(arrowstyle=']->',facecolor=fcol,edgecolor=ecol))
            elif mdata.loc[j,'rel_age_upper']!=0:
                ax.scatter(mdata.loc[j,'rel_age_upper'],midpoints[2*i],marker='|',color=fcol)
                ax.arrow(mdata.loc[j,'rel_age_upper'],midpoints[2*i],-mdata.loc[j,'rel_age_upper']*0.3,0,head_width=0.5,head_length=mdata.loc[j,'rel_age_upper']*0.05,capstyle='butt',edgecolor=ecol,facecolor=fcol)
            else:
                raise ValueError('Both values are 0')
    #overall figure things
    ax.set_yticks(midpoints[::2],mclass)
    ax.set_xlabel('Time after CAIs /Myr')
    ax.set_xlim([0,xlim])
    ax.set_title(title)
    
plt.savefig('../Plots/CoS/paleomag_record.pdf',dpi=450,bbox_inches='tight')

#%% Errorbar version
fig, axes = plt.subplots(ncols=2,sharey='row',gridspec_kw={'width_ratios': [3,1]},figsize=[15,5],tight_layout=True)
xlim = [300,12] #x axis limits
title = ['Full record',f'First {xlim[1]} Myr']

for ax, xlim, title in zip(axes,xlim,title):
    #plot all classes
    for i, met in enumerate(mclass):
        mdata = paleo2.loc[paleo2['Classification']==met,:] #filter by class
        mdata.reset_index(inplace=True)
        
        for j in range(len(mdata)): #plot each value in class - probably want this to be fill between 
            if mdata.loc[j,'magnetised']==True:
                if mdata.loc[j,'radiometric']==True:
                    fcol = '#5b93eb'
                    ecol = '#0000af'
                    linestyle = '-'
                else: 
                    fcol = '#b5b5b1'
                    ecol = '#2a2728'
                    linestyle = '--'
            else:
                if mdata.loc[j,'radiometric']==True:
                    fcol = 'white'
                    ecol = '#0000af'
                    hatch = '/'
                else:
                    fcol = 'white'
                    ecol = '#2a2728'
                    hatch = '/'
            
            if (mdata.loc[j,'rel_age_lower']!=0)&(mdata.loc[j,'rel_age_upper']!=0): #if range of values fill between
                if (mdata.loc[j,'up_age_source']=='both'): #if upper and lower bounds use different methods
                    #lower limit is radiometric
                    ax.errorbar(mdata.loc[j,'midpoints'],yplot[2*i],xerr=[[mdata.loc[j,'error']],[0]],markerfacecolor=fcol,markeredgecolor=ecol,ecolor=ecol,capsize=5,marker='o')
                    #upper limit from modelling
                    ax.errorbar(mdata.loc[j,'midpoints'],yplot[2*i],xerr=[[0],[mdata.loc[j,'error']]],ecolor='#2a2728',capsize=5)
                else:
                    ax.errorbar(mdata.loc[j,'midpoints'],yplot[2*i],xerr=mdata.loc[j,'error'],markerfacecolor=fcol,markeredgecolor=ecol,ecolor=ecol,capsize=5,marker='o')
            
            elif mdata.loc[j,'rel_age_lower']!=0: #plot single point for now
                ax.scatter(mdata.loc[j,'rel_age_lower'],yplot[2*i],marker='|',color=fcol)
                ax.arrow(mdata.loc[j,'rel_age_lower'],yplot[2*i],mdata.loc[j,'rel_age_lower']*0.3,0,head_width=0.5,head_length=mdata.loc[j,'rel_age_lower']*0.05,capstyle='butt',edgecolor=ecol,facecolor=fcol)
                if xlim > 100: #for big plot
                    ax.annotate(' > 870 Myr',(xlim-10,midpoints[2*i]),(xlim-50,midpoints[2*i]),arrowprops=dict(arrowstyle=']->',facecolor=fcol,edgecolor=ecol))
            elif mdata.loc[j,'rel_age_upper']!=0:
                ax.scatter(mdata.loc[j,'rel_age_upper'],midpoints[2*i],marker='|',color=fcol)
                ax.arrow(mdata.loc[j,'rel_age_upper'],midpoints[2*i],-mdata.loc[j,'rel_age_upper']*0.3,0,head_width=0.5,head_length=mdata.loc[j,'rel_age_upper']*0.05,capstyle='butt',edgecolor=ecol,facecolor=fcol)
            else:
                raise ValueError('Both values are 0')
    #overall figure things
    ax.set_yticks(midpoints[::2],mclass)
    ax.set_xlabel('Time after CAIs /Myr')
    ax.set_xlim([0,xlim])
    ax.set_title(title)
    
plt.savefig('../Plots/CoS/paleomag_record_dots.png',dpi=450,bbox_inches='tight')