#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot the paleomagnetic record
In Inkscape then convert non-magnetised values to dashes
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

#%% Arrow version
fig, axes = plt.subplots(ncols=2,sharey='row',gridspec_kw={'width_ratios': [1,3]},figsize=[15,5],tight_layout=True)
xlim_up = [12,300] #x axis upper limits
xlim_low = [0.8,12]
title = [f'First {xlim_up[0]} Myr',f'{xlim_low[1]} to {xlim_up[1]} Myr']
k = 0 #enumerates through axes
hwidth = 1.5 #headwidth
hlengths = [0.2,2] #head length - shorter in first plot
for ax, xlim_up, xlim_low, title, hlength in zip(axes,xlim_up,xlim_low,title,hlengths):
    #plot all classes
    for i, met in enumerate(mclass):
        mdata = paleo2.loc[paleo2['Classification']==met,:] #filter by class
        mdata.reset_index(inplace=True)
        
        for j in range(len(mdata)): #plot each value in class - probably want this to be fill between 
            if mdata.loc[j,'magnetised']==True:
                if mdata.loc[j,'radiometric']==True:
                    fcol = '#007aaf'#'#5b93eb'
                    ecol = '#007aaf'
                    ls = 'solid'
                else: 
                    fcol = '#2a2728'#'#b5b5b1'
                    ecol = '#2a2728'
                    ls = 'solid'
            else:
                if mdata.loc[j,'radiometric']==True:
                    fcol = 'white'
                    ecol = '#007aaf'
                    ls = (5,(3,2))
                else:
                    fcol = 'white'
                    ecol = '#2a2728'
                    ls = (3.5,(2,1.5))
            
            if (mdata.loc[j,'rel_age_lower']!=0)&(mdata.loc[j,'rel_age_upper']!=0): #if range of values fill between
                if (mdata.loc[j,'up_age_source']=='both'): #if upper and lower bounds use different methods
                    #lower limit is radiometric
                    ax.arrow(mdata.loc[j,'midpoints'],midpoints[2*i],mdata.loc[j,'error'],0,head_width=hwidth,head_length=hlength,edgecolor=ecol,facecolor=fcol,linestyle=ls,linewidth=hwidth)
                    #upper limit from modelling
                    ax.arrow(mdata.loc[j,'midpoints'],midpoints[2*i],-mdata.loc[j,'error'],0,head_width=hwidth,head_length=hlength,edgecolor='#2a2728',facecolor=fcol,linestyle=ls,linewidth=hwidth)
                else: #two overlapping arrows of same color
                    ax.arrow(mdata.loc[j,'midpoints'],midpoints[2*i],mdata.loc[j,'error'],0,head_width=hwidth,head_length=hlength,edgecolor=ecol,facecolor=fcol,linestyle=ls,linewidth=hwidth)
                    ax.arrow(mdata.loc[j,'midpoints'],midpoints[2*i],-mdata.loc[j,'error'],0,head_width=hwidth,head_length=hlength,edgecolor=ecol,facecolor=fcol,linestyle=ls,linewidth=hwidth)

            elif (mdata.loc[j,'rel_age_upper']==0)&(mdata.loc[j,'rel_age_lower']!=0): #plot single point for now
                ax.scatter(mdata.loc[j,'rel_age_lower'],midpoints[2*i],marker='|',color=fcol)
                ax.arrow(mdata.loc[j,'rel_age_lower'],midpoints[2*i],mdata.loc[j,'rel_age_lower']*0.3,0,head_width=hwidth,head_length=hlength,capstyle='butt',edgecolor=ecol,facecolor=fcol)
                if xlim_up > 100: #for big plot
                    ax.annotate(' > 870 Myr',(xlim_up-10,midpoints[2*i]),(xlim_up-50,midpoints[2*i]),arrowprops=dict(arrowstyle=']->',facecolor=fcol,edgecolor=ecol))
            elif (mdata.loc[j,'rel_age_upper']!=0)&(mdata.loc[j,'rel_age_lower']==0): #extend arrow all the way to zero
                ax.scatter(mdata.loc[j,'rel_age_upper'],midpoints[2*i],marker='|',color=fcol)
                ax.arrow(mdata.loc[j,'rel_age_upper'],midpoints[2*i],-mdata.loc[j,'rel_age_upper']+1,0,head_width=hwidth,head_length=hlength,capstyle='butt',edgecolor=ecol,facecolor=fcol)
            else:
                raise ValueError('Both values are 0')
    
    #overall figure things
    ax.set_yticks(midpoints[::2],mclass)
    ax.set_xlabel('Time /Ma')
    ax.set_xlim([xlim_low,xlim_up])
    ax.set_title(title)
    
#plt.savefig('../Plots/CoS/paleomag_record.pdf',dpi=450,bbox_inches='tight')