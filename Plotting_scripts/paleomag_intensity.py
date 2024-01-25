#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot the paleomagnetic record intensity vs time for the first 6Ma after CAIs
Only includes data from discussion section (omits IVAs and R chondrites)
Meteorite groups are grouped together
In Inkscape then add boxes for nebula and dynamo field
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%% Import data
paleo = pd.read_csv("../meteorite_paleomag_intensity.csv")
savefolder = 'EPSL_paper/'
save = True
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
#remove anything older than 12 Ma after CAI formation
paleo2 = paleo2[paleo2['midpoints']<12]
paleo2 = paleo2[paleo2['Classification']!='R matrix'] #remove R chondrites
mclass = paleo2['Classification'].unique()

#%% Intensity axis
fig, axes = plt.subplots(ncols=1,figsize=[7.5,5],tight_layout=True)
xlim_low = [0.8,12]
hwidth = 10 #headwidth
hlengths = 0.2 #head length - shorter in first plot

#colour dictionary
col = {'chondrule':'black','bulk':'#007aaf','QDM':'#007aaf','angrite':'#A00143'}
cs = 5 #capsize
for i, met in enumerate(mclass):
    mdata = paleo2.loc[paleo2['Classification']==met,:] #filter by class
    mdata.reset_index(inplace=True)
    
    for j in range(len(mdata)): #plot each value in class 
        if ((mdata.loc[j,'strength_lower']>0)&(mdata.loc[j,'strength_upper']>0))|(mdata.loc[j,'Meteorite']=='CR chondrules'): #if both plot as an error bar
            pcol = col[mdata.loc[j,'measured_phase']]
            if mdata.loc[j,'CRM_problem'] == True: #set data transparency
                trans = 0.3
            else:
                trans = 1
            axes.errorbar(mdata.loc[j,'midpoints'],mdata.loc[j,'strength_middle'],xerr=mdata.loc[j,'error'],
                          yerr=[[mdata.loc[j,'strength_middle']-mdata.loc[j,'strength_lower']],
                                [mdata.loc[j,'strength_upper']-mdata.loc[j,'strength_middle']]],color=pcol,
                          alpha = trans)
           
        elif (mdata.loc[j,'strength_lower']<=0)&(mdata.loc[j,'strength_middle']<=0): #if only upper limit
            pcol = col[mdata.loc[j,'measured_phase']]
            if mdata.loc[j,'CRM_problem'] == True: #set data transparency
                trans = 0.3
            else:
                trans = 1
            axes.errorbar(mdata.loc[j,'midpoints'],mdata.loc[j,'strength_upper'],xerr=mdata.loc[j,'error'],yerr=1,
                          uplims=True,color=pcol,alpha = trans)
        
        elif (mdata.loc[j,'Meteorite']=='Winchcombe')|(mdata.loc[j,'Meteorite']=='Kaba'): #only central value for these meteorite
            if mdata.loc[j,'CRM_problem'] == True: #set data transparency
                trans = 0.3
            else:
                trans = 1
            pcol = col[mdata.loc[j,'measured_phase']]
            axes.errorbar(mdata.loc[j,'midpoints'],mdata.loc[j,'strength_middle'],xerr=mdata.loc[j,'error'],
                          color=pcol,alpha = trans)
        
        elif (mdata.loc[j,'Meteorite']=='Kaba'): #all values are zero
            pcol = col[mdata.loc[j,'measured_phase']]
            if mdata.loc[j,'CRM_problem'] == True: #set data transparency
                trans = 0.3
            else:
                trans = 1
            axes.errorbar(mdata.loc[j,'midpoints'],mdata.loc[j,'strength_middle'],xerr=mdata.loc[j,'error'],
                          color=pcol,alpha = trans)
        
        elif (mdata.loc[j,'strength_upper']<=0)&(mdata.loc[j,'strength_middle']<=0): #if only lower limit
            pcol = col[mdata.loc[j,'measured_phase']]
            if mdata.loc[j,'CRM_problem'] == True: #set data transparency
                trans = 0.3
            else:
                trans = 1
            axes.errorbar(mdata.loc[j,'midpoints'],mdata.loc[j,'strength_upper'],xerr=mdata.loc[j,'error'],
                          lolims=True,color=pcol,alpha = trans)
        
        elif (mdata.loc[j,'strength_upper']<=0)&(mdata.loc[j,'strength_middle']!=0)&(mdata.loc[j,'strength_lower']!=0): #if only lower limit and middle value (add arrows in inkscape)
            pcol = col[mdata.loc[j,'measured_phase']]
            if mdata.loc[j,'CRM_problem'] == True: #set data transparency
                trans = 0.3
            else:
                trans = 1
            axes.errorbar(mdata.loc[j,'midpoints'],mdata.loc[j,'strength_middle'],xerr=mdata.loc[j,'error'],yerr=[[mdata.loc[j,'strength_middle']-mdata.loc[j,'strength_lower']],[10]],
                          color=pcol,alpha = trans)
        
        #add extra lines for no nebula correction
        if (mdata.loc[j,'nebula_correction']==True)&(mdata.loc[j,'measured_phase']=='QDM'): #Allende measurement
            axes.errorbar(mdata.loc[j,'midpoints'],mdata.loc[j,'strength_middle']/2,xerr=mdata.loc[j,'error'],
                          yerr=[[(mdata.loc[j,'strength_middle']-mdata.loc[j,'strength_lower'])/2],[10]],
                          color='grey')
        elif (mdata.loc[j,'Meteorite']=='Winchcombe'): #Winchcombe measurement
            axes.errorbar(mdata.loc[j,'midpoints'],mdata.loc[j,'strength_middle']/2,xerr=mdata.loc[j,'error'],
                          color='grey')
axes.fill_betweenx([0.1,150],1,3.8,color='#FFA07A',alpha=0.5) #nebula field
axes.fill_betweenx([0.1,60],4,12,color='#afe2c6') #previous dynamo
axes.fill_betweenx([0.1,60],1,4,color='#afe2c6',alpha=0.5) #current dynamo
axes.set_xlabel('Time /Ma')
axes.set_ylabel('Paleointensity /$\\mu$T')
axes.set_yscale('log')

if save == True:   
    plt.savefig(f'../Plots/{savefolder}paleomag_intensity_log.pdf',dpi=450,bbox_inches='tight')