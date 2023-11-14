#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create the filled timings plot for all variables
N.B. It is long because of complicated label positioning
"""
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

#%% Set up folders etc.
folder = 'Paper_run4/'
subfolders = {'rcmf':1,'eta0':2,'beta':3,'etal':4,'Xs_0':5,'Fe0':6,'alpha_n':7,'r':8}
labels = {'rcmf':'$\\phi_{{RCMF}}$','eta0':'$\\eta_0$','beta':'$\\beta$','etal':'$\\eta_l$ ','Xs_0':'$X_{{s,0}}$','Fe0':'$^{{60}}Fe/^{{56}}Fe$','alpha_n':'$\\alpha_n$','r':'radius'}
units = {'rcmf':'','eta0':'Pas','beta':'$K^{-1}$','etal':'Pas','Xs_0':'wt %','Fe0':'','alpha_n':'','r':'km'}
logs =[False,True,False,True,False,True,False,False]
Myr = 365*24*3600*1e6 #number of s in Myr
save = True

#%% Load a given variable
variables = ['rcmf','eta0','beta','etal','Xs_0','Fe0','alpha_n','r']

for i, var in enumerate(variables):

    unit = units[var]
    varlab = labels[var]
    logvar = logs[i]
    path = '../Results_combined/'+folder+f"params_{subfolders[var]}/"
    
    #find run numbers
    var_data = pd.read_csv(path+'auto_params.csv',skiprows=[1])
    var_results = pd.read_csv(f'../Results_combined/{folder}/all_sucess_info.csv',skiprows=[1])
    minrun = min(var_data['run'])
    maxrun = max(var_data['run'])
    nrun = len(var_data)
    data = var_results[(var_results['run']>=minrun)&(var_results['run']<=maxrun)].copy(deep=True)
    data.reset_index(inplace=True,drop=True)
    
    if var == 'r':
        data[var] = data[var]/1e3 #convert to km
        data2 = data[data['magon_1']>0] #remove off dynamo generation periods

    if var == 'Fe0':
        data.loc[data['Fe0']==0,'Fe0']=1e-10
        
    data['fcond_t']=data['fcond_t'].fillna(data['tsolid']) #if mantle doesnt stop convecting fill until end of thermal evolution
    data1 = data[data['Bn1']>1] #filter by number of generation periods
    
    #%% Colour options for all plots
    dcol = 'gray' #undifferentiated
    pmcol = 'firebrick' #peak mantle temp
    pccol = 'dimgray' #peak core temp
    scol = 'yellow' #stratification 
    srcol ='yellowgreen' #stratification being removed
    bcol = 'cornflowerblue' #magnetic field
    mccol = 'coral' #convecting mantle
    cscol = 'mediumpurple' #core solidifying
    ecol ='gray' #edgecolor for generic label box
    
    asval = 0.5 #alpha for shading
    alval=0.3 #alpha for label box
    
    #hatching
    chatch = '/' #solidifying core
    shatch = '*' #stratified core
    bhatch = '.' #magnetic field
    

    #%%make the plot
    with sns.plotting_context('talk'):
        plt.figure(figsize=[20,7.5])
        plt.title(f'Variation in thermal history and dynamo timings \n as a function of {varlab}')
        plt.plot(data['tmax'],data[var],label='peak mantle temp',color=pmcol,alpha=asval)
        plt.plot(data['tcoremax'],data[var],label='peak core temp',color=pccol,alpha=asval)
        plt.fill_betweenx(data[var],0.8,data['diff_time'],color=dcol,alpha=asval)
        plt.fill_betweenx(data[var],data['tstrat_start'],data['tstrat_remove'],color=scol,alpha=asval,hatch=shatch,edgecolor='lightgrey')
        plt.fill_betweenx(data[var],data['tstrat_remove'],data['terode'],color=srcol,alpha=asval,hatch=shatch)
        if var =='r':
            plt.fill_betweenx(data2[var],data2['magon_1'],data2['magoff_1'],color=bcol,alpha=alval,hatch=bhatch)
        else:
            plt.fill_betweenx(data[var],data['magon_1'],data['magoff_1'],color=bcol,alpha=alval,hatch=bhatch)
        plt.fill_betweenx(data1[var],data1['magon_2'],data1['magoff_2'],color=bcol,alpha=alval,hatch=bhatch)
        plt.fill_betweenx(data[var],data['diff_time'],data['fcond_t'],color=mccol,alpha=asval)
        plt.fill_betweenx(data[var],data['tsolid_start'],data['tsolid'],color=cscol,alpha=alval,hatch=chatch)
        plt.xlabel('Time/Myr')
        plt.xscale('log')
        plt.ylabel(f'{varlab}/ {unit}')
        if logvar == True:
            plt.yscale('log')
        
        #specific variable labels
        if var == 'rcmf':
            plt.text(x=1,y=0.3,s='Undifferentiated',rotation='vertical',bbox=dict(edgecolor=ecol,facecolor=dcol,alpha=alval))
            plt.annotate('Peak mantle temp',(data.loc[2,'tmax'],0.4),(data.loc[2,'tmax']+1.5,0.4),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=pmcol,alpha=alval))
            plt.annotate('Peak core temp',(data.loc[0,'tcoremax'],0.32),(data.loc[0,'tcoremax']+1,0.32),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=pccol,alpha=alval))
            plt.annotate('Erosion of core \n thermal stratification',(data.loc[2,'terode']-0.1,0.27),(data.loc[2,'terode']+1,0.27),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=srcol,alpha=asval))
            plt.text(x=1.4,y=0.35,s='Core \n thermally \n stratified',bbox=dict(edgecolor=ecol,facecolor=scol,alpha=asval))
            plt.text(x=data.loc[1,'magon_1']*5,y=0.29,s='Dynamo on',bbox=dict(edgecolor=ecol,facecolor=bcol,alpha=alval,hatch=bhatch))
            plt.text(x=110,y=0.3,s='Dynamo on',bbox=dict(edgecolor=ecol,facecolor=bcol,alpha=alval,hatch=bhatch))
            plt.text(x=data.loc[1,'fcond_t']*0.3,y=0.35,s='Mantle convecting',bbox=dict(edgecolor=ecol,facecolor=mccol))
            plt.text(x=140,y=min(data[var])*1.2,s='Core \n solidifying',bbox=dict(edgecolor=ecol,facecolor=cscol,alpha=alval,hatch=chatch))

        if var == 'eta0':
            plt.text(x=1,y=max(data[var])*1e-6,s='Undifferentiated',rotation='vertical',bbox=dict(edgecolor=ecol,facecolor=dcol,alpha=alval))
            plt.annotate('Peak mantle temp',(data.loc[10,'tmax'],1e23),(data.loc[10,'tmax']+1.5,1e23),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=pmcol,alpha=alval))
            plt.annotate('Peak core temp',(data.loc[1,'tcoremax'],1e15),(data.loc[1,'tcoremax']+1,1e15),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=pccol,alpha=alval))
            plt.annotate('Erosion of core \n thermal stratification',(data.loc[2,'terode'],max(data[var])*1e-3),(data.loc[2,'terode']+1,max(data[var])*5e-4),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=srcol,alpha=asval))
            plt.text(x=1.4,y=1e18,s='Core \n thermally \n stratified',bbox=dict(edgecolor=ecol,facecolor=scol,alpha=asval))
            plt.text(x=data.loc[2,'magon_1']*4,y=max(data[var])*1e-5,s='Dynamo on',bbox=dict(edgecolor=ecol,facecolor=bcol,alpha=alval,hatch=bhatch))
            plt.text(x=110,y=1e21,s='Dynamo on',bbox=dict(edgecolor=ecol,facecolor=bcol,alpha=alval,hatch=bhatch))
            plt.text(x=20,y=max(data[var])*1e-2,s='Mantle convecting',bbox=dict(edgecolor=ecol,facecolor=mccol))
            plt.text(x=100,y=1e16,s='Core \n solidifying',bbox=dict(edgecolor=ecol,facecolor=cscol,alpha=alval,hatch=chatch))
            
        if var == 'beta':
            plt.text(x=1,y=max(data[var]/2),s='Undifferentiated',rotation='vertical',bbox=dict(edgecolor=ecol,facecolor=dcol,alpha=alval))
            plt.annotate('Peak mantle temp',(data.loc[2,'tmax'],0.03),(data.loc[2,'tmax']+1.5,0.032),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=pmcol,alpha=alval))
            plt.annotate('Peak core temp',(data.loc[0,'tcoremax'],0.012),(data.loc[0,'tcoremax']+1,0.012),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=pccol,alpha=alval))
            plt.annotate('Erosion of core \n thermal stratification',(data.loc[2,'tstrat_remove']+0.1,max(data[var])*0.5),(data.loc[2,'terode']+1,max(data[var])*0.47),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=srcol,alpha=asval))
            plt.text(x=1.4,y=0.0225,s='Core \n thermally \n stratified',bbox=dict(edgecolor=ecol,facecolor=scol,alpha=asval))
            plt.text(x=10,y=0.03,s='Dynamo on',bbox=dict(edgecolor=ecol,facecolor=bcol,alpha=alval,hatch=bhatch))
            plt.text(x=110,y=0.03,s='Dynamo on',bbox=dict(edgecolor=ecol,facecolor=bcol,alpha=alval,hatch=bhatch))
            plt.text(x=10,y=0.0225,s='Mantle convecting',bbox=dict(edgecolor=ecol,facecolor=mccol))
            plt.text(x=140,y=min(data[var])*1.2,s='Core \n solidifying',bbox=dict(edgecolor=ecol,facecolor=cscol,alpha=alval,hatch=chatch))
       
        if var == 'etal':
            plt.text(x=1,y=1,s='Undifferentiated',rotation='vertical',bbox=dict(edgecolor=ecol,facecolor=dcol,alpha=alval))
            plt.annotate('Peak mantle temp',(data.loc[2,'tmax'],6),(data.loc[2,'tmax']+1.5,6),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=pmcol,alpha=alval))
            plt.annotate('Peak core temp',(data.loc[0,'tcoremax'],3),(data.loc[0,'tcoremax']+1,3),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=pccol,alpha=alval))
            plt.annotate('Erosion of core \n thermal stratification',(1.9,10),(3,12),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=srcol,alpha=asval))
            plt.text(x=1.3,y=25,s='Core \n thermally \n stratified',bbox=dict(edgecolor=ecol,facecolor=scol,alpha=asval))
            plt.text(x=20,y=50,s='Dynamo on',bbox=dict(edgecolor=ecol,facecolor=bcol,alpha=alval,hatch=bhatch))
            plt.text(x=110,y=50,s='Dynamo on',bbox=dict(edgecolor=ecol,facecolor=bcol,alpha=alval,hatch=bhatch))
            plt.text(x=20,y=20,s='Mantle convecting',bbox=dict(edgecolor=ecol,facecolor=mccol))
            plt.text(x=140,y=10,s='Core \n solidifying',bbox=dict(edgecolor=ecol,facecolor=cscol,alpha=alval,hatch=chatch))
            
        if var == 'Xs_0':
            plt.text(x=1,y=29,s='Undifferentiated',rotation='vertical',bbox=dict(edgecolor=ecol,facecolor=dcol,alpha=alval))
            plt.annotate('Peak mantle temp',(data.loc[2,'tmax'],31),(data.loc[2,'tmax']+1.5,31),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=pmcol,alpha=alval))
            plt.annotate('Peak core temp',(data.loc[0,'tcoremax'],28),(data.loc[0,'tcoremax']+1,28),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=pccol,alpha=alval))
            plt.annotate('Erosion of core \n thermal stratification',(1.9,28.5),(4,29),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=srcol,alpha=asval))
            plt.text(x=1.25,y=29,s='Core \n thermally \n stratified',bbox=dict(edgecolor=ecol,facecolor=scol,alpha=asval))
            plt.text(x=20,y=31,s='Dynamo on',bbox=dict(edgecolor=ecol,facecolor=bcol,alpha=alval,hatch=bhatch))
            plt.text(x=110,y=31,s='Dynamo on',bbox=dict(edgecolor=ecol,facecolor=bcol,alpha=alval,hatch=bhatch))
            plt.text(x=20,y=29,s='Mantle convecting',bbox=dict(edgecolor=ecol,facecolor=mccol))
            plt.text(x=140,y=28,s='Core \n solidifying',bbox=dict(edgecolor=ecol,facecolor=cscol,alpha=alval,hatch=chatch))
        
        if var == 'alpha_n':
            plt.text(x=1,y=32.5,s='Undifferentiated',rotation='vertical',bbox=dict(edgecolor=ecol,facecolor=dcol,alpha=alval))
            plt.annotate('Peak mantle temp',(data.loc[2,'tmax'],42),(data.loc[2,'tmax']+1.5,42),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=pmcol,alpha=alval))
            plt.annotate('Peak core temp',(data.loc[0,'tcoremax'],27),(data.loc[0,'tcoremax']+1,27),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=pccol,alpha=alval))
            plt.annotate('Erosion of core \n thermal stratification',(1.8,31),(3,31),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=srcol,alpha=asval))
            plt.text(x=1.25,y=35,s='Core \n thermally \n stratified',bbox=dict(edgecolor=ecol,facecolor=scol,alpha=asval))
            plt.text(x=20,y=40,s='Dynamo on',bbox=dict(edgecolor=ecol,facecolor=bcol,alpha=alval,hatch=bhatch))
            plt.text(x=100,y=40,s='Dynamo on',bbox=dict(edgecolor=ecol,facecolor=bcol,alpha=alval,hatch=bhatch))
            plt.text(x=20,y=35,s='Mantle convecting',bbox=dict(edgecolor=ecol,facecolor=mccol))
            plt.text(x=140,y=27,s='Core \n solidifying',bbox=dict(edgecolor=ecol,facecolor=cscol,alpha=alval,hatch=chatch))

        if var == 'Fe0':
            plt.text(x=1,y=1e-9,s='Undifferentiated',rotation='vertical',bbox=dict(edgecolor=ecol,facecolor=dcol,alpha=alval))
            plt.annotate('Peak mantle temp',(data.loc[4,'tmax'],1e-7),(data.loc[4,'tmax']+1.5,1e-7),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=pmcol,alpha=alval))
            plt.annotate('Peak core temp',(data.loc[2,'tcoremax'],1e-8),(data.loc[2,'tcoremax']+1,1e-8),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=pccol,alpha=alval))
            plt.annotate('Erosion of core \n thermal stratification',(5,1e-9),(12,3e-10),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=srcol,alpha=asval))
            plt.text(x=1.5,y=5e-10,s='Core \n thermally \n stratified',bbox=dict(edgecolor=ecol,facecolor=scol,alpha=asval))
            plt.text(x=20,y=1e-9,s='Dynamo on',bbox=dict(edgecolor=ecol,facecolor=bcol,alpha=alval,hatch=bhatch))
            plt.text(x=110,y=1e-9,s='Dynamo on',bbox=dict(edgecolor=ecol,facecolor=bcol,alpha=alval,hatch=bhatch))
            plt.text(x=20,y=max(data[var])*1e-2,s='Mantle convecting',bbox=dict(edgecolor=ecol,facecolor=mccol))
            plt.text(x=140,y=1e-7,s='Core \n solidifying',bbox=dict(edgecolor=ecol,facecolor=cscol,alpha=alval,hatch=chatch))
            plt.yticks(ticks=[1e-10,1e-9,1e-8,1e-7,6e-7],labels=[0,'$10^{-9}$','$10^{-8}$','$10^{-7}$','$6\\times10^{-7}$'])
    
        if var == 'r':
            plt.text(x=1,y=200,s='Undifferentiated',rotation='vertical',bbox=dict(edgecolor=ecol,facecolor=dcol,alpha=alval))
            plt.annotate('Peak mantle temp',(data.loc[3,'tmax'],470),(data.loc[3,'tmax']+1.5,470),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=pmcol,alpha=alval))
            plt.annotate('Peak core temp',(data.loc[0,'tcoremax'],110),(data.loc[0,'tcoremax']+1,110),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=pccol,alpha=alval))
            plt.annotate('Erosion of core \n thermal stratification',(2,200),(3,200),arrowprops=dict(facecolor='black',edgecolor='black'),bbox=dict(edgecolor=ecol,facecolor=srcol,alpha=asval))
            plt.text(x=1.25,y=300,s='Core \n thermally \n stratified',bbox=dict(edgecolor=ecol,facecolor=scol,alpha=asval))
            plt.text(x=20,y=400,s='Dynamo on',bbox=dict(edgecolor=ecol,facecolor=bcol,alpha=alval,hatch=bhatch))
            plt.text(x=160,y=400,s='Dynamo on',bbox=dict(edgecolor=ecol,facecolor=bcol,alpha=alval,hatch=bhatch))
            plt.text(x=20,y=250,s='Mantle convecting',bbox=dict(edgecolor=ecol,facecolor=mccol))
            plt.text(x=50,y=150,s='Core \n solidifying',bbox=dict(edgecolor=ecol,facecolor=cscol,alpha=alval,hatch=chatch))

    

    if save == True:
        plt.savefig(f'../Plots/{folder}all_timings_{var}.png',dpi=450,bbox_inches='tight') 