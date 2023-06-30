#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 11:53:18 2023

@author: exet5460
"""
import pandas as pd
import matplotlib.pyplot as plt

file='../Results_combined/Timestep_test/timestep_test.csv'
save = False #decide whether to save
data = pd.read_csv(file,delimiter=',',skiprows=[1])
n =len(data)

#scale dt and dr by original step
data['dt']=data['dt']/data.loc[0,'dt']
data['dr']=data['dr']/500

#magnetic field generation as a fraction of total time of field generation
data['MAC_length']=data['MAC_stop']-data['MAC_start']
data['CIA_length']=data['CIA_stop']-data['CIA_start']
data['comp_length']=data['comp_stop']-data['comp_start']

#find mean and standard deviation of each column
mean = data.mean()
std = data.std()

#normalise data so can plot it all
data_norm = abs(data - mean) /mean
data_diff = data - mean
data_range = data.max()-data.min()
diff = data_range/mean

data_dt = data[data['dr']==1] #runs where only dt was changed
data_dt2 = data_dt[data_dt['dt']>0.45] #only up to half the original timestep
data_dr = data[data['dr']!=1] #runs where dr was changed so dt changed
data_dr = pd.concat([data_dr,data[data['run']==1]]) #add default run
data_norm_dr = (data_dr - data_dr.mean()) /data_dr.mean()
data_norm_dt = (data_dt - data_dt.mean()) /data_dt.mean()
data_norm_dt2 = (data_dt2 - data_dt2.mean()) /data_dt2.mean()

###################### Range in values/ mean ############################
drange = (data_dt2.max()-data_dt2.min())/data_dt2.mean()*100 #range as percentage of mean

plt.figure(figsize=[10,5])
plt.scatter(data_dt2.columns[3:],drange[3:])
plt.xticks(rotation=90)
#plt.yscale('log')
plt.ylabel('Range/mean /%')

################ Variation in mean for all data ########################
plt.figure(figsize=[10,5],tight_layout=True)
for col in data_norm_dt.columns[3:]:
    plt.scatter([col]*8,abs(data_norm_dt[col]*100),color='navy')
plt.xticks(rotation=90)
plt.yscale('log')
plt.hlines(1,'diff_time','comp_length',color='black',linestyle='--')
plt.fill_between(['diff_time','comp_length'],1e-4,1, alpha=0.3,color='black')
plt.ylabel('Deviation from mean %')
plt.title('Variation from the mean for variation in dt')
if save == True:
    plt.savefig('../Plots/timestep_test_dt.png',dpi=300)

plt.figure(figsize=[10,5],tight_layout=True)
for col in data_norm_dr.columns[3:]:
    plt.scatter([col]*4,abs(data_norm_dr[col]*100),color='navy')
plt.xticks(rotation=90)
plt.yscale('log')
plt.hlines(1,'diff_time','comp_length',color='black',linestyle='--')
plt.fill_between(['diff_time','comp_length'],1e-4,1, alpha=0.3,color='black')
plt.ylabel('Deviation from mean %')
plt.title('Variation from the mean for variation in dr')
if save == True:
    plt.savefig('../Plots/timestep_test_dr.png',dpi=300)

plt.figure(figsize=[10,5],tight_layout=True)
for col in data_norm_dt2.columns[3:]:
    plt.scatter([col]*4,abs(data_norm_dt2[col]*100),color='navy')
plt.xticks(rotation=90)
plt.yscale('log')
plt.hlines(1,'diff_time','comp_length',color='black',linestyle='--')
plt.fill_between(['diff_time','comp_length'],1e-4,1, alpha=0.3,color='black')
plt.ylabel('Deviation from mean %')
plt.title('Variation from the mean for dt$\\geq$0.5dt$_0$')
if save == True:
    plt.savefig('../Plots/timestep_test_dt05.png',dpi=300)
#################### Percentage difference plot, colour by dt ###########################
col1 = ['diff_T','diff_time','peakT','tmax','super_ad_end']
col2 = ['MAC_start','MAC_stop','CIA_start','CIA_stop','comp_start','comp_stop']
col3 = ['tstrat_remove','terode','super_ad_start','cond_t','max_Rem_mac','max_Remt_mac']
col4 = ['MAC_length','CIA_length','comp_length']

plt.figure(tight_layout=True,figsize=[10,10])
plt.suptitle('Deviation across change in dt')
plt.subplot(4,1,1)
for col in col1:
    plt.scatter([col]*n,data_norm[col]*100,label=f'{col}',c=data['dt'])
    plt.ylabel('deviation from mean %')
plt.colorbar(label='dt/dt$_0$',ticks=[0.25,0.5,0.64,1])

plt.subplot(4,1,2)
for col in col2:
    plt.scatter([col]*n,data_norm[col]*100,label=f'{col}',c=data['dt'])
    plt.ylabel('deviation from mean %')

plt.subplot(4,1,3)
for col in col3:
    plt.scatter([col]*n,data_norm[col]*100,label=f'{col}',c=data['dt'])
    plt.ylabel('deviation from mean %')

plt.subplot(4,1,4)
for col in col4:
    plt.scatter([col]*n,data_norm[col]*100,label=f'{col}',c=data['dt'])
    plt.ylabel('deviation from mean %')
if save == True:
    plt.savefig('Results_combined/Timestep_test/perc_diff.png',dpi=300)

############### Deviation for change in dt only ############################    
l = len(data_dt)
plt.figure(tight_layout=True,figsize=[10,10])
plt.suptitle('Deviation across change in dt,dr=500m')
plt.subplot(4,1,1)
for col in col1:
    plt.scatter([col]*l,data_norm_dt[col]*100,label=f'{col}',c=data_dt['dt'])
    plt.ylabel('deviation from mean %')
plt.colorbar(label='dt/dt$_0$',ticks=[0.25,0.5,1])

plt.subplot(4,1,2)
for col in col2:
    plt.scatter([col]*l,data_norm_dt[col]*100,label=f'{col}',c=data_dt['dt'])
    plt.ylabel('deviation from mean %')

plt.subplot(4,1,3)
for col in col3:
    plt.scatter([col]*l,data_norm_dt[col]*100,label=f'{col}',c=data_dt['dt'])
    plt.ylabel('deviation from mean %')

plt.subplot(4,1,4)
for col in col4:
    plt.scatter([col]*l,data_norm_dt[col]*100,label=f'{col}',c=data_dt['dt'])
    plt.ylabel('deviation from mean %')
if save == True:
    plt.savefig('Results_combined/Timestep_test/perc_diff_dt.png',dpi=300)

############### Deviation for change in dt only for dt>0.5dt0 ############################    
l = len(data_dt2)
tick_lab = [0.5,0.6,0.8,1]
plt.figure(tight_layout=True,figsize=[10,10])
plt.suptitle('Deviation across change in dt,dr=500m')
plt.subplot(4,1,1)
for col in col1:
    plt.scatter([col]*l,data_norm_dt2[col]*100,label=f'{col}',c=data_dt2['dt'])
    plt.ylabel('deviation from mean %')
plt.colorbar(label='dt/dt$_0$',ticks=tick_lab)

plt.subplot(4,1,2)
for col in col2:
    plt.scatter([col]*l,data_norm_dt2[col]*100,label=f'{col}',c=data_dt2['dt'])
    plt.ylabel('deviation from mean %')

plt.subplot(4,1,3)
for col in col3:
    plt.scatter([col]*l,data_norm_dt2[col]*100,label=f'{col}',c=data_dt2['dt'])
    plt.ylabel('deviation from mean %')

plt.subplot(4,1,4)
for col in col4:
    plt.scatter([col]*l,data_norm_dt2[col]*100,label=f'{col}',c=data_dt2['dt'])
    plt.ylabel('deviation from mean %')
if save == True:
    plt.savefig('Results_combined/Timestep_test/perc_diff_dt2.png',dpi=300)
#################### Percentage difference plot, colour by dr ###########################
m = len(data_dr)
plt.figure(tight_layout=True,figsize=[10,10])
plt.suptitle('Deviation across change in dr, dt=0.1$\\tau_{cond,core}$')
plt.subplot(4,1,1)
for col in col1:
    plt.scatter([col]*m,data_norm_dr[col]*100,label=f'{col}',c=data_dr['dr'])
    plt.ylabel('deviation from mean %')
plt.colorbar(label='dr/dr$_0$',ticks=[0.5,0.8,1])

plt.subplot(4,1,2)
for col in col2:
    plt.scatter([col]*m,data_norm_dr[col]*100,label=f'{col}',c=data_dr['dr'])
    plt.ylabel('deviation from mean %')

plt.subplot(4,1,3)
for col in col3:
    plt.scatter([col]*m,data_norm_dr[col]*100,label=f'{col}',c=data_dr['dr'])
    plt.ylabel('deviation from mean %')

plt.subplot(4,1,4)
for col in col4:
    plt.scatter([col]*m,data_norm_dr[col]*100,label=f'{col}',c=data_dr['dr'])
    plt.ylabel('deviation from mean %')
if save == True:
    plt.savefig('Results_combined/Timestep_test/perc_diff_dr.png',dpi=300)

######### Difference in magnetic field start and stop times as a fraction of total generation time #############
plt.figure()
for i, col in enumerate(col2):
    plt.scatter([col]*n,data_diff[col]*100/data[col4[int(i/2)]],label=f'{col}')
    plt.ylabel('t-$\\bar{t}$/$\\Delta$t %')
    
######################### Change in each value with timestep only #################
col5 = ['diff_T','diff_time','peakT','tmax']
col6 = ['MAC_start','MAC_stop','CIA_start','CIA_stop','comp_start','comp_stop','max_Rem_mac','max_Remt_mac']
col7 = ['tstrat_remove','terode','super_ad_start','super_ad_end','cond_t']
col8 = ['MAC_length','CIA_length','comp_length']

ylabels5 = ['T/K','t/Myr','T/K','t/Myr']
ylabels6 = ['t/Myr','t/Myr','t/Myr','t/Myr','t/Myr','t/Myr','','t/Myr']
ylabels7 = ['t/Myr','t/Myr','t/Myr','t/Myr','t/Myr']
ylabels8 = ['t/Myr','t/Myr','t/Myr']

plt.figure(tight_layout=True,figsize = [10,7])
plt.suptitle('Differentiation and peak temperature')
for i, col in enumerate(col5):
    plt.subplot(2,2,i+1)
    plt.scatter(data_dt['dt'],data_dt[col])
    plt.xlabel('dt/dt$_0$')
    plt.ylabel(ylabels5[i])
    plt.title(col)
    
plt.figure(tight_layout=True,figsize = [10,10])
plt.suptitle('Magnetic field generation')
for i, col in enumerate(col6):
    plt.subplot(4,2,i+1)
    plt.scatter(data_dt['dt'],data_dt[col])
    plt.xlabel('dt/dt$_0$')
    plt.ylabel(ylabels6[i])
    plt.title(col)
    
plt.figure(tight_layout=True,figsize = [10,10])
for i, col in enumerate(col7):
    plt.subplot(3,2,i+1)
    plt.scatter(data_dt['dt'],data_dt[col])
    plt.xlabel('dt/dt$_0$')
    plt.ylabel(ylabels7[i])
    plt.title(col)
    
plt.figure(tight_layout=True,figsize = [10,10])
for i, col in enumerate(col8):
    plt.subplot(3,1,i+1)
    plt.scatter(data_dt['dt'],data_dt[col])
    plt.xlabel('dt/dt$_0$')
    plt.ylabel(ylabels8[i])
    plt.title(col)
    
################ Change with timestep only - plots to indicate plateau #######

col6 = ['MAC_start','MAC_stop','CIA_start','max_Rem_mac']
col7 = ['tstrat_remove','super_ad_start','cond_t','MAC_length','CIA_length','comp_length']

ylabels6 = ['t/Myr','t/Myr','t/Myr','']
ylabels7 = ['t/Myr','t/Myr','t/Myr','t/Myr','t/Myr','t/Myr']
    
plt.figure(tight_layout=True,figsize = [15,10])
plt.suptitle('Magnetic field generation')
for i, col in enumerate(col6):
    plt.subplot(2,2,i+1)
    plt.scatter(data_dt['dt'],data_dt[col])
    plt.xlabel('dt/dt$_0$')
    plt.ylabel(ylabels6[i])
    plt.title(col)
if save == True:
    plt.savefig('../Plots/Timestep_test/mag_field_dt.png',dpi=600)
    
plt.figure(tight_layout=True,figsize = [15,15])
plt.suptitle('Cooling times and magnetic field duration')
for i, col in enumerate(col7):
    plt.subplot(3,2,i+1)
    plt.scatter(data_dt['dt'],data_dt[col])
    plt.xlabel('dt/dt$_0$')
    plt.ylabel(ylabels7[i])
    plt.title(col)
if save == True:
    plt.savefig('../Plots/Timestep_test/cooling_mag_dt.png',dpi=600)
    
######################### Change in each value with dr and timestep  #################
col5 = ['diff_T','diff_time','peakT','tmax']
col6 = ['MAC_start','MAC_stop','CIA_start','CIA_stop','comp_start','comp_stop','max_Rem_mac','max_Remt_mac']
col7 = ['tstrat_remove','terode','super_ad_start','super_ad_end','cond_t']
col8 = ['MAC_length','CIA_length','comp_length']

ylabels5 = ['T/K','t/Myr','T/K','t/Myr']
ylabels6 = ['t/Myr','t/Myr','t/Myr','t/Myr','t/Myr','t/Myr','','t/Myr']
ylabels7 = ['t/Myr','t/Myr','t/Myr','t/Myr','t/Myr']
ylabels8 = ['t/Myr','t/Myr','t/Myr']

plt.figure(tight_layout=True,figsize = [10,7])
plt.suptitle('Differentiation and peak temperature')
for i, col in enumerate(col5):
    plt.subplot(2,2,i+1)
    plt.scatter(data_dr['dr'],data_dr[col])
    plt.xlabel('dr /dr$_0$')
    plt.ylabel(ylabels5[i])
    plt.title(col)
    
plt.figure(tight_layout=True,figsize = [10,10])
plt.suptitle('Magnetic field generation')
for i, col in enumerate(col6):
    plt.subplot(4,2,i+1)
    plt.scatter(data_dr['dr'],data_dr[col])
    plt.xlabel('dr /dr$_0$')
    plt.ylabel(ylabels6[i])
    plt.title(col)
if save == True:
    plt.savefig('Results_combined/Timestep_test/mag_field_dr.png',dpi=300)
    
plt.figure(tight_layout=True,figsize = [10,10])
for i, col in enumerate(col7):
    plt.subplot(3,2,i+1)
    plt.scatter(data_dr['dr'],data_dr[col])
    plt.xlabel('dr /dr$_0$')
    plt.ylabel(ylabels7[i])
    plt.title(col)
    
plt.figure(tight_layout=True,figsize = [10,10])
for i, col in enumerate(col8):
    plt.subplot(3,1,i+1)
    plt.scatter(data_dr['dr'],data_dr[col])
    plt.xlabel('dr /dr$_0$')
    plt.ylabel(ylabels8[i])
    plt.title(col)