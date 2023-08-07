#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Intermittence analysis - investigate which runs have multiple periods of dynamo generation
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

folder = 'Minirun3' #original run
folder2='Minirun4' #run with smaller solidification timestep
data = pd.read_csv(f'../Results_combined/{folder}/all_sucess_info.csv',delimiter=',',skiprows=[1],header=0,index_col=False)
data_s = pd.read_csv(f'../Results_combined/{folder2}/all_sucess_info.csv',delimiter=',',skiprows=[1],header=0,index_col=False)
data_s['comp_diff']=data['comp_n']-data_s['comp_n']
data['r']=data['r']/1e3 #rescale to km
data_s['r']=data_s['r']/1e3 #rescale to km
#filter for more than 1 dynamo generation period
data2 = data[(data['mac_n']>1)|(data['cia_n']>1)|(data['comp_n']>1)]
data_s2 = data_s[(data_s['mac_n']>1)|(data_s['cia_n']>1)|(data_s['comp_n']>1)]

data4 = data[data['comp_n']>1]
#filter for ones with gaps >10Myr 
data3=data2[(data2['cia_ngl100']>0)|(data2['comp_ngl100']>0)|((data2['cia_n']-data2['cia_ngl100']-data2['cia_ngl10']-1)>0)]
print(f'There are {len(data3)} dynamos with gaps in generation > 10Myr')

print('There are',len(data2[data2['mac_n']>1]),'intermittent MAC dynamos')
print('There are',len(data2[data2['cia_n']>1]),'intermittent CIA dynamos')
print('There are',len(data2[data2['comp_n']>1]),'intermittent compositional dynamos')

#percent on
data2.loc[:,'comp_perc_on']=data2['comp_dur']/(data2['comp_off']-data2['comp_on'])
data_s2.loc[:,'comp_perc_on']=data_s2['comp_dur']/(data_s2['comp_off']-data_s2['comp_on'])
#difference between datasets
data_s.loc[:,'comp_diff']=data['comp_n']-data_s['comp_n'] #difference in number of on periods for same run
data_s2.loc[:,'comp_perc_diff']=data_s2['comp_perc_on']-data2['comp_perc_on'] #difference in percentage on for same run


######### try numpy histogram instead #####################################
w=1/len(data)
fig, ax = plt.subplots(nrows=2,ncols=1,sharex=True,tight_layout=True)
ax[0].hist((data2[data2['r']==100]['comp_perc_on'],data2[data2['r']==200]['comp_perc_on'],data2[data2['r']==300]['comp_perc_on'],data2[data2['r']==400]['comp_perc_on']),range=(0,1),weights=(np.ones([len(data2[data2['r']==100])])*w,np.ones([len(data2[data2['r']==200])])*w,np.ones([len(data2[data2['r']==300])])*w,np.ones([len(data2[data2['r']==400])])*w),stacked=True)
ax[0].set_ylabel('Fraction of total runs')
ax[0].set_title('Same timestep in solidification')
ax[0].set_ylim([0,0.6])
ax[1].hist((data_s2[data_s2['r']==100]['comp_perc_on'],data_s2[data_s2['r']==200]['comp_perc_on'],data_s2[data_s2['r']==300]['comp_perc_on'],data_s2[data_s2['r']==400]['comp_perc_on']),range=(0,1),weights=(np.ones([len(data_s2[data_s2['r']==100])])*w,np.ones([len(data_s2[data_s2['r']==200])])*w,np.ones([len(data_s2[data_s2['r']==300])])*w,np.ones([len(data_s2[data_s2['r']==400])])*w),stacked=True)
ax[1].set_xlabel('Fraction of time dynamo is on \n between first and last generation times')
ax[1].set_ylabel('Fraction of total runs')
ax[1].set_title('Smaller timestep in solidification')
ax[1].set_ylim([0,0.6])
ax[1].legend(labels=['100km','200km','300km','400km'])
fig.suptitle('Compositional dynamo')
#plt.savefig('../Plots/intermittent_hist.png',bbox_inches='tight')
########## How many have fractional on times < 80% ###########################
frac80=len(data2[data2['comp_perc_on']<0.8])*100/len(data)
frac70=len(data2[data2['comp_perc_on']<0.7])*100/len(data)
print(f"Originally {frac80:.0f} % of compositional dynamos are on less than 80% of the time ")
frac80=len(data_s2[data_s2['comp_perc_on']<0.8])*100/len(data_s)
frac70=len(data_s2[data_s2['comp_perc_on']<0.7])*100/len(data_s)
print(f"After the change {frac80:.0f} % of compositional dynamos are on less than 80% of the time ")
#### how much did changing the timestep improve things
impn = len(data_s2[data_s2['comp_diff']>0]) #reduction in on periods
impperc = len(data_s2[data_s2['comp_perc_diff']>0]) #increase in percentage of time on

print(f'{impn/len(data)*100:.1f}% of runs had a reduction in on periods')
print(f'{impperc*100/len(data):.1f}% of runs had an increase in fractional duration')

