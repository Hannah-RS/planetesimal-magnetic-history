#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Intermittence analysis - investigate which runs have multiple periods of dynamo generation
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

folder = 'Fullrun2'
data = pd.read_csv(f'../Results_combined/{folder}/all_sucess_info.csv',delimiter=',',skiprows=[1],header=0,index_col=False)
data['r']=data['r']/1e3 #rescale to km

#filter for more than 1 dynamo generation period
data2 = data[(data['mac_n']>1)|(data['cia_n']>1)|(data['comp_n']>1)]
data4 = data[data['comp_n']>1]
#filter for ones with gaps >10Myr 
data3=data2[(data2['cia_ngl100']>0)|(data2['comp_ngl100']>0)|((data2['cia_n']-data2['cia_ngl100']-data2['cia_ngl10']-1)>0)]
print(f'There are {len(data3)} dynamos with gaps in generation > 10Myr')

print('There are',len(data2[data2['mac_n']>1]),'intermittent MAC dynamos')
print('There are',len(data2[data2['cia_n']>1]),'intermittent CIA dynamos')
print('There are',len(data2[data2['comp_n']>1]),'intermittent compositional dynamos')

#percent on
data2.loc[:,'comp_perc_on']=data2['comp_dur']/(data2['comp_off']-data2['comp_on'])

data2.loc[:,'cia_perc_on']=data2['cia_dur']/(data2['cia_off']-data2['cia_on'])

######### try numpy histogram instead #####################################
w=1/len(data)
plt.figure()
plt.hist((data2[data2['r']==100]['comp_perc_on'],data2[data2['r']==200]['comp_perc_on'],data2[data2['r']==300]['comp_perc_on'],data2[data2['r']==400]['comp_perc_on']),range=(0,1),weights=(np.ones([len(data2[data2['r']==100])])*w,np.ones([len(data2[data2['r']==200])])*w,np.ones([len(data2[data2['r']==300])])*w,np.ones([len(data2[data2['r']==400])])*w),stacked=True)
plt.legend(labels=['100km','200km','300km','400km'])
plt.xlabel('Fraction of time dynamo is on \n between first and last generation times')
plt.ylabel('Fraction of total runs')
plt.title('Compositional dynamo')
plt.savefig('../Plots/intermittent_hist.png',bbox_inches='tight')
########## How many have fractional on times < 80% ###########################
frac80=len(data2[data2['comp_perc_on']<0.8])*100/len(data)
frac70=len(data2[data2['comp_perc_on']<0.7])*100/len(data)
print(f"{frac80:.0f} % of compositional dynamos are on less than 80% of the time ")
