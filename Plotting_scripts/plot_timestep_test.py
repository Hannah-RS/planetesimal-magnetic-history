#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare outputs of two models by plotting their fluxes and temperatures with time overlaid,
as well as magnetic Reynolds number and solid core fraction
"""
#import modules
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys
# setting path
sys.path.append('../')
from load_info import load_run_info
#choose your runs
run1 = 24
run2 = 26 
run3 = 14
run4 = 28
#choose model labels
model1 = 'No $^{56}$Fe, switch'
model2 = 'No $^{56}$Fe, no switch'
model3 = '$10^{-7}$ $^{56}$Fe/$^{60}$Fe, switch'
model4 = '$10^{-7}$ $^{56}$Fe/$^{60}$Fe, no switch'

#scale time to Myr
from plotting_constants import Myr

path = '../Results_combined/'
#import data from npz file - run1
npzfile = np.load(f'{path}run_{run1}.npz')
Ra1 = npzfile['Ra']
Racrit1 = npzfile['Racrit']
t1 = npzfile['t'] #time in s
Flux1 = npzfile['Flux']
f1 = npzfile['f']
Rem_t1 = npzfile['Rem_t'][:2,:] # magnetic Reynolds number from thermal convection based on MAC balance
Rem_c1 = npzfile['Rem_c'] # magnetic Reynolds number from compositional convection 
B1 = npzfile['B']

Fs1 = Flux1[0]
Fcmb1 = Flux1[1]
Fad1 = Flux1[2]
Frad1 = Flux1[3]

t_plot1 = t1/Myr

#import data from npz file - run2
npzfile = np.load(f'{path}run_{run2}.npz')
Ra2 = npzfile['Ra']
Racrit2 = npzfile['Racrit']
t2 = npzfile['t'] #time in s
Flux2 = npzfile['Flux']
f2 = npzfile['f']
Rem_t2 = npzfile['Rem_t'][:2,:] # magnetic Reynolds number from thermal convection  based on MAC balance
Rem_c2 = npzfile['Rem_c'] # magnetic Reynolds number from compositional convection
B2 = npzfile['B']
 
Fs2 = Flux2[0]
Fcmb2 = Flux2[1]
Fad2 = Flux2[2]
Frad2 = Flux2[3]

t_plot2 = t2/Myr

#import data from npz file - run3
npzfile = np.load(f'{path}run_{run3}.npz')
Ra3 = npzfile['Ra']
Racrit3 = npzfile['Racrit'] 
t3 = npzfile['t'] #time in s
Flux3 = npzfile['Flux']
f3 = npzfile['f']
Rem_t3 = npzfile['Rem_t'][:2,:] # magnetic Reynolds number from thermal convection  based on MAC balance
Rem_c3 = npzfile['Rem_c'] # magnetic Reynolds number from compositional convection
B3 = npzfile['B']

Fs3 = Flux3[0]
Fcmb3 = Flux3[1]
Fad3 = Flux3[2]
Frad3 = Flux3[3]

t_plot3 = t3/Myr

#import data from npz file - run4
npzfile = np.load(f'{path}run_{run4}.npz')
Ra4 = npzfile['Ra']
Racrit4 = npzfile['Racrit']  
t4 = npzfile['t'] #time in s
Flux4 = npzfile['Flux']
f4 = npzfile['f']
Rem_t4 = npzfile['Rem_t'][:2,:] # magnetic Reynolds number from thermal convection  based on MAC balance
Rem_c4 = npzfile['Rem_c'] # magnetic Reynolds number from compositional convection
B4 = npzfile['B']

Fs4 = Flux4[0]
Fcmb4 = Flux4[1]
Fad4 = Flux4[2]
Frad4 = Flux4[3]

t_plot4 = t4/Myr

# #import label info - read in from correct row in csv
r, dr, tstart, tstep, viscosity = load_run_info(1,'../Results_combined/Timestep_test/run_info.csv')
 
threshold = 10

################### Just Ra ####################################################
plt.figure(tight_layout=True)
plt.plot(t_plot1,Ra1,color='mediumblue',label=f'{model1}')
plt.plot(t_plot2,Ra2,color='cornflowerblue',label=f'{model2}',linestyle='dashed')
plt.plot(t_plot3,Ra3,color='forestgreen',label=f'{model3}')
plt.plot(t_plot4,Ra4,color='limegreen',label=f'{model4}',linestyle='dashed')
plt.plot(t_plot1,Racrit1,color='navy',label='critical Ra',linestyle='dotted')
#plt.plot(t_plot2,Racrit2,color='black',label='critical Ra',linestyle='dotted')
#plt.plot(t_plot3,Racrit3,color='black',label='critical Ra',linestyle='dotted')
plt.plot(t_plot4,Racrit4,color='darkgreen',label=f'critical Ra -{model4}',linestyle='dotted')
plt.yscale('log')
plt.xscale('log')
plt.xlim(right=500)
plt.ylim([1e6,1e11])
plt.ylabel('Ra')
plt.xlabel('t/Myr')
plt.legend(ncols=2)
#plt.savefig('../Plots/Xs_r_tests/Ra_switch.png',dpi=600)

################# Rayleigh numbers and surface fluxes ##########################################

fig, ax1 = plt.subplots(tight_layout=True)
plt.plot(t_plot1,Ra1,color='navy',label=f'Ra {model1}')
plt.plot(t_plot1,Racrit1,color='lightcoral',label='critical Ra')
plt.plot(t_plot3,Ra3,color='mediumblue',label=f'{model3}',linestyle='dotted')
plt.plot(t_plot3,Racrit3,color='maroon',linestyle='dotted')
plt.plot(t_plot4,Ra4,color='dodgerblue',label=f'{model4}',linestyle='-.')
plt.plot(t_plot4,Racrit4,color='firebrick',linestyle='-.')
plt.plot(t_plot2,Ra2,color='cornflowerblue',label=f'{model2}',linestyle='dashed')
plt.plot(t_plot2,Racrit2,color='indianred',linestyle='dashed')
plt.yscale('log')
plt.xlim([1,12])
plt.ylim([1e6,1e13])
plt.ylabel('Ra')
plt.xlabel('t/Myr')
plt.legend()
ax2 = ax1.twinx()
plt.plot(t_plot1,Fs1,color='darkgreen',label=f'F$_s$ {model1}')
plt.plot(t_plot3,Fs3,color='forestgreen',linestyle='dotted',label=f'F$_s$ {model3}')
plt.plot(t_plot4,Fs4,color='limegreen',linestyle='-.',label=f'F$_s$ {model4}')
plt.plot(t_plot2,Fs2,color='lightgreen',linestyle='dashed',label=f'F$_s$ {model2}')
plt.ylabel('F$_{s}$/W$m^{-2}$')
plt.ylim([0.1,10])
plt.legend(loc='lower left')
#plt.savefig('../Plots/Timestep_test/Ra_Fs.png',dpi=600)

################### Flux plots ################################################

plt.figure(tight_layout=True)
plt.suptitle('Thermal evolution of a {:.0f}km asteroid '.format(r/1e3))#\n Tm0 = {}K, Tsolidus ={}K 

xmin=tstart
#fluxes as function of time

plt.loglog(t_plot1,Fs1,label=f'$F_s$, {model1}',color='royalblue')
plt.loglog(t_plot1,Fcmb1,label='$F_{CMB}$',color='firebrick')
plt.loglog(t_plot1,Fad1,label='$F_{ad}$',color='forestgreen')
plt.loglog(t_plot1,Frad1,label='$F_{rad}$',color='orange')

plt.loglog(t_plot2,Fs2,color='royalblue',linestyle='--',label=f'{model2}')
plt.loglog(t_plot2,Fcmb2,color='firebrick',linestyle='--')
plt.loglog(t_plot2,Fad2,color='forestgreen',linestyle='--')
plt.loglog(t_plot2,Frad2,color='orange',linestyle='--')

plt.loglog(t_plot3,Fs3,color='royalblue',linestyle='dotted',label=f'{model3}')
plt.loglog(t_plot3,Fcmb3,color='firebrick',linestyle='dotted')
plt.loglog(t_plot3,Fad3,color='forestgreen',linestyle='dotted')
plt.loglog(t_plot3,Frad3,color='orange',linestyle='dotted')

plt.loglog(t_plot4,Fs4,color='royalblue',linestyle='-.',label=f'{model4}')
plt.loglog(t_plot4,Fcmb4,color='firebrick',linestyle='-.')
plt.loglog(t_plot4,Fad4,color='forestgreen',linestyle='-.')
plt.loglog(t_plot4,Frad4,color='orange',linestyle='-.')

plt.xlabel('Time/ Myr')
plt.ylim([1e-3,1e2])   #use these limits when comparing runs
plt.ylabel('Flux/ W$m^{-2}$')
plt.legend(loc='upper right',ncol=2,fontsize='small')
#plt.savefig('../Plots/Tflux_comp.png',dpi=450)



############# Compositional and thermal on same plot ##########################
with sns.plotting_context('talk',font_scale=0.8):
    plt.figure(tight_layout=True,figsize=[10,3.5])
    plt.title('Magnetic field generation in a {:.0f}km asteroid '.format(r/1e3))
    plt.plot(t_plot1,Rem_t1[0,:],label=model1,color='cornflowerblue',alpha=0.7)
    plt.plot(t_plot2,Rem_t2[0,:],label=model2,color='mediumblue',linestyle='dashed')
    plt.plot(t_plot3,Rem_t3[0,:],label=model3,color='navy',linestyle='dotted')
    plt.plot(t_plot4,Rem_t4[0,:],label=model4,color='royalblue',linestyle='-.')
    plt.plot(t_plot1,Rem_t1[1,:],label='CIA',color='lightcoral',alpha=0.7)
    plt.plot(t_plot2,Rem_t2[1,:],color='indianred',linestyle='dashed')
    plt.plot(t_plot3,Rem_t3[1,:],color='maroon',linestyle='dotted')
    plt.plot(t_plot4,Rem_t4[1,:],color='firebrick',linestyle='-.')
    plt.plot(t_plot1[Rem_c1>10],Rem_c1[Rem_c1>10],color='lightgreen',alpha=0.7,label='compositional')
    plt.plot(t_plot2[Rem_c2>10],Rem_c2[Rem_c2>10],color='limegreen',linestyle='dashed')
    plt.plot(t_plot3[Rem_c3>10],Rem_c3[Rem_c3>10],color='darkgreen',linestyle='dotted')
    plt.plot(t_plot4[Rem_c4>10],Rem_c4[Rem_c4>10],color='turquoise',linestyle='-.')
    plt.yscale('log')
    #plt.xscale('log')
    plt.hlines(10,xmin=0,xmax=t_plot2[-1],color='k',linestyle='--')
    plt.xlabel('Time/Myr')
    plt.ylabel('Rem')
    plt.legend(loc='upper right',ncol=2)
    plt.ylim([7,30])
    plt.xlim([0,20])
    #plt.savefig('../Plots/Timestep_test/cia_start2.png',dpi=600)

########### MAC Rem and CMB flux on same plot #############################
fig, ax1 = plt.subplots(tight_layout=True,figsize=[10,3.5])
plt.title('Magnetic field generation in a {:.0f}km asteroid '.format(r/1e3))
plt.plot(t_plot1,Rem_t1[0,:],label=model1,color='cornflowerblue',alpha=0.7)
plt.plot(t_plot2,Rem_t2[0,:],label=model2,color='mediumblue',linestyle='dashed')
plt.plot(t_plot3,Rem_t3[0,:],label=model3,color='navy',linestyle='dotted')
plt.plot(t_plot4,Rem_t4[0,:],label=model4,color='royalblue',linestyle='-.')
plt.xlabel('Time/Myr')
plt.ylabel('Rem - MAC')
plt.legend(loc='upper right',ncol=2)
plt.ylim([7,20])
plt.xlim([0,20])
ax2 = ax1.twinx()
plt.plot(t_plot1,Fcmb1,color='lightgreen')
plt.plot(t_plot2,Fcmb2,color='limegreen',linestyle='dashed')
plt.plot(t_plot3,Fcmb3,color='darkgreen',linestyle='dotted')
plt.plot(t_plot4,Fcmb4,color='turquoise',linestyle='-.')
plt.ylabel('F$_{cmb}$/W$m^{-2}$')
plt.ylim([0,0.4])
#ax1yscale('log')
#plt.xscale('log')
#plt.hlines(10,xmin=0,xmax=t_plot2[-1],color='k',linestyle='--')



############### Field strength but thermal only ############################

with sns.plotting_context('talk',font_scale=0.8):
    plt.figure(tight_layout=True,figsize=[10,3.5])
    plt.title('Field strengths for supercritical Re$_m$')
    plt.plot(t_plot1[Rem_t1[0,:]>threshold],B1[0,Rem_t1[0,:]>threshold]/1e-6,label=f'{model1} - ml',color='lightgreen',alpha=0.7)
    plt.plot(t_plot1[Rem_t1[0,:]>threshold],B1[1,Rem_t1[0,:]>threshold]/1e-6,label=f'{model1} - MAC',color='cornflowerblue',alpha=0.7)
    plt.plot(t_plot1[Rem_t1[1,:]>threshold],B1[2,Rem_t1[1,:]>threshold]/1e-6,label=f'{model1} - CIA',color='lightcoral',alpha=0.7)
    plt.plot(t_plot2[Rem_t2[0,:]>threshold],B2[0,Rem_t2[0,:]>threshold]/1e-6,label=f'{model2} - ml',color='limegreen',linestyle='dashed',alpha=0.7)
    plt.plot(t_plot2[Rem_t2[0,:]>threshold],B2[1,Rem_t2[0,:]>threshold]/1e-6,color='mediumblue',linestyle='dashed')
    plt.plot(t_plot2[Rem_t2[1,:]>threshold],B2[2,Rem_t2[1,:]>threshold]/1e-6,color='indianred',linestyle='dashed')
    plt.plot(t_plot3[Rem_t3[0,:]>threshold],B3[0,Rem_t3[0,:]>threshold]/1e-6,label=f'{model3} - ml',color='darkgreen',linestyle='dotted')
    plt.plot(t_plot3[Rem_t3[0,:]>threshold],B3[1,Rem_t3[0,:]>threshold]/1e-6,color='navy',linestyle='dotted')
    plt.plot(t_plot3[Rem_t3[1,:]>threshold],B3[2,Rem_t3[1,:]>threshold]/1e-6,color='maroon',linestyle='dotted')
    plt.plot(t_plot4[Rem_t4[0,:]>threshold],B4[0,Rem_t4[0,:]>threshold]/1e-6,label=f'{model4} - ml',color='turquoise',linestyle='-.')
    plt.plot(t_plot4[Rem_t4[0,:]>threshold],B4[1,Rem_t4[0,:]>threshold]/1e-6,color='royalblue',linestyle='-.')
    plt.plot(t_plot4[Rem_t4[1,:]>threshold],B4[2,Rem_t4[1,:]>threshold]/1e-6,color='firebrick',linestyle='-.')
    #plt.yscale('log')
    plt.xlim([3,12])
    plt.xlabel('Time/Myr')
    plt.ylabel('B/$\mu$T')
    plt.legend(loc='upper left')
    #plt.savefig('../Plots/Timestep_test/maxB.png',dpi=600)


