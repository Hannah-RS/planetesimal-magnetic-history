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
run1 = 1
run2 = 5 
run3 = 18
#choose model labels
model1 = 'dt,dr'
model2 = '0.5dt,dr'
model3 = '0.8dt,dr'

#scale time to Myr
from plotting_constants import Myr

path = '../Results_combined/Timestep_test/'
#import data from npz file - run1
npzfile = np.load(f'{path}run_{run1}.npz')
Tm1 = npzfile['Tm_mid']  
Tc1= npzfile['Tc'] 
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
Tm2 = npzfile['Tm_mid'] 
Tc2= npzfile['Tc'] 
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
Tm3 = npzfile['Tm_mid'] 
Tc3= npzfile['Tc'] 
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

# #import label info - read in from correct row in csv
r, dr, tstart, tstep, viscosity = load_run_info(1,'../Results_combined/Timestep_test/run_info.csv')

################# Core and mantle temp plot for transfer ##########################################
with sns.plotting_context('talk'):
    plt.figure(tight_layout=True,figsize=[10,7])
    plt.title(f'Thermal evolution of a {r/1e3:.0f}km asteroid ')
    xmin=tstart
    #temperatures as function of time
    plt.semilogx(t_plot1,Tm1,label='mantle',color='red')  
    plt.semilogx(t_plot1,Tc1,label='core',color='black')
    plt.semilogx(t_plot2,Tm2,color='red',linestyle='--')
    plt.semilogx(t_plot2,Tc2,color='black',linestyle='--',label=model2)
    plt.semilogx(t_plot3,Tm3,color='red',linestyle='dotted',label=model3)
    plt.semilogx(t_plot3,Tc3,color='black',linestyle='dotted',label=model1)
    plt.semilogx(t_plot1,Tc1,label=model1,color='black')
    plt.ylabel('T/K')
    plt.xlabel('Time/Myr')
    plt.legend(loc='upper right',ncol=2,fontsize='small')
    #plt.savefig('../Plots/viscosity_bryson_dodds_talk.png',dpi=600)

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

plt.xlabel('Time/ Myr')
plt.ylim([1e-3,1e2])   #use these limits when comparing runs
plt.ylabel('Flux/ W$m^{-2}$')
plt.legend(loc='upper right',ncol=2,fontsize='small')
#plt.savefig('../Plots/Tflux_comp.png',dpi=450)


############# Just Rem ###################################################
with sns.plotting_context('talk',font_scale=0.8):
    plt.figure(tight_layout=True,figsize=[10,7])
    plt.suptitle('Magnetic field generation in a {:.0f}km asteroid '.format(r/1e3))
    plt.subplot(2,1,1)
    plt.title('Thermal dynamo')
    plt.loglog(t_plot1,Rem_t1[0,:],label=model1,color='cornflowerblue',alpha=0.7)
    plt.loglog(t_plot2,Rem_t2[0,:],label=model2,color='mediumblue',linestyle='dashed')
    plt.loglog(t_plot3,Rem_t3[0,:],label=model3,color='navy',linestyle='dotted')
    plt.yscale('log')
    plt.xscale('log')
    plt.hlines(10,xmin=0,xmax=t_plot2[-1],color='k',linestyle='--')
    plt.xlabel('Time/Myr')
    plt.ylabel('Rem')
    plt.legend(loc='upper right',ncol=2)
    plt.ylim([7,100])
    plt.xlim([xmin,max(t_plot1)])
    
    plt.subplot(2,1,2)
    plt.title('Compositional dynamo')
    plt.loglog(t_plot1[Rem_c1>10],Rem_c1[Rem_c1>10],label=model1,color='cornflowerblue',alpha=0.7)
    plt.loglog(t_plot2[Rem_c2>10],Rem_c2[Rem_c2>10],label=model2,color='mediumblue',linestyle='dashed')
    plt.loglog(t_plot3[Rem_c3>10],Rem_c3[Rem_c3>10],label=model3,color='navy',linestyle='dotted')
    plt.yscale('log')
    plt.xscale('log')
    plt.hlines(10,xmin=0,xmax=t_plot2[-1],color='k',linestyle='--')
    plt.xlabel('Time/Myr')
    plt.ylabel('Rem')
    plt.legend(loc='upper left',ncol=2)
    plt.ylim([7,100])
    plt.xlim([xmin,max(t_plot1)])
    #plt.savefig('../Plots/Rem_ukpf.png',dpi=600)

############# Compositional and thermal on same plot ##########################
with sns.plotting_context('talk',font_scale=0.8):
    plt.figure(tight_layout=True,figsize=[10,3.5])
    plt.title('Magnetic field generation in a {:.0f}km asteroid '.format(r/1e3))
    plt.loglog(t_plot1,Rem_t1[0,:],label=model1,color='cornflowerblue',alpha=0.7)
    plt.loglog(t_plot2,Rem_t2[0,:],label=model2,color='mediumblue',linestyle='dashed')
    plt.loglog(t_plot3,Rem_t3[0,:],label=model3,color='navy',linestyle='dotted')
    plt.loglog(t_plot1,Rem_t1[1,:],label=model1,color='cornflowerblue',alpha=0.7)
    plt.loglog(t_plot2,Rem_t2[1,:],label=model2,color='mediumblue',linestyle='dashed')
    plt.loglog(t_plot3,Rem_t3[1,:],label=model3,color='navy',linestyle='dotted')
    plt.loglog(t_plot1[Rem_c1>10],Rem_c1[Rem_c1>10],color='cornflowerblue',alpha=0.7)
    plt.loglog(t_plot2[Rem_c2>10],Rem_c2[Rem_c2>10],color='mediumblue',linestyle='dashed')
    plt.loglog(t_plot3[Rem_c3>10],Rem_c3[Rem_c3>10],color='navy',linestyle='dotted')
    #plt.yscale('log')
    plt.xscale('log')
    plt.hlines(10,xmin=0,xmax=t_plot2[-1],color='k',linestyle='--')
    plt.xlabel('Time/Myr')
    plt.ylabel('Rem')
    plt.legend(loc='upper left',ncol=2)
    plt.ylim([10,100])
    plt.xlim([xmin,max(t_plot1)])
    #plt.savefig('../Plots/Rem_ukpf2.png',dpi=600)
    
############# Magnetic field strength ########################################
threshold = 10 #critical Rem

with sns.plotting_context('talk',font_scale=0.8):
    plt.figure(tight_layout=True,figsize=[10,3.5])
    plt.title('Field strengths for supercritical Re$_m$')
    plt.loglog(t_plot1[Rem_t1>threshold],B1[0,Rem_t1>threshold]/1e-6,label=f'{model1} - thermal',color='cornflowerblue',alpha=0.7)
    plt.loglog(t_plot1[Rem_c1>threshold],B1[3,Rem_c1>threshold]/1e-6,label=f'{model1} - compositional',color='mediumblue',alpha=0.7)
    plt.loglog(t_plot2[Rem_t2>threshold],B2[0,Rem_t2>threshold]/1e-6,label=f'{model2} - thermal',color='seagreen',alpha=0.7,linestyle='dashed')
    plt.loglog(t_plot2[Rem_c2>threshold],B2[3,Rem_c2>threshold]/1e-6,color='forestgreen',alpha=0.7,linestyle='dashed')
    plt.loglog(t_plot3[Rem_t3>threshold],B3[0,Rem_t3>threshold]/1e-6,label=f'{model3} - thermal',color='mediumpurple',alpha=0.7,linestyle='dotted')
    plt.loglog(t_plot3[Rem_c3>threshold],B3[3,Rem_c3>threshold]/1e-6,color='darkorchid',alpha=0.7,linestyle='dotted')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Time/Myr')
    plt.ylabel('B/$\mu$T')
    plt.legend(loc='upper right')
    plt.xlim([xmin,max(t_plot1)])

  
##################### Just f ################################################
with sns.plotting_context('talk',font_scale=0.8):
    plt.figure(tight_layout=True,figsize=[10,3.5])
    plt.semilogx(t_plot1,f1,label=model1,color='cornflowerblue',alpha=0.7)
    plt.semilogx(t_plot2,f2,label=model2,color='mediumblue',linestyle='dashed')
    plt.semilogx(t_plot3,f3,label=model3,color='navy',linestyle='dotted')
    #plt.xlim([xmin,max(t_plot)])
    plt.xlabel('Time/ Myr')
    plt.ylabel('f')
    plt.legend()
    #plt.savefig('../Plots/f_transfer.png',dpi=600)

##################### Flux plot ######################################    
plt.figure()
plt.loglog(t_plot3,Fcmb3,color='cornflowerblue',linestyle='dotted',label=model3)
plt.loglog(t_plot2,Fcmb2,color='blue',linestyle='--',label=model2)
plt.loglog(t_plot1,Fcmb1,label=model1,color='black')
plt.loglog(t_plot1,Fad1,label='$F_{ad}$',color='forestgreen',linestyle='dashdot')
plt.xlabel('Time/ Myr')
plt.ylim([1e-4,1e2])   #use these limits when comparing runs
plt.ylabel('CMB Flux/ W$m^{-2}$')
plt.legend(loc='upper right')
