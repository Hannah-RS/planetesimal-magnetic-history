#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare outputs of two models by plotting their fluxes and temperatures with time overlaid,
as well as magnetic Reynolds number and solid core fraction
"""
#import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from load_info import load_run_info
#choose your runs
run1 = 2
run2 = 3 
run3 = 4
#choose model labels
model1 = 'Bryson et. al. (2019)'
model2 = 'Bryson $\phi_{RCMF}$=0.3'
model3 = 'Bryson $\eta_0$=10$^{14}$Pas'
# model1= 'Dodds et. al. (2021)'
# model2 = 'Bryson et al. (2019)'
# model3='Sterenborg & Crowley (2013)'
#scale time to Myr
from parameters import Myr

#import data from npz file - run1
npzfile = np.load('Results_combined/run_{}.npz'.format(run1))
Tm1 = npzfile['Tm_mid'] 
Tm_surf1 = npzfile['Tm_surf'] 
Tc1= npzfile['Tc'] 
t1 = npzfile['t'] #time in s
Flux1 = npzfile['Flux']
f1 = npzfile['f']
Rem_t1 = npzfile['Rem_t'][0,:] # magnetic Reynolds number from thermal convection based on MAC balance
Rem_c1 = npzfile['Rem_c'] # magnetic Reynolds number from compositional convection 
#Rem1 = Rem_t1
#Rem1[Rem_c1>Rem_t1]=Rem_c1[Rem_c1>Rem_t1]

Fs1 = Flux1[0]
Fcmb1 = Flux1[1]
Fad1 = Flux1[2]
Frad1 = Flux1[3]

t_plot1 = t1/Myr

#import data from npz file - run2
npzfile = np.load('Results_combined/run_{}.npz'.format(run2))
Tm2 = npzfile['Tm_mid'] 
Tm_surf2 = npzfile['Tm_surf'] 
Tc2= npzfile['Tc'] 
t2 = npzfile['t'] #time in s
Flux2 = npzfile['Flux']
f2 = npzfile['f']
Rem_t2 = npzfile['Rem_t'][0,:] # magnetic Reynolds number from thermal convection 
Rem_c2 = npzfile['Rem_c'] # magnetic Reynolds number from compositional convection
#Rem2 = Rem_t2
#Rem2[Rem_c2>Rem_t2]=Rem_c2[Rem_c2>Rem_t2]

Fs2 = Flux2[0]
Fcmb2 = Flux2[1]
Fad2 = Flux2[2]
Frad2 = Flux2[3]

t_plot2 = t2/Myr

#import data from npz file - run3
npzfile = np.load('Results_combined/run_{}.npz'.format(run3))
Tm3 = npzfile['Tm_mid'] 
Tm_surf3 = npzfile['Tm_surf'] 
Tc3= npzfile['Tc'] 
t3 = npzfile['t'] #time in s
Flux3 = npzfile['Flux']
f3 = npzfile['f']
Rem_t3 = npzfile['Rem_t'][0,:] # magnetic Reynolds number from thermal convection 
Rem_c3 = npzfile['Rem_c'] # magnetic Reynolds number from compositional convection
#Rem3 = Rem_t3
#Rem3[Rem_c3>Rem_t3]=Rem_c3[Rem_c3>Rem_t3]

Fs3 = Flux3[0]
Fcmb3 = Flux3[1]
Fad3 = Flux3[2]
Frad3 = Flux3[3]

t_plot3 = t3/Myr

# #import label info - read in from correct row in csv
r, dr, tstart, tend, tstep, tsolid, cond_t, dt, viscosity = load_run_info(2,'run_info.csv',True)

################### Main Plot #########################################

# make  log-log plot - similar to Bryson 2019

plt.figure(tight_layout=True)
plt.suptitle('Thermal evolution of a {:.0f}km asteroid '.format(r/1e3))#\n Tm0 = {}K, Tsolidus ={}K 

xmin=tstart

#temperatures as function of time
plt.subplot(2,1,1)
plt.semilogx(t_plot1,Tm1,label='mantle',color='red')
plt.semilogx(t_plot1,Tc1,label='core',color='black')
plt.semilogx(t_plot2,Tm2,color='red',linestyle='--',label=model2)
plt.semilogx(t_plot2,Tc2,color='black',linestyle='--')
plt.vlines(cond_t,ymin=min(Tm1),ymax=1600,color='black',linestyle='dotted',label='conduction')
#plt.xlim([xmin,max(t_plot)])
#plt.xlabel('Time/ Myr')
#plt.xlim([xmin,500])  #use these limits when comparing runs
#plt.ylim([1400,1650]) #use these limits when comparing runs
plt.ylabel('T/K')
plt.legend(loc='upper right',ncol=2,fontsize='small')

#fluxes as function of time
plt.subplot(2,1,2)
plt.loglog(t_plot1,Fs1,label='$F_s$',color='royalblue')
plt.loglog(t_plot1,Fcmb1,label='$F_{CMB}$',color='firebrick')
plt.loglog(t_plot1,Fad1,label='$F_{ad}$',color='forestgreen')
plt.loglog(t_plot1,Frad1,label='$F_{rad}$',color='orange')

plt.loglog(t_plot2,Fs2,color='royalblue',linestyle='--')
plt.loglog(t_plot2,Fcmb2,color='firebrick',linestyle='--')
plt.loglog(t_plot2,Fad2,color='forestgreen',linestyle='--')
plt.loglog(t_plot2,Frad2,color='orange',linestyle='--')
plt.vlines(cond_t,ymin=min(Fad1),ymax=max(Frad1),color='black',linestyle='dotted',label='conduction')
plt.xlabel('Time/ Myr')
plt.ylim([1e-3,1e2])   #use these limits when comparing runs
# plt.xlim([xmin,500])   #use these limits when comparing runs
plt.ylabel('Flux/ W$m^{-2}$')
plt.legend(loc='upper right',ncol=2,fontsize='small')
#plt.savefig('Plots/Tflux_comp.png',dpi=450)

### Do same thing for Rem and f 
# with sns.plotting_context('talk',font_scale=0.8):
#     plt.figure(tight_layout=True,figsize=[10,7])
#     plt.suptitle('Thermal evolution of a {:.0f}km asteroid '.format(r/1e3))
#     plt.subplot(2,1,1)
#     plt.loglog(t_plot1,Rem1,label=model1,color='cornflowerblue',alpha=0.7)
#     plt.loglog(t_plot2,Rem2,label=model2,color='mediumblue',linestyle='dashed')
#     plt.loglog(t_plot3,Rem3,label=model3,color='navy',linestyle='dotted')
#     #plt.xlim([xmin,max(t_plot)])
#     plt.hlines(10,xmin=0,xmax=t_plot2[-1],color='k',linestyle='--')
#     #plt.xlabel('Time/Myr')
#     plt.ylabel('Rem')
#     plt.legend(loc='upper left',ncol=2)
#     plt.ylim([1,100])
    
#     plt.subplot(2,1,2)
#     plt.semilogx(t_plot1,f1,label=model1,color='cornflowerblue',alpha=0.7)
#     plt.semilogx(t_plot2,f2,label=model2,color='mediumblue',linestyle='dashed')
#     plt.semilogx(t_plot3,f3,label=model3,color='navy',linestyle='dotted')
#     #plt.xlim([xmin,max(t_plot)])
#     plt.xlabel('Time/ Myr')
#     plt.ylabel('f')
#     plt.legend()
#     #plt.savefig('Plots/Remf_comp_bbb.png',dpi=600)

############# Just Rem ###################################################
tsolid3 = t_plot3[f3<0.985]
Rem_solid3 = Rem_c3[f3<0.985]
with sns.plotting_context('talk',font_scale=0.8):
    plt.figure(tight_layout=True,figsize=[10,7])
    plt.title('Magnetic field generation in a {:.0f}km asteroid '.format(r/1e3))
    plt.subplot(2,1,1)
    plt.title('Thermal dynamo')
    plt.loglog(t_plot1,Rem_t1,label=model1,color='cornflowerblue',alpha=0.7)
    plt.loglog(t_plot2,Rem_t2,label=model2,color='mediumblue',linestyle='dashed')
    plt.loglog(t_plot3,Rem_t3,label=model3,color='navy',linestyle='dotted')
    plt.yscale('log')
    plt.xscale('log')
    plt.hlines(10,xmin=0,xmax=t_plot2[-1],color='k',linestyle='--')
    plt.xlabel('Time/Myr')
    plt.ylabel('Rem')
    plt.legend(loc='upper right',ncol=2)
    plt.ylim([7,100])
    plt.xlim([xmin,max(t_plot1)])
    #plt.savefig('Plots/Rem_transfer.png',dpi=600)
    
    plt.subplot(2,1,2)
    plt.title('Compositional dynamo')
    plt.loglog(t_plot1[Rem_c1>10],Rem_c1[Rem_c1>10],label=model1,color='cornflowerblue',alpha=0.7)
    plt.loglog(t_plot2[Rem_c2>10],Rem_c2[Rem_c2>10],label=model2,color='mediumblue',linestyle='dashed')
    plt.loglog(tsolid3[Rem_solid3>10],Rem_solid3[Rem_solid3>10],label=model3,color='navy',linestyle='dotted')
    # plt.scatter(t_plot1[Rem_c1>10],Rem_c1[Rem_c1>10],label=model1,color='cornflowerblue',alpha=0.7,s=3)
    # plt.scatter(t_plot2[Rem_c2>10],Rem_c2[Rem_c2>10],label=model2,color='mediumblue',linestyle='dashed',s=3)
    # plt.scatter(tsolid3[Rem_solid3>10],Rem_solid3[Rem_solid3>10],label=model3,color='navy',linestyle='dotted',s=3)
    #plt.xlim([xmin,max(t_plot)])
    plt.yscale('log')
    plt.xscale('log')
    plt.hlines(10,xmin=0,xmax=t_plot2[-1],color='k',linestyle='--')
    plt.xlabel('Time/Myr')
    plt.ylabel('Rem')
    plt.legend(loc='upper left',ncol=2)
    plt.ylim([7,100])
    plt.xlim([xmin,max(t_plot1)])
    #plt.savefig('Plots/Rem_ukpf.png',dpi=600)

with sns.plotting_context('talk',font_scale=0.8):
    plt.figure(tight_layout=True,figsize=[10,3.5])
    plt.title('Magnetic field generation in a {:.0f}km asteroid '.format(r/1e3))
    plt.loglog(t_plot1,Rem_t1,label=model1,color='cornflowerblue',alpha=0.7)
    plt.loglog(t_plot2,Rem_t2,label=model2,color='mediumblue',linestyle='dashed')
    plt.loglog(t_plot3,Rem_t3,label=model3,color='navy',linestyle='dotted')
    plt.loglog(t_plot1[Rem_c1>10],Rem_c1[Rem_c1>10],color='cornflowerblue',alpha=0.7)
    plt.loglog(t_plot2[Rem_c2>10],Rem_c2[Rem_c2>10],color='mediumblue',linestyle='dashed')
    plt.loglog(tsolid3[Rem_solid3>10],Rem_solid3[Rem_solid3>10],color='navy',linestyle='dotted')
    plt.yscale('log')
    plt.xscale('log')
    plt.hlines(10,xmin=0,xmax=t_plot2[-1],color='k',linestyle='--')
    plt.xlabel('Time/Myr')
    plt.ylabel('Rem')
    plt.legend(loc='upper left',ncol=2)
    plt.ylim([7,100])
    plt.xlim([xmin,max(t_plot1)])
    #plt.savefig('Plots/Rem_ukpf2.png',dpi=600)
tsolid3[Rem_solid3>10][0]   
##################### Just f ################################################
with sns.plotting_context('talk',font_scale=0.8):
    plt.figure(tight_layout=True,figsize=[10,3.5])
    plt.title('Thermal evolution of a {:.0f}km asteroid '.format(r/1e3))
    plt.semilogx(t_plot1,f1,label=model1,color='cornflowerblue',alpha=0.7)
    plt.semilogx(t_plot2,f2,label=model2,color='mediumblue',linestyle='dashed')
    plt.semilogx(t_plot3,f3,label=model3,color='navy',linestyle='dotted')
    #plt.xlim([xmin,max(t_plot)])
    plt.xlabel('Time/ Myr')
    plt.ylabel('f')
    plt.legend()
    #plt.savefig('Plots/f_transfer.png',dpi=600)


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
    #plt.savefig('Plots/viscosity_bryson_dodds_talk.png',dpi=600)

##################### Flux plot ######################################    
plt.figure()
plt.subplot(2,1,1)
plt.loglog(t_plot3,Fcmb3,color='cornflowerblue',linestyle='dotted',label=model3)
plt.loglog(t_plot2,Fcmb2,color='blue',linestyle='--',label=model2)
plt.loglog(t_plot1,Fcmb1,label=model1,color='black')
plt.xlabel('Time/ Myr')
plt.ylim([1e-4,1e2])   #use these limits when comparing runs
plt.ylabel('CMB Flux/ W$m^{-2}$')
plt.legend(loc='upper right')

plt.subplot(2,1,2)
plt.semilogx(t_plot2,f2,label=model2,color='blue',linestyle='--')
plt.semilogx(t_plot1,f1,label=model1,color='black')
plt.xlabel('Time/ Myr')
plt.ylabel('Fractional of \n core solidified')
plt.legend()
#plt.savefig('Plots/viscosity_flux_comp.pdf',dpi=600)