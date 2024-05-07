#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for solving the thermal evolution of an asteroid. 
"""
# import modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time #use this to time the integration

#import time constants and initial conditions
from parameters import  run, t_acc_m, t_end_m, dr, automated, Myr, Ts, f0, r, rc, rcr,\
    kappa_c, save_interval_d, save_interval_t, save_interval_mag, km, Vm, As, rhom, step_m, Xs_0, \
    default, rcmf, Fe0, full_save, B_save, eta0, etal, w, alpha_n, beta, n_cells, icfrac
from viscosity_def import viscosity

if automated == True: 
    import sys
    folder = sys.argv[1]
else:
    folder = 'Results_combined/' #folder where you want to save the results
    ind = None #no index for csv
#set flag for run started
if automated == True:
    from parameters import ind
    auto = pd.read_csv(f'{folder}auto_params.csv')
    auto.loc[ind+1,'status']=0 #indicates started
    auto.to_csv(f'{folder}auto_params.csv',index=False)
else: #save run parameters in run_info file
    run_info = {"run":[run],"r":[r],"rcr":[rcr],"default":[default],"rcmf":[rcmf],"eta0":[eta0], 
                "beta":[beta],"w":[w],"etal":[etal],"alpha_n":[alpha_n],"Xs_0":[Xs_0], 
                "Fe0":[Fe0], "t_acc_m":[t_acc_m], "t_end_m":[t_end_m], "dr":[dr],
                "step_m":[step_m],"icfrac":[icfrac]}
    run_info = pd.DataFrame(run_info)
    run_info.to_csv(f'{folder}run_info.csv',index=False,mode='a',header=False)

#calculate the stencil for the conductive profile, save so can be reloaded in later steps
from stencil import cond_stencil_core, cond_stencil_mantle
from diff_function import differentiation
from full_euler import thermal_evolution

import scipy.sparse as sp

dT_mat_c = cond_stencil_core(r,rc,dr,kappa_c) 
dT_mat_m = cond_stencil_mantle(r,rc,dr,km/rhom)  
sparse_mat_m = sp.dia_matrix(dT_mat_m)
sparse_mat_c = sp.dia_matrix(dT_mat_c)


# define the run number, start and end times
t_acc = t_acc_m*Myr #Accretion time
t_end=t_end_m*Myr #end time 
#step_m=0.1*t_cond_core  #max timestep must be smaller than conductive timestep
n_save_d = int(save_interval_d/step_m)
n_save_t = int(save_interval_t/step_m)

# set initial temperature profile
Tint = np.ones([n_cells])*Ts#first element in the array is at r=0, accrete cold at surface temp 
Tint[-1]=Ts

#Check viscosity profile is monotonically decreasing before start
Ttest = np.linspace(1200,1900,200)
#calculate viscosity
eta_test = viscosity(Ttest)
eta_diff = np.diff(eta_test) #calculate differences with sucessive elements
if np.all(eta_diff<=0):
    print('Viscosity profile is monotonically decreasing - proceeding')
else: #put marker in csv
    if automated == True: #no need to reimport ind as will have been imported earlier
        auto = pd.read_csv(f'{folder}auto_params.csv')
        auto.loc[ind+1,'status']=-1 #indicates viscosity error
        auto.to_csv(f'{folder}auto_params.csv',index=False)
    raise ValueError('Invalid viscosity model')
    
print(f'Beginning run {run}')
print('Initial conditions set')
#%%
########################### Differentiation ###################################
tic = time.perf_counter()
Tdiff, Xfe, Xsi, cp, Ra, Ra_crit, convect, d0, t_diff, H  = differentiation(Tint,
                                                                            t_acc,r, 
                                                                            dr, step_m)
toc = time.perf_counter()
int_time1 = toc - tic  

# update user on progress and plot differentiated temperature profile 
if automated == False:
    rplot= np.arange(0,r+dr,dr)/1e3
    
    plt.figure()
    plt.scatter(rplot,Tdiff[:,-1])
    plt.xlabel('r/km')
    plt.ylabel('Temperature/K')
    plt.title('Temperature profile post differentiation')

#rescale data and save here in case thermal evolution crashes
#relabel so don't change input to model on next step
Tdiffs = Tdiff[:,0::n_save_d]
Xfes = Xfe[:,0::n_save_d]
Xsis = Xsi[:,0::n_save_d]
cps = cp[:,0::n_save_d]
Ras = Ra[0::n_save_d]
Ra_crits = Ra_crit[0::n_save_d]
convects = convect[0::n_save_d]
d0s = d0[0::n_save_d]
t_diffs = t_diff[0::n_save_d]
Hs = H[0::n_save_d]

if full_save == True:
    np.savez_compressed(f'{folder}run_{run}_diff', Tdiff = Tdiffs, Xfe = Xfes, 
                        Xsi = Xsis, cp = cps, Ra = Ras, Ra_crit = Ra_crits, 
                        convect = convects, d0=d0s, t_diff = t_diffs, H=Hs)

print('Differentiation complete. It took', time.strftime("%Hh%Mm%Ss", time.gmtime(int_time1)))
#%%
######################## Thermal evolution ####################################

#integrate
tic = time.perf_counter()
Tc, Tc_conv, Tcmb, Tm_mid, Tm_conv, Tm_surf, Tprofile, f, Xs, dl, dc, d0,  \
    min_unstable, Ur, Ra, RaH, RanoH, Racrit, Fs, Flid, Fad, Fcmb, Rem, B, \
        buoyr, t, fcond_t = thermal_evolution(t_diff[-1],t_end,step_m,
                                              Tdiff[:,-1],f0,sparse_mat_c,sparse_mat_m) 
toc = time.perf_counter()
int_time2 = toc - tic    

#update on progress
if automated == False:
    plt.figure()
    plt.scatter(rplot,Tprofile[-1,:])
    plt.xlabel('r/km')
    plt.ylabel('Temperature/K')
    plt.title('Temperature profile post thermal evolution')

print('Thermal evolution complete', time.strftime("%Hh%Mm%Ss", time.gmtime(int_time2)))
#%%
############################# Process data ####################################
################ all processes which happen once  #############################
int_time = int_time1+int_time2 #total time for the two scripts
nmantle = int((r/dr)/2)
diff_time = t_diff[-1]/Myr
diff_T = Tdiff[int(nmantle),-1]
peakT = np.amax(Tprofile[:,nmantle+1:])
loc_max1 = np.where(Tprofile[:,nmantle+1:]==peakT)[0][0] #take the set of time coordinates and first value (they should all be the same)
tmax = t[loc_max1]/Myr
peak_coreT = np.amax(Tprofile[:,:nmantle])
loc_max2 = np.where(Tprofile[:,:nmantle]==peak_coreT)[0][0] #take the set of time coordinates and first value (they should all be the same)
tcoremax = t[loc_max2]/Myr
if np.any(f<f0):
    tsolid_start = t[f<f0][0]/Myr #start of core solidification
    tsolid = t[-1]/Myr #time of core solidification
else:
    tsolid_start = np.nan
    tsolid = np.nan

#erosion of stratification
tstrat = t[min_unstable!=0] #all times when the core is stratified

if len(tstrat)==0: #core never stratified
    tstrat_start = t[0]/Myr
    tstrat_remove = t[0]/Myr
    strat_end = t[0]/Myr
else:
    tstrat_start = tstrat[0]/Myr
    if tstrat[-1]==t[-1]: #stratification never removed
        tstrat_remove = np.inf
        strat_end = np.inf
    else:
        strat_end = tstrat[-1]/Myr+0.1 #end of stratification erosion
        max_strat = round(n_cells/2)-1 #max height of stratification
        if len(t[(min_unstable<max_strat)&(min_unstable>0)])==0: #stratification removed instantaneously
            tstrat_remove = strat_end
        else:
            tstrat_remove = t[(min_unstable<max_strat)&(min_unstable>0)][0]/Myr #beginning of stratification erosion

#switch to conduction

if np.any((t/Myr)<fcond_t):
    fcond_T = Tprofile[(t/Myr)<fcond_t,nmantle+1][-1] #temperature when mantle stops convecting
else:
    fcond_T = np.nan
    
# Frad - radiogenic heat flux, normalised to surface of body
from heating import al_heating
h = al_heating(t) 
Frad = h*rhom*Vm/As #radiogenic heatflux

#combine these in a single array
Flux = [Fs, Fcmb, Fad, Frad, Flid]
Fdrive = Fcmb - Fad        

################# Threshold for magnetic field being on ########################
threshold1 = 10 
threshold2 = 40
threshold3 = 100

########### Dynamo maxima - must be greater than threshold1 to be on ##########
#thermal dynamo
Remtherm = Rem[f>=f0]
Btherm = B[f>=f0]

if np.any(Remtherm>threshold1):
    max_Btherm = max(Btherm)
    max_Bthermt = t[B==max_Btherm][0]/Myr
else:
    max_Btherm = 0
    max_Bthermt = np.nan
    
max_Rtherm = max(Remtherm)
max_Rthermt = t[Rem==max_Rtherm][0]/Myr

#compositional - turn into df for rolling average
#no max time as averaging
Remcomp = pd.Series(Rem[f<f0])
Bcomp = pd.Series(B[f<f0])
wfrac = 30 #fractional window width
wn = int(len(Remcomp)/wfrac)
if wn < 1: #check if window is wider than data
    wn = int(len(Remcomp)/5)
    
Remav = Remcomp.rolling(window=wn,center=True).mean()
Bav = Bcomp.rolling(window=wn,center=True).mean()
if np.any(Remav[wn:-wn]>threshold1):
    max_Bcomp = max(Bav[wn:-wn]) #exclude nan
else:
    max_Bcomp = 0
    
max_Rcomp = max(Remav[wn:-wn]) 

########################## on and off times - calculate and save ####################
from duration_calc import on_off_test

#Rem > 10  
on, off, dur = on_off_test(t/Myr,Rem,threshold1,save_interval_mag/Myr) 
Bn1 = len(on) #number of on periods
if on.size > 0:
    if (on[0] < strat_end) & (off[0]>strat_end): #dynamo generation extends through stratification
        on[0] = strat_end
    elif (on[0] < strat_end) & (off[0]<strat_end): #dynamo only on before stratification
        on[0] = 0
        off[0] = 0
    else:
        pass        
    magon_1 = on[0]
    magoff_1 = off[0]
else:
    magon_1 = 0
    magoff_1 = 0    
if len(on) > 1:
    magon_2 = on[1]
    magoff_2 = off[1]
else:
    magon_2 = 0
    magoff_2 = 0
if len(on) > 2:
    magon_3 = on[2]
    magoff_3 = off[2]
else:
    magon_3 = 0
    magoff_3 = 0

#Rem > 40
on, off, dur = on_off_test(t/Myr,Rem,threshold2,save_interval_mag/Myr) 
Bn2 = len(on) #number of on periods
if on.size > 0: #i.e. there is one on value > nan
    if (on[0] < strat_end) & (off[0]>strat_end): #dynamo generation extends through stratification
        on[0] = strat_end
    elif (on[0] < strat_end) & (off[0]<strat_end): #dynamo only on before stratification
        on[0] = 0
        off[0] = 0
    else:
        pass 
    magon_4 = on[0]
    magoff_4 = off[0]
else:
    magon_4 = 0
    magoff_4 = 0
    
if len(on) > 1:
    magon_5 = on[1]
    magoff_5 = off[1]
else:
    magon_5 = 0
    magoff_5 = 0
    
# Rem > 100
on, off, dur = on_off_test(t/Myr,Rem,threshold3,save_interval_mag/Myr) 
Bn3 = len(on) #number of on periods
if on.size > 0:
    if (on[0] < strat_end) & (off[0]>strat_end): #dynamo generation extends through stratification
        on[0] = strat_end
    elif (on[0] < strat_end) & (off[0]<strat_end): #dynamo only on before stratification
        on[0] = 0
        off[0] = 0
    else:
        pass 
    magon_6 = on[0]
    magoff_6 = off[0]
else:
    magon_6 = 0
    magoff_6 = 0
    
if len(on) > 1:
    magon_7 = on[1]
    magoff_7 = off[1]
else:
    magon_7 = 0
    magoff_7 = 0

############################ Save results #####################################
# save variables to file
if full_save == True:
    np.savez_compressed(f'{folder}run_{run}', Tc = Tc, Tc_conv = Tc_conv, 
                        Tcmb = Tcmb,  Tm_mid = Tm_mid, Tm_conv = Tm_conv, 
                        Tm_surf = Tm_surf, T_profile = Tprofile, Flid = Flid, 
                        f=f, Xs = Xs, dl = dl, dc=dc, d0 = d0, min_unstable=min_unstable, 
                        Ur=Ur, Ra = Ra, RaH= RaH, RanoH = RanoH, Racrit = Racrit, 
                        t=t, Rem = Rem, B=B, buoyr = buoyr, Flux = Flux) 

if B_save == True:
    np.savez_compressed(f'{folder}run_{run}_B', B=B, Rem = Rem, t = t)
    
#write parameters to the run file
from csv import writer
  
var_list = [run,tsolid,int_time,diff_time, diff_T, peakT, tmax, peak_coreT, tcoremax, 
            tstrat_start, tstrat_remove, strat_end, fcond_t, fcond_T, tsolid_start, 
            max_Rtherm, max_Rthermt, max_Btherm, max_Bthermt, max_Rcomp, max_Bcomp, 
            Bn1, magon_1, magoff_1, magon_2, magoff_2, magon_3, magoff_3, Bn2, 
             magon_4, magoff_4, magon_5, magoff_5, Bn3, magon_6, magoff_6, magon_7, magoff_7]

with open(f'{folder}run_results.csv','a') as f_object:
     writer_object = writer(f_object) #pass file object to csv.writer
     writer_object.writerow(var_list) # pass list as argument into write row
     f_object.close() #close file

print('Results and run parameters saved. Run sucessful')

#add done flag to run
if automated == True: #no need to reimport ind as will have been imported earlier
    auto = pd.read_csv(f'{folder}auto_params.csv')
    auto.loc[ind+1,'status']=1 #indicates completed
    auto.to_csv(f'{folder}auto_params.csv',index=False)

