#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for solving the thermal evolution of an asteroid. Doesn't include differentiation (i.e commented out)'
"""
# import modules
import numpy as np
import matplotlib.pyplot as plt
import time #use this to time the integration

#import time constants and initial conditions
from parameters import  Myr, Ts, f0, r, rc, dr, kappa_c, save_interval_d, save_interval_t, km, Vm, As, rhom, default, rcmf, Xs_0, Fe0, t_cond_core


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
run = 8

t_acc=0.8*Myr  #Accretion time
t_end_m=10#end time in Myr

t_end=t_end_m*Myr
step_m=0.1*t_cond_core  #max timestep must be smaller than conductive timestep
n_save_d = int(save_interval_d/step_m)
n_save_t = int(save_interval_t/step_m)

# set initial temperature profile
n_cells = int(r/dr) #number of cells needed to span body
Tint = np.ones([n_cells])*Ts#first element in the array is at r=0, accrete cold at surface temp 
Tint[-1]=Ts

print('Initial conditions set')

###############  sintering and compaction code will go here  ##################

########################### Differentiation ###################################
tic = time.perf_counter()
Tdiff, Xfe, Xsi, cp, Ra, Ra_crit, convect, d0, t_diff, H  = differentiation(Tint,t_acc,r, dr, step_m)
toc = time.perf_counter()
diff_time = toc - tic  

# update user on progress and plot differentiated temperature profile 
rplot= np.arange(0,r,dr)/1e3

plt.figure()
plt.scatter(rplot,Tdiff[:,-1])
plt.xlabel('r/km')
plt.ylabel('Temperature/K')
plt.title('Temperature profile post differentiation')

#rescale data and save here in case thermal evolution crashes
Tdiff = Tdiff[:,0::n_save_d]
Xfe = Xfe[:,0::n_save_d]
Xsi = Xsi[:,0::n_save_d]
cp = cp[:,0::n_save_d]
Ra = Ra[0::n_save_d]
Ra_crit = Ra_crit[0::n_save_d]
convect = convect[0::n_save_d]
d0 = d0[0::n_save_d]
t_diff = t_diff[0::n_save_d]
H = H[0::n_save_d]

np.savez(f'Results_combined/run_{run}_diff', Tdiff = Tdiff, Xfe = Xfe, Xsi = Xsi, cp = cp, Ra = Ra, Ra_crit = Ra_crit, convect = convect, d0=d0, t_diff = t_diff, H=H)

print('Differentiation complete. It took', time.strftime("%Hh%Mm%Ss", time.gmtime(diff_time)))
######################## Thermal evolution ####################################

#integrate
tic = time.perf_counter()
Tc, Tc_conv, Tcmb, Tm_mid, Tm_conv, Tm_surf, Tprofile, f, Xs, dl, dc, d0, min_unstable, Ur, Ra, RaH, RanoH, RaRob, Racrit, Fs, Flid, Fad, Fcmb, Rem_c, Bcomp, t, cond_t = thermal_evolution(t_diff[-1],t_end,step_m,Tdiff[:,-1],f0,sparse_mat_c,sparse_mat_m) 
toc = time.perf_counter()
int_time = toc - tic    

#update on progress
plt.figure()
plt.scatter(rplot,Tprofile[-1,:])
plt.xlabel('r/km')
plt.ylabel('Temperature/K')
plt.title('Temperature profile post thermal evolution')

print('Thermal evolution complete', time.strftime("%Hh%Mm%Ss", time.gmtime(int_time)))

############################# Process data ####################################
########## Current comparitive parameters ####################

nmantle = int((r/dr)/2)
diff_time = t_diff[-1]/Myr
diff_T = Tdiff[int(nmantle),-1]
peakT = np.amax(Tprofile[:,nmantle:])
loc_max = np.where(Tprofile[:,nmantle:]==peakT)[1][0] #take the set of time coordinates and first value (they should all be the same)
tmax = t[loc_max]/Myr
if np.all(Tprofile[:,int(nmantle)-2]<Tcmb):
    tstrat_remove = np.inf
else:
    tstrat_remove = t[Tcmb < Tprofile[:,int(nmantle)-2]][0]/Myr
    
if np.any(min_unstable==0):
     strat_end = t[np.where(min_unstable==0)[0]][0]/Myr
else:
    strat_end = np.inf
    
if np.all(Fcmb < Fad):
    super_ad_start = np.inf
    super_ad_end = np.inf
else:
    super_ad_start = t[np.where(Fcmb>Fad)[0]][0]/Myr
    super_ad_end = t[np.where(Fcmb>Fad)[0]][-1]/Myr

# Frad - radiogenic heat flux, normalised to surface of body
from heating import Al_heating
h = Al_heating(t) 
Frad = h*rhom*Vm/As #radiogenic heatflux

#combine these in a single array
Flux = [Fs, Fcmb, Fad, Frad]

# calculate thermal magnetic reynolds number
from Rem_calc import Rem_therm, B_flux_therm
from parameters import Xs_eutectic
# need Fdrive for Rem_therm
#only calculate this for Fdrive >0
Fdrive = Fcmb - Fad
Fdrive_nn = Fdrive.copy()
Fdrive_nn[Fdrive<0]=0
Fdrive_nn[Xs>=Xs_eutectic]=0 #no dynamo in eutectic solidification
Rem_t = Rem_therm(Fdrive_nn,f,min_unstable) # magnetic Reynolds number for thermal convection - tuple of MAC and CIA balance
Bml, Bmac, Bcia = B_flux_therm(Fdrive_nn,f,min_unstable) # field strength for thermal convection based on energy flux scaling 
B = np.array([Bml, Bmac, Bcia, Bcomp])
Rem = np.array([Rem_t[0],Rem_t[0],Rem_t[1],Rem_c]) #use conservative MAC Rem for ML

#calculate maximum field strengths
m = 4
threshold = 10 
max_B = np.zeros([m])
max_Bt = np.zeros([m])
for i in range(m):
    if np.any(Rem[i,:]>threshold):
        max_B[i] = np.max(B[i,Rem[i,:]>threshold])
        max_Bt[i] = t[B[i,:]==max_B[i]][0]/Myr
        
# maximum Rem
max_Rem = np.zeros([m-1])
max_Remt = np.zeros([m-1])
for i in range(1,m):
    if np.any(Rem[i,:]>threshold):
        max_Rem[i] = np.max(Rem[i,Rem[i,:]>threshold])
        max_Remt[i] = t[Rem[i,:]==max_Rem[i]][0]/Myr
        
#on and off times
t_plot_t = t/Myr
t_therm1 = t_plot_t[Rem[0]>threshold]
t_therm11 = t_therm1[t_therm1<100] #distinguish between early and late dynamo
t_therm12 = t_therm1[t_therm1>100]
t_therm2 = t_plot_t[Rem[2]>10]
t_therm21 = t_therm2[t_therm2<100]
t_therm22 = t_therm2[t_therm2>100]
t_comp = t_plot_t[Rem_c>10]
if np.any(Rem[0]>threshold):
    if np.any(t_therm1<100):
        MAC_start = t_therm11[0]
        MAC_stop = t_therm11[-1]
        print(f'MAC balance predicts there is a thermal dynamo between {MAC_start:.1f} and {MAC_stop:.1f}Myr')
    if np.any(t_therm1>100):
        MAC_start = t_therm12[0]
        MAC_stop = t_therm12[-1]
        print(f'MAC balance predicts there is a thermal dynamo between {MAC_start:.1f} and {MAC_stop:.1f}Myr')
else:
    MAC_start = np.nan
    MAC_stop = np.nan
    print('MAC balance predicts there is no thermal dynamo')
if np.any(Rem[2]>threshold):
    if np.any(t_therm2<100):
        CIA_start = t_therm21[0]
        CIA_stop = t_therm21[-1]
        print(f'CIA balance predicts there is a thermal dynamo between {CIA_start:.1f} and {CIA_stop:.1f}Myr')
    if np.any(t_therm2>threshold):
        CIA_start = t_therm22[0]
        CIA_stop = t_therm22[-1]
        print(f'CIA balance predicts there is a thermal dynamo between {CIA_start:.1f} and {CIA_stop:.1f}Myr')
else:
    CIA_start = np.nan
    CIA_stop = np.nan
    print('CIA balance predicts there is no thermal dynamo')
if np.any(Rem_c>threshold):
    comp_start = t_comp[0]
    comp_stop = t_comp[-1]
    print(f'There is a compositional dynamo between {comp_start:.1f} and {comp_stop:.1f} Myr')
else:
    comp_start = np.nan
    comp_stop = np.nan
    print('There is no compositional dynamo')
print('Fluxes and magnetic Reynolds number calculated.')

############################ Save results #####################################
# save variables to file
np.savez('Results_combined/run_{}'.format(run), Tc = Tc, Tc_conv = Tc_conv, Tcmb = Tcmb,  Tm_mid = Tm_mid, Tm_conv = Tm_conv, Tm_surf = Tm_surf, 
         T_profile = Tprofile, Flid = Flid, f=f, Xs = Xs, dl = dl, dc=dc, d0 = d0, min_unstable=min_unstable, Ur=Ur, 
         Ra = Ra, RaH= RaH, RanoH = RanoH, RaRob = RaRob, Racrit = Racrit, t=t, Rem_t = Rem_t, B = B, Rem_c = Rem_c, Flux = Flux) 

#write parameters to the run file
from csv import writer

var_list = [run, r, dr, t_acc/Myr, t_end_m, step_m/Myr, max(t)/Myr, cond_t, int_time,  save_interval_d/Myr, save_interval_t/Myr, 
            default, rcmf, Xs_0, Fe0]

    
with open('Results_combined/run_info.csv','a') as f_object:         
    writer_object = writer(f_object) #pass file object to csv.writer
    writer_object.writerow(var_list) # pass list as argument into write row
    f_object.close() #close file
  
var_list2 = [run,step_m,dr,diff_time, diff_T, peakT, tmax, tstrat_remove, 
             strat_end, super_ad_start, super_ad_end, cond_t, max_Rem[0], max_Remt[0], max_Rem[1], max_Remt[1], max_Rem[2], 
             max_Remt[2],MAC_start,MAC_stop, CIA_start, CIA_stop, comp_start, comp_stop,
             max_B[0],max_Bt[0],max_B[1],max_Bt[1],max_B[2],max_Bt[2],max_B[3],max_Bt[3]]
from csv import writer
with open('Results_combined/timestep_test.csv','a') as f_object:
     writer_object = writer(f_object) #pass file object to csv.writer
     writer_object.writerow(var_list2) # pass list as argument into write row
     f_object.close() #close file
# print('Results and run parameters saved')
