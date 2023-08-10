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
from parameters import  run, t_acc_m, t_end_m, dr, automated, Myr, Ts, f0, r, rc, kappa_c, save_interval_d, save_interval_t, km, Vm, As
from parameters import rhom, step_m, Xs_0, default, rcmf, Fe0, full_save, B_save, conv_tol, frht, eta0, etal, w, alpha_n

if automated == True: #should say true am just testing
    import sys
    folder = sys.argv[1]
else:
    folder = 'Results_combined/' #folder where you want to save the results

#set flag for run started
if automated == True:
    from parameters import ind
    auto = pd.read_csv(f'{folder}auto_params.csv')
    auto.loc[ind+1,'status']=0 #indicates started
    auto.to_csv(f'{folder}auto_params.csv',index=False)
else: #save run parameters in run_info file
    run_info = {"run":[run],"r":[r],"default":[default],"rcmf":[rcmf],"eta0":[eta0],"frht":[frht],"w":[w],"etal":[etal],"alpha_n":[alpha_n],"Xs_0":[Xs_0], "Fe0":[Fe0], "t_acc_m":[t_acc_m], "t_end_m":[t_end_m], "dr":[dr],"step_m":[step_m]}
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
n_cells = int(r/dr)+1 #number of cells needed to span body - add one to include the centre
Tint = np.ones([n_cells])*Ts#first element in the array is at r=0, accrete cold at surface temp 
Tint[-1]=Ts

print(f'Beginning run {run}')
print('Initial conditions set')

########################### Differentiation ###################################
tic = time.perf_counter()
Tdiff, Xfe, Xsi, cp, Ra, Ra_crit, convect, d0, t_diff, H  = differentiation(Tint,t_acc,r, dr, step_m)
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

if full_save == True:
    np.savez_compressed(f'{folder}run_{run}_diff', Tdiff = Tdiff, Xfe = Xfe, Xsi = Xsi, cp = cp, Ra = Ra, Ra_crit = Ra_crit, convect = convect, d0=d0, t_diff = t_diff, H=H)

print('Differentiation complete. It took', time.strftime("%Hh%Mm%Ss", time.gmtime(int_time1)))
######################## Thermal evolution ####################################

#integrate
tic = time.perf_counter()
Tc, Tc_conv, Tcmb, Tm_mid, Tm_conv, Tm_surf, Tprofile, f, Xs, dl, dc, d0, min_unstable, Ur, Ra, RaH, RanoH, RaRob, Racrit, Fs, Flid, Fad, Fcmb, Rem_c, Bcomp, t = thermal_evolution(t_diff[-1],t_end,step_m,Tdiff[:,-1],f0,sparse_mat_c,sparse_mat_m) 
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

############################# Process data ####################################
########## Current comparitive parameters ####################
from duration_calc import on_off_save

################ all processes which happen once  #############################
int_time = int_time1+int_time2 #total time for the two scripts
nmantle = int((r/dr)/2)
diff_time = t_diff[-1]/Myr
diff_T = Tdiff[int(nmantle),-1]
peakT = np.amax(Tprofile[:,nmantle:])
loc_max = np.where(Tprofile[:,nmantle:]==peakT)[0][0] #take the set of time coordinates and first value (they should all be the same)
tmax = t[loc_max]/Myr
peak_coreT = np.amax(Tprofile[:,:nmantle])
loc_max = np.where(Tprofile[:,:nmantle]==peak_coreT)[0][0] #take the set of time coordinates and first value (they should all be the same)
tcoremax = t[loc_max]/Myr
tsolid = t[-1]/Myr #time of core solidification
if np.all(Tprofile[:,int(nmantle)-2]<Tcmb):
    tstrat_remove = np.inf
else:
    tstrat_remove = t[Tcmb < Tprofile[:,int(nmantle)-2]][0]/Myr
    
if np.any(min_unstable==0):
     strat_end = t[np.where(min_unstable==0)[0]][0]/Myr
else:
    strat_end = np.inf

#switch to conduction
if np.any(Ra/Racrit<0.5):
    fcond_t = t[Ra/Racrit<0.5][0]/Myr #half the critical value (end of buffering)
    if np.any(d0>(r-rc)):
        fcond_t2 = t[d0>(r-rc)][0]/Myr #check convection not ended by stagnant lid thickening
        if fcond_t2 < fcond_t:
            fcond_t = fcond_t2 #lid thickening shuts off convection first
else:
    fcond_t = np.nan
if np.any(Ra/Racrit>(2-conv_tol)):
    lconv_t = t[Ra/Racrit>(2-conv_tol)][-1]/Myr #last supercritical time (start of buffering)
    lconv_T = Tprofile[Ra/Racrit>(2-conv_tol),nmantle+1][-1] #temperature when Ra first starts buffering
else:
    lconv_t = np.nan
    lconv_T = np.nan
    
# Frad - radiogenic heat flux, normalised to surface of body
from heating import Al_heating
h = Al_heating(t) 
Frad = h*rhom*Vm/As #radiogenic heatflux

#combine these in a single array
Flux = [Fs, Fcmb, Fad, Frad]

# calculate thermal magnetic reynolds number
from Rem_calc import Rem_therm, B_flux_therm
# need Fdrive for Rem_therm
#only calculate this for Fdrive >0
Fdrive = Fcmb - Fad
Fdrive_nn = Fdrive.copy()
Fdrive_nn[Fdrive<0]=0
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
        
# maximum Rem - ignore first entry as Rem_mac is duplicated
max_Rem = np.zeros([m-1])
max_Remt = np.zeros([m-1])
for i in range(m-1):
    if np.any(Rem[i+1,:]>threshold):
        max_Rem[i] = np.max(Rem[i+1,Rem[i+1,:]>threshold])
        max_Remt[i] = t[Rem[i+1,:]==max_Rem[i]][0]/Myr
        
########################## on and off times - calculate and save ####################
#set 10*save_interval (1Myr) as on off tolerance to smooth on-off periods smaller than 1Myr
t_plot_t = t/Myr

mac_on, mac_off, mac_dur, mac_n, mac_ngl10, mac_ngl100 = on_off_save(t_plot_t, Rem_t[0],threshold,save_interval_t/Myr, f'{folder}MAC_onoff.csv', 'MAC', run) #MAC on off
cia_on, cia_off, cia_dur, cia_n, cia_ngl10, cia_ngl100 = on_off_save(t_plot_t, Rem_t[1],threshold,10*save_interval_t/Myr, f'{folder}CIA_onoff.csv', 'CIA', run) #CIA on off
comp_on, comp_off, comp_dur, comp_n, comp_ngl10, comp_ngl100 = on_off_save(t_plot_t, Rem_c, threshold,10*save_interval_t/Myr, f'{folder}comp_onoff.csv', 'comp', run) #comp on off
coreconv_on, coreconv_off, coreconv_dur, coreconv_n, coreconv_ngl10, coreconv_ngl100 = on_off_save(t_plot_t, Fdrive,0,10*save_interval_t/Myr, f'{folder}coreconv_onoff.csv', 'core_conv', run) #core convection on off

############################ Save results #####################################
# save variables to file
if full_save == True:
    np.savez_compressed(f'{folder}run_{run}', Tc = Tc, Tc_conv = Tc_conv, Tcmb = Tcmb,  Tm_mid = Tm_mid, Tm_conv = Tm_conv, Tm_surf = Tm_surf, 
             T_profile = Tprofile, Flid = Flid, f=f, Xs = Xs, dl = dl, dc=dc, d0 = d0, min_unstable=min_unstable, Ur=Ur, 
             Ra = Ra, RaH= RaH, RanoH = RanoH, RaRob = RaRob, Racrit = Racrit, t=t, Rem_t = Rem_t, B = B, Rem_c = Rem_c, Flux = Flux) 

if B_save == True:
    np.savez_compressed(f'{folder}run_{run}_B', t=t, Rem_t = Rem_t, B = B, Rem_c = Rem_c)
    
#write parameters to the run file
from csv import writer
  
var_list = [run,tsolid,int_time,diff_time, diff_T, peakT, tmax, peak_coreT, tcoremax, tstrat_remove, 
             strat_end, fcond_t, lconv_t,lconv_T, max_Rem[0], max_Remt[0], max_Rem[1], max_Remt[1], max_Rem[2], 
             max_Remt[2],max_B[0],max_Bt[0],max_B[1],max_Bt[1],max_B[2],max_Bt[2],max_B[3],max_Bt[3], mac_on, 
             mac_off, mac_dur, mac_n, mac_ngl10, mac_ngl100, cia_on, cia_off, cia_dur, cia_n, cia_ngl10, 
             cia_ngl100,comp_on, comp_off, comp_dur, comp_n, comp_ngl10, comp_ngl100,coreconv_on, coreconv_off, coreconv_dur, coreconv_n, coreconv_ngl10, coreconv_ngl100]

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

