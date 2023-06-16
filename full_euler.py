#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for solving thermal evolution for the entire evolution of the body
By convention Tprofile[int(n_cells/2)] is last core cell
Flow:
    1. Obtain conductive profile for whole body
    2. Calculate stagnant lid thickness and Rayleigh number
    3. If Ra > Rac, replace rc < r < R-d0 with convective profile for mantle
    4. Replace r < rc with isothermal core profile
    5. Calculate df/dt and f    
"""
#import modules
import numpy as np
from parameters import Ts, Myr, dr, out_interval, save_interval_t, km, kc, alpha_c, r, rc, rhoc, gc, Vm, rhom, As
from parameters import cpc, Xs_0, default, Xs_eutectic, Acmb, Lc, Pc

#import required functions
from T_cond import Tm_cond_calc, Tc_cond_calc
from dTmdt_def import dTmdt_calc
from dTcdt_def import dTcdt_calc, dTcdt_calc_solid 
from Rayleigh_def import Rayleigh_calc, Rayleigh_crit
from cmb_bl import delta_l, delta_c
from q_funcs import Qr
from fe_fes_liquidus import fe_fes_liquidus_bw 
from stratification import volume_average
from heating import Al_heating

def thermal_evolution(tstart,tend,dt,T0,f0,sparse_mat_c,sparse_mat_m):
    """
    

    Parameters
    ----------
    tstart : float  
        start time [s]
    tend : float
        end time [s]
    dt : float
        time step [s]
    T0 : array
        initial temperature profile [K]
    f0: float
        initial fractional inner core radius
    sparse_mat_c: sparse matrix
        stencil for conductive temperature evolution of core
    sparse_mat_m: sparse matrix
        stencil for conductive temperature evolution of mantle
        
    Returns
    -------
    Tc: array
        core temperature, measured at centre [K]
    Tc_conv: array
        temperature of convecting core (0 if not convecting) [K]
    Tcmb: array
        CMB temperature [K]
    Tm_mid: array
        mantle temp in the middle cell [K]
    Tm_conv: array
        temperature of convecting mantle (0 if not convecting) [K]
    Tm_surf: array
       mantle temperature one cell below the surface [K]
    Tprofile: array
        temperature profile as a function of time, output every 10Myr
    f: array
        fractional inner core radius [m]
    Xs: array
        wt % S in core 
    dl: array
        thickness of mantle CMB boundary layer (0 if mantle not convecting) [m]
    dc: array
        thickness of core CMB boundary layer (0 if core not thermally convecting) [m]
    d0: array
        stagnant lid thickness for convecting mantle [m]
    min_unstable : int
        index of base of convecting core
    Ur : float
        Urey ratio 
    Ra: array
        Rayleigh number for convecting mantle
    RaH: array
        radiogenic Rayleigh number for convecting mantle
    RanoH: array
        non-radiogenic Rayleigh number for convecting mantle
    RaRob : array
        modified Rayleigh Roberst number for convecting mantle (eqn 25 of Sturtz 2022b)
    Racrit: array
        critical Rayleigh number for convecting mantle
    Flid: array
        heat flux across base of stagnant lid [W m^-2]
    Fs: array
        surface heat flux [W m^-2]
    Fad: array
        adiabatic CMB heat flux [W m^-2]
    Fcmb: array
        CMB heat flux [W m^-2]
    Rem_c : array
        compositional magnetic Reynolds number
    Bcomp : float
        compositional magnetic field strength [T]
    tsolve: array
        time points corresponding to each of the values above [s]
    cond_t: float
        time when the mantle switched from conduction to convection, nan if didn't switch [s]
         
    """

    #initialise arrays for output
    m = round((tend-tstart)/save_interval_t)+1 #add one so always enough space
    n_cells = len(T0) #number of cells
    i_core = round(n_cells/2)-1 # index in array of last core cell (-1 as indexing starts at 0)
    cond_t = 'nan' #set as string and then reset if switches to conduction
    mantle_conv = False #flag for mantle convection
    core_conv = False #flag for core convection

    
    #output variables
    Xs = np.ones([m])*Xs_0 #core sulfur fraction
    Ra = np.zeros([m])
    Racrit = np.zeros([m])
    RaH = np.zeros([m])
    RanoH = np.zeros([m])
    RaRob = np.zeros([m])
    d0 = np.zeros([m])
    dl = np.zeros([m])
    dc = np.zeros([m])
    min_unstable = np.ones([m],dtype=int)*(i_core-1) #smallest index of cells in the core that are convectively unstable - as a minimum it is the one below the CMB
    Ur = np.zeros([m])
    Tprofile= np.zeros([m,n_cells])
    Tc = np.zeros([m])
    Tc_conv = np.zeros([m])
    Tcmb = np.zeros([m])
    Tm_conv = np.zeros([m])
    Tm_mid = np.zeros([m])
    Tm_surf = np.zeros([m])
    f = np.ones([m])*f0 #set as initial core size by default and only override if core starts to change
    Fs = np.zeros([m])
    Flid = np.zeros([m])
    Fad = np.zeros([m])
    Fcmb = np.zeros([m])
    Rem_c = np.zeros([m])
    Bcomp = np.zeros([m])
    tsolve = np.zeros([m])

    
    #Step 0. Calculate time, get two separate temperature arrays
    # the last cell of the core array is the same as the first cell of the mantle array
    tsolve_new = tstart + dt
    T0_core = T0[:i_core+1] #include last core cell
    T0_mantle = T0[i_core:]
    nmantle_cells = len(T0_mantle)
    ncore_cells = len(T0_core)
    i=0 #counter for saving and printing
    
    # Initialise values that might not get calculated
    f_new = f0
    Xs_new = Xs_0
    dc_new = 0
    dl_new = 0
    min_unstable_new = int((i_core-1))

    
    ##################    Initial step    #########################
    # Step 1. Calculate conductive profile for mantle
    T_new_mantle = Tm_cond_calc(tsolve_new,dt,T0_mantle,sparse_mat_m)
    dTdt_mantle_new = (T_new_mantle-T0_mantle)[1] #default use temp change at base above CMB
    Fs_new = -km*(T_new_mantle[-1]-T_new_mantle[-2])/dr
    Flid_new = Fs_new #by default lid flux is same as surface
    h = Al_heating(tsolve_new)
    Ur_new = rhom*Vm*h/abs(Fs_new*As) #calculate Urey ratio
    
    # Step 2. Is the mantle convecting? Calculate stagnant lid thickness, base thickness and Rayleigh number
    Ra_new, d0_new, RaH_new, RanoH_new, RaRob_new = Rayleigh_calc(tsolve_new,T0_mantle[1],dTdt_mantle_new,Ur_new,default) #use temp at base of mantle 
    Racrit_new = Rayleigh_crit(T0_mantle[1])   
      
    # in initial step use balance of eqn 23 and 24 to find Tcmb
    Tcmb_new = (km/kc*T0_mantle[0]+T0_core[0])/(km/kc+1)
    T_new_mantle[0] = Tcmb_new
    Fcmb_new = -km*(T0_mantle[0]-Tcmb_new)/dr # CMB heat flux eqn 23 in Dodds 2020
    
    if (Ra_new < Racrit_new) | (d0_new > (r-rc)):# not convecting
        Tm_conv_new = 0 # convective mantle temperature is 0 if mantle not convecting
       
    else: #mantle is convecting replace mantle below stagnant lid with isothermal convective profile
        dl_new = delta_l(T0_mantle[1],T0_mantle[0],Ur_new)
        nlid_cells = round(d0_new/dr)
        if nlid_cells ==0:
            lid_start = nmantle_cells -2
        else:
            lid_start = nmantle_cells - nlid_cells - 1 #index in temp array where lid starts
        if d0_new < dr:
            Flid_new = -km*(Ts-T_new_mantle[lid_start])/d0_new #if less than grid thickness choose d0 so don't overestimate thickness
            Fs_new = Flid_new #lid determines flux out of surface
        else:
            Flid_new = -km*(T_new_mantle[lid_start+1]-T_new_mantle[lid_start])/dr 
        dTdt_mantle_new = dTmdt_calc(tsolve_new,T0_mantle[1],d0_new,Flid_new,Fcmb_new)
        Tm_conv_new = T0_mantle[1] + dTdt_mantle_new*dt #take temp one above CMB
        T_new_mantle[:lid_start+1] = Tm_conv_new
        T_new_mantle[-1] = Ts #pin surface to 200K

        
        
    #store values   
    Tm_mid_new = T_new_mantle[round(nmantle_cells/2)] # temperature at the mid mantle
    Tm_surf_new = T_new_mantle[-2] #temperature one cell above surface pinned to 200K
                
    # Step 3. Calculate conductive profile for the core
    T_new_core = Tc_cond_calc(tsolve_new,dt,T0_core,sparse_mat_c,True)
    Tc_conv_new = 0 #0 by default 
    Rem_c_new = 0 # by default
    Bcomp_new = 0
    
    Fad_new = kc*T0_core[-2]*alpha_c*gc/cpc
    # Step 4. Is the core solidifying? 
    # is the core solidifying?
    Tliquidus = fe_fes_liquidus_bw(Xs_0,Pc)
    if T0_core[-2] < Tliquidus: #core solidifies outside in - convecting so isothermal beneath CMB
        core_conv = False
        if Xs_0>= Xs_eutectic:
            dTcdt = 0 # whilst undergoing eutectic solidification there is no temp change
            dfdt = -(Fcmb_new*Acmb-Qr(tsolve_new))/(4*np.pi*rc**3*f0**2*Lc*rhoc)                    
            f_new = f0 + dfdt*dt
            Xs_new = Xs_0 #sulfur concentration unchanged in eutectic solidification
            #find new convective temperature
            Tc_conv_new = T0_core[0] 
            T_new_core[:] = Tc_conv_new #replace everything with the convective temperature
            Rem_c_new = 0
            Bcomp_new = 0
        else:

            min_unstable_new = 0               
            dTcdt, f_new, Rem_c_new, Bcomp_new = dTcdt_calc_solid(tsolve_new,Fcmb_new, T0_core, f0, Xs_0, dt) 
            #find new convective temperature
            Tc_conv_new = T0_core[0] + dTcdt*dt 
            T_new_core[:] = Tc_conv_new #replace everything with the convective temperature
            Xs_new = Xs_0/(f_new**3) #update sulfur content
            
    # Is the core convecting  without solidification?
    elif (Fcmb_new > Fad_new) and (min_unstable_new==0): #super adiabatic and no stratification, core convects
        
        core_conv = True
        min_unstable_new = 0
        nbl_cells = round(dc_new/dr)
        bl_start = ncore_cells - nbl_cells - 1 #index in temp array where lid starts

        dTcdt = dTcdt_calc(tsolve_new,Fcmb_new, T0_core, f0)
        
        #find new convective temperature
        Tc_conv_new = T0_core[bl_start-1] + dTcdt*dt 
        
        T_new_core[:bl_start] = Tc_conv_new #replace everything with the convective temperature up to b.l
        
    else: # don't have whole core convection, for first step don't check for stratification as haven't recalculated Tcmb yet
        pass
        
    
    # Step 5. Recalculate the CMB temperature so now it is determined by the end of the step as it will be in subsequent steps
    #use balance of eqn 23 and 24 to find Tcmb
    Tcmb_new = (km/kc*T_new_mantle[1]+T_new_core[-2])/(km/kc+1)
    T_new_mantle[0] = Tcmb_new
    T_new_core[-1] = Tcmb_new
    Fcmb_new = -km*(T_new_mantle[0]-Tcmb_new)/dr # CMB heat flux eqn 23 in Dodds 2020
    Tc_new = T_new_core[0] #temperature of core is always taken at centre
    
    # Step 6. Replace old array with new ready for next step
    T_old_core = T_new_core
    T_old_mantle = T_new_mantle
    Tprofile_old = np.hstack((T_new_core,T_new_mantle[1:]))
    Tc_old = Tc_new
    Tc_conv_old = Tc_conv_new
    Tcmb_old = Tcmb_new
    Tm_mid_old = Tm_mid_new
    Tm_conv_old = Tm_conv_new
    Tm_surf_old = Tm_surf_new
    dTdt_mantle_old = dTdt_mantle_new
    f_old = f_new
    Xs_old = Xs_new
    dl_old = dl_new
    dc_old = dc_new
    d0_old = d0_new
    min_unstable_old = min_unstable_new
    Ur_old = Ur_new
    Ra_old = Ra_new
    RaH_old = RaH_new
    RanoH_old = RanoH_new
    RaRob_old = RaRob_new
    Racrit_old = Racrit_new
    Flid_old = Flid_new
    Fs_old = Fs_new
    Fad_old = Fad_new
    Fcmb_old = Fcmb_new
    tsolve_old = tsolve_new
    
    while tsolve_new < tend:

        #Step 0. Calculate time
        tsolve_new = tsolve_old + dt
        
        # Step 1. Calculate conductive profile for mantle
        T_new_mantle = Tm_cond_calc(tsolve_new,dt,T_old_mantle,sparse_mat_m)
        Fs_new = -km*(T_new_mantle[-1]-T_new_mantle[-2])/dr
        Flid_new = Fs_new #by default lid flux is same as surface
        dTdt_mantle_new = (T_new_mantle-T_old_mantle)[1] #default use temp change at base just above CMB
        
        # Step 2. Is the mantle convecting? Calculate stagnant lid thickness, base thickness and Rayleigh number
        Ra_new, d0_new, RaH_new, RanoH_new, RaRob_new = Rayleigh_calc(tsolve_new,T_old_mantle[1],dTdt_mantle_old,Ur_old,default) #use temp at base of mantle 
        Racrit_new = Rayleigh_crit(T_old_mantle[1])
       
        if d0_new<(r-rc): #other convection criteria are meaningless if lid thickness is greater than mantle thickness
            nlid_cells = round(d0_old/dr)
            if nlid_cells ==0:
                lid_start = nmantle_cells -2
            else:
                lid_start = nmantle_cells - nlid_cells - 1  #index in temp array where lid starts
                          
            dTdt_mantle_conv = dTmdt_calc(tsolve_old,T_old_mantle[1],d0_old,Flid_old,Fcmb_old) #convective dTdt - use base mantle temp as convective temp
            
            if (dTdt_mantle_conv<0): #mantle is cooling
                if (Ra_new <= Racrit_new) | (Tm_conv_old ==0): #check for Ra and switch to conduction
                    mantle_conv = False
                    Tm_conv_new = 0 # convective mantle temperature is 0 if mantle not convecting
                    
                else: #mantle is convecting replace mantle below stagnant lid with isothermal convective profile 
                    mantle_conv = True
                    dTdt_mantle_new = dTdt_mantle_conv
                    Tm_conv_new = T_old_mantle[1] + dTdt_mantle_new*dt #temperature of convecting region 
                    T_new_mantle[:lid_start+1] = Tm_conv_new
                    T_new_mantle[-1] = Ts #pin surface to 200K in case d0 < 1 cell thick
                    if d0_new < dr:
                        Flid_new = -km*(Ts-T_new_mantle[lid_start])/d0_new #flux from convecting region to stagnant lid
                        Fs_new = Flid_new #lid determines flux out of surface
                    else:
                        Flid_new = -km*(T_new_mantle[lid_start+1]-T_new_mantle[lid_start])/dr
                        
            else: #mantle still heating up
                if Ra_new <= Racrit_new: #check Ra
                    mantle_conv = False
                    Tm_conv_new = 0 # convective mantle temperature is 0 if mantle not convecting
    
                else: #mantle is convecting replace mantle below stagnant lid with isothermal convective profile 
                    mantle_conv = True               
                    dTdt_mantle_new = dTdt_mantle_conv
                    Tm_conv_new = T_old_mantle[1] + dTdt_mantle_new*dt #temperature of convecting region 
                    T_new_mantle[:lid_start+1] = Tm_conv_new
                    T_new_mantle[-1] = Ts #pin surface to 200K in case d0 < 1 cell thick
                    if d0_new < dr:
                        Flid_new = -km*(Ts-T_new_mantle[lid_start])/d0_new #flux from convecting region to stagnant lid
                        Fs_new = Flid_new #lid determines flux out of surface
                    else:
                        Flid_new = -km*(T_new_mantle[lid_start+1]-T_new_mantle[lid_start])/dr
             
        else:
            mantle_conv = False
            Tm_conv_new = 0

        if (Tm_conv_old !=0) & (Tm_conv_new ==0): #record time for switch to conduction
            cond_t = tsolve_new
            print(f"Mantle conducting, time is {tsolve_new/Myr}")
            
            
        # Calculate Urey ratio with correct heat flux
        h = Al_heating(tsolve_new)
        Ur_new = rhom*Vm*h/abs(Fs_new*As)
        
        #store values
        Tm_mid_new = T_new_mantle[round(nmantle_cells/2)] # temperature at the mid mantle
        Tm_surf_new = T_new_mantle[-2] #temperature one cell above surface pinned to 200K
          
        # Step 3. Calculate conductive profile for the core
        T_new_core = Tc_cond_calc(tsolve_new,dt,T_old_core,sparse_mat_c,True)
        Tc_conv_new = 0 #by default, over write if core convects
        Rem_c_new = 0 #by default assume no compositional convection
        Bcomp_new = 0
        f_new = f_old #by default overwrite if solidifies
        Xs_new = Xs_old
        min_unstable_new = min_unstable_old #continuity of mixed layer thickness by default
        
        # Step 4. Is the core solidifying? 
        # is the core solidifying?
        Tliquidus = round(fe_fes_liquidus_bw(Xs_old,Pc)) #round as Buono & Walker data to nearest integer
        if T_old_core[-2] < Tliquidus: #core solidifies from outside in- convecting so isothermal beneath CMB
            core_conv = False #core convects but b.l. thickness set by rho not T so use conductive Fcmb
            if Xs_old>= Xs_eutectic:
                dTcdt = 0 # whilst undergoing eutectic solidification there is no temp change
                dfdt = -(Fcmb_old*Acmb-Qr(tsolve_new))/(4*np.pi*rc**3*f_old**2*Lc*rhoc)                    
                f_new = f_old + dfdt*dt
                Xs_new = Xs_old #sulfur concentration unchanged in eutectic solidification
                #find new convective temperature
                Tc_conv_new = T_old_core[0] 
                T_new_core[:] = Tc_conv_new #replace everything with the convective temperature
                
            else:
                min_unstable_new = 0               
                dTcdt, f_new, Rem_c_new, Bcomp_new = dTcdt_calc_solid(tsolve_new,Fcmb_old, T_old_core, f_old, Xs_old, dt) 
                #find new convective temperature
                Tc_conv_new = T_old_core[0] + dTcdt*dt 
                T_new_core[:] = Tc_conv_new #replace everything with the convective temperature
                Xs_new = Xs_0/(f_new**3) #update sulfur content
                
            
        
        # Is the core convecting  without solidification?
        elif (Fcmb_old > Fad_old) and (min_unstable_old==0): #super adiabatic and no stratification, core convects
            core_conv = True
            min_unstable_new = 0
            nbl_cells = round(dc_old/dr)
            bl_start = ncore_cells - nbl_cells - 1 #index in temp array where lid starts

            dTcdt = dTcdt_calc(tsolve_new,Fcmb_old, T_old_core, f_old)
            
            #find new convective temperature
            Tc_conv_new = T_old_core[bl_start-1] + dTcdt*dt 
            
            T_new_core[:bl_start] = Tc_conv_new #replace everything with the convective temperature up to b.l
            
        else: # don't have whole core convection, 
            #check if there is thermal stratification
            core_conv = False 
            
            if min_unstable_old>0:  #there is thermal stratification 
                
                if np.all(T_old_core[:-1] < Tcmb_old):
                    min_unstable_new = min_unstable_old  #continuity of stratification
                        # scenario 1 - just conduction in the core
                        # use already calculated condctive profile and keep core in current state

                else: #scenario 2 - erosion of stratification, convective layer at top of core
                    core_conv = True
                    b_ind = np.where( T_old_core[:-1] >= Tcmb_old)[0] #indices of unstable layer as array
                    min_unstable_new = b_ind[0]
                    
                    dTcdt = dTcdt_calc(tsolve_new,Fcmb_old, T_old_core, f_old, stratification = [True, min_unstable_old])
                    Tc_conv_new = T_old_core[min_unstable_old]+dTcdt*dt #replace convecting layer from last timestep with new temp - in later steps use i-1 and i
                    T_new_core[min_unstable_old:-1] = Tc_conv_new
                                             
                    #now perform volume average over unstable layer
                    Tlayer = volume_average(T_new_core, b_ind,dr)
                    T_new_core[min_unstable_new:-1] = Tlayer #replace unstable layer with average temp
    
            else: #there is no stratification (min_unstable ==0) and the core is not thermally convecting 
                pass #keep conductive profile
                 
        Tc_new = T_new_core[0] #temperature of core is always taken at centre
        Fad_new = kc*T_new_core[-2]*alpha_c*gc/cpc   

        # Step 5. Find the CMB temperature, how you do this depends on if core and mantle are convecting
        if mantle_conv == True:
            if core_conv == True:
                # eqn 26 and 29 in Dodds 2020 - use bl from previous timestep
                # check if non zero if just switched to convection b.l. at previous timestep = 0
                if dl_old == 0 and dc_old == 0:
                    dc_old = delta_c(Tc_conv_new,Tcmb_old) #find approximate core cmb b.l. thickness
                    dl_old = delta_l(Tm_conv_new,Tcmb_old,Ur_old) #find approximate mantle cmb b.l. thickness
                    factor = (kc*dl_old)/(km*dc_old)
                    
                elif dl_old != 0 and dc_old == 0:
                     dc_old = delta_c(Tc_conv_new,Tcmb_old) #find approximate core cmb b.l. thickness
                     factor = (kc*dl_old)/(km*dc_old)
                     
                elif dl_old == 0 and dc_old != 0:
                     dl_old = delta_l(Tm_conv_new,Tcmb_old,Ur_old) #find approximate mantle cmb b.l. thickness
                     factor = (kc*dl_old)/(km*dc_old)
                     
                else: #both convective b.l. are non zero
                    factor = (kc*dl_old)/(km*dc_old)
                    
                Tcmb_new = (Tm_conv_new + factor*Tc_conv_new)/(1+factor) 
                dc_new = delta_c(Tc_conv_new,Tcmb_new) #find core cmb b.l. thickness
                dl_new = delta_l(Tm_conv_new,Tcmb_new,Ur_old) #find mantle cmb b.l. thickness
                Fcmb_new = -km*(Tm_conv_new-Tcmb_new)/dr#dl_new
            else: #eqn 23 = 24  
                
                Tcmb_new = (T_new_mantle[1]+kc/km*T_new_core[-2])/(1+kc/km)
                dl_new = delta_l(Tm_conv_new,Tcmb_new,Ur_old) #find mantle cmb b.l. thickness
                Fcmb_new = -km*(T_new_mantle[1]-Tcmb_new)/dr # CMB heat flux eqn 23 in Dodds 2020 - until core starts to convect heat transport is modelled as diffusive
        else:
            if core_conv == True:
                #mantle sets CMB temp
                #Tcmb_new = T_new_mantle[0]
                if dc_old == 0:
                    dc_old = delta_c(Tc_conv_new,Tcmb_old) #find approximate core cmb b.l. thickness
                factor = (kc*dr)/(km*dc_old)
                Tcmb_new = (T_new_mantle[1]+factor*T_new_core[-2])/(1+factor)
                dc_new = delta_c(Tc_conv_new,Tcmb_new) #find core cmb b.l. thickness
                Fcmb_new = -km*(T_new_mantle[1]-Tcmb_new)/dr # CMB heat flux eqn 23 in Dodds 2020
            else: # eqn 23 = 24
                Tcmb_new = (T_new_mantle[1]+kc/km*T_new_core[-2])/(1+kc/km)
                Fcmb_new = -km*(T_new_mantle[1]-Tcmb_new)/dr # CMB heat flux eqn 23 in Dodds 2020

        #replace CMB nodes
        T_new_mantle[0] = Tcmb_new
        T_new_core[-1] = Tcmb_new
 
        # Step 6. Replace old values with new ready for next step
        T_old_core = T_new_core
        T_old_mantle = T_new_mantle
        Tprofile_old = np.hstack((T_new_core,T_new_mantle[1:])) #top cell of core and bottom cell of mantle are Tcmb
        Tc_old = Tc_new
        Tc_conv_old = Tc_conv_new
        Tcmb_old = Tcmb_new
        Tm_mid_old = Tm_mid_new
        Tm_conv_old = Tm_conv_new
        Tm_surf_old = Tm_surf_new
        dTdt_mantle_old = dTdt_mantle_new
        f_old = f_new
        Xs_old = Xs_new
        dl_old = dl_new
        dc_old = dc_new
        d0_old = d0_new
        min_unstable_old = min_unstable_new
        Ur_old = Ur_new
        Ra_old = Ra_new
        RaH_old = RaH_new
        RanoH_old = RanoH_new
        RaRob_old = RaRob_new
        Racrit_old = Racrit_new
        Flid_old = Flid_new
        Fs_old = Fs_new
        Fad_old = Fad_new
        Fcmb_old = Fcmb_new
        tsolve_old = tsolve_new
        i+=1 #use for saving
        
        ### add things here
        if i%int(save_interval_t/dt)==0: #save parameters
            save_ind = int(i*dt/save_interval_t)
            Tc[save_ind] = Tc_old
            Tc_conv[save_ind] = Tc_conv_old
            Tcmb[save_ind] = Tcmb_old
            Tm_mid[save_ind] = Tm_mid_old
            Tm_conv[save_ind] = Tm_conv_old
            Tm_surf[save_ind] = Tm_surf_old
            Tprofile[save_ind,:] = Tprofile_old
            f[save_ind] = f_old
            Xs[save_ind] = Xs_old
            dl[save_ind] = dl_old
            dc[save_ind] = dc_old
            d0[save_ind] = d0_old
            min_unstable[save_ind] = min_unstable_old
            Ur[save_ind] = Ur_old
            Ra[save_ind] = Ra_old
            RaH[save_ind] = RaH_old
            RanoH[save_ind] = RanoH_old
            RaRob[save_ind] = RaRob_old
            Racrit[save_ind] = Racrit_old
            Fs[save_ind] = Fs_old
            Flid[save_ind] = Flid_old
            Fad[save_ind] = Fad_old
            Fcmb[save_ind] = Fcmb_old
            Rem_c[save_ind] = Rem_c_new
            Bcomp[save_ind] = Bcomp_new
            tsolve[save_ind] = tsolve_old
        
        if i%int((tstart-tend)/(dt*out_interval))==0: #every 1/out_interval fraction of the run print the time
            print('t={:.2f}Myr'.format(tsolve_new/Myr)) #useful to track progress of simulation
        else: 
            pass
 
        if f_new<=0.001: #stop integration if core is solid i.e.inner liquid core radius < 0.1%, 
            #save values at this point

            if save_ind < len(Tc)-3: #if on final save then don't add one more as not space in array
                save_ind =+ save_ind
                Tc[save_ind] = Tc_old
                Tc_conv[save_ind] = Tc_conv_old
                Tcmb[save_ind] = Tcmb_old
                Tm_mid[save_ind] = Tm_mid_old
                Tm_conv[save_ind] = Tm_conv_old
                Tm_surf[save_ind] = Tm_surf_old
                Tprofile[save_ind,:] = Tprofile_old
                f[save_ind] = f_old
                Xs[save_ind] = Xs_old
                dl[save_ind] = dl_old
                dc[save_ind] = dc_old
                d0[save_ind] = d0_old
                min_unstable[save_ind] = min_unstable_old
                Ur[save_ind] = Ur_old
                Ra[save_ind] = Ra_old
                RaH[save_ind] = RaH_old
                RanoH[save_ind] = RanoH_old
                RaRob[save_ind] = RaRob_old
                Racrit[save_ind] = Racrit_old
                Fs[save_ind] = Fs_old
                Flid[save_ind] = Flid_new
                Fad[save_ind] = Fad_old
                Fcmb[save_ind] = Fcmb_old
                Rem_c[save_ind] = Rem_c_new
                Bcomp[save_ind] = Bcomp_new
                tsolve[save_ind] = tsolve_old
            
            #truncate arrays to only return non-zero values
            Tc = Tc[:save_ind+1]
            Tc_conv = Tc_conv[:save_ind+1]
            Tcmb = Tcmb[:save_ind+1]
            Tm_mid = Tm_mid[:save_ind+1]
            Tm_conv = Tm_conv[:save_ind+1]
            Tm_surf = Tm_surf[:save_ind+1]
            Tprofile = Tprofile[:save_ind+1]
            f=f[:save_ind+1]
            Xs = Xs[:save_ind+1]
            dl = dl[:save_ind+1]
            dc = dc[:save_ind+1]
            d0 = d0[:save_ind+1]
            min_unstable = min_unstable[:save_ind+1]
            Ur = Ur[:save_ind+1]
            Ra = Ra[:save_ind+1]
            RaH = RaH[:save_ind+1]
            RanoH = RanoH[:save_ind+1]
            RaRob = RaRob[:save_ind+1]
            Racrit = Racrit[:save_ind+1]
            Fs = Fs[:save_ind+1]
            Flid = Flid[:save_ind+1]
            Fad = Fad[:save_ind+1]
            Fcmb = Fcmb[:save_ind+1]
            Rem_c = Rem_c[:save_ind+1]
            Bcomp = Bcomp[:save_ind+1]
            tsolve = tsolve[:save_ind+1]
            
            break
        else: 
            pass
              
    #truncate arrays to only return non-zero values in case core hasn't solidified and there is one cell empty
    Tc = Tc[:save_ind+1]
    Tc_conv = Tc_conv[:save_ind+1]
    Tcmb = Tcmb[:save_ind+1]
    Tm_mid = Tm_mid[:save_ind+1]
    Tm_conv = Tm_conv[:save_ind+1]
    Tm_surf = Tm_surf[:save_ind+1]
    Tprofile = Tprofile[:save_ind+1]
    f=f[:save_ind+1]
    Xs = Xs[:save_ind+1]
    dl = dl[:save_ind+1]
    dc = dc[:save_ind+1]
    d0 = d0[:save_ind+1]
    min_unstable = min_unstable[:save_ind+1]
    Ur = Ur[:save_ind+1]
    Ra = Ra[:save_ind+1]
    RaH = RaH[:save_ind+1]
    RanoH = RanoH[:save_ind+1]
    Racrit = Racrit[:save_ind+1]
    RaRob = RaRob[:save_ind+1]
    Fs = Fs[:save_ind+1]
    Flid = Flid[:save_ind+1]
    Fad = Fad[:save_ind+1]
    Fcmb = Fcmb[:save_ind+1]
    Rem_c = Rem_c[:save_ind+1]
    Bcomp = Bcomp[:save_ind+1]
    tsolve = tsolve[:save_ind+1]
            
    return Tc, Tc_conv, Tcmb, Tm_mid, Tm_conv, Tm_surf, Tprofile, f, Xs, dl, dc, d0, min_unstable, Ur, Ra, RaH, RanoH, RaRob, Racrit, Fs, Flid, Fad, Fcmb, Rem_c, Bcomp, tsolve, cond_t