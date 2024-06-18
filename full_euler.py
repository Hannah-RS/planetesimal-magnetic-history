#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#import modules
import numpy as np
from parameters import Ts, Myr, dr, out_interval, save_interval_t, km, kc, alpha_c,\
    r, rc, rhoc, gc, Vm, rhom, As, cpc, Xs_0, default, Xs_eutectic, Acmb, Lc, \
        Pc, automated, conv_tol, n_cells, temp_tol, rcr, nccells, nmcells

#import required functions
from temp_cond import Tm_cond_calc, Tc_cond_calc
from dTmdt_def import dTmdt_calc
from dTcdt_def import dTcdt_calc, dTcdt_calc_solid 
from rayleigh_def import rayleigh_calc, rayleigh_crit
from cmb_bl import delta_l, delta_c
from q_funcs import qr
from fe_fes_liquidus import fe_fes_liquidus_bw 
from stratification import volume_average
from heating import al_heating

def thermal_evolution(tstart,tend,dt,T0,f0,sparse_mat_c,sparse_mat_m):
    """
    Script for solving thermal evolution for the entire evolution of the body after differentiation.
    Workflow:
        1. Obtain conductive profile for whole body
        2. Calculate stagnant lid thickness and Rayleigh number
        3. If d0 + dl < R - rc replace rc < r < R-d0 with convective profile for mantle
        4. Replace r < rc with isothermal core profile if the core is convecting
        5. Calculate change in inner core radius, magnetic Reynolds number and field strength
    Calculation ends when the core has solidified or specified end time has been exceeded.

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
        core sulfur content [wt %]
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
    Racrit: array
        critical Rayleigh number for convecting mantle
    eta : array
        mantle viscosity [Pas]
    Flid: array
        heat flux across base of stagnant lid [W m^-2]
    Fs: array
        surface heat flux [W m^-2]
    Fad: array
        adiabatic CMB heat flux [W m^-2]
    Fcmb: array
        CMB heat flux [W m^-2]
    Rem : array
        magnetic Reynolds number
    B : float
         magnetic field strength [T]
    buoyr : ndarray
         compositional buoyr[0,:] and thermal buoyr[1,:] buoyancy fluxes [kg/s]
    qcore : np.ndarray
            heat sources in the core, qr=radiogenic, qs=secular cooling
            ql= latent heat release, qg=gpe release [qr, qs, ql, qg] [W]
    tsolve: array
        time points corresponding to each of the values above [s]
    fcond_t : float
        time for cessation of convection [Myr]
         
    """

    #initialise arrays for output
    m = round((tend-tstart)/save_interval_t)+1 #add one so always enough space
    i_core = round((n_cells-3)*(rcr))+1 # index in array of last core cell 
    mantle_conv = False #flag for mantle convection
    core_conv = False #flag for core convection
    fcond_t = np.nan #end time of convection - default nan
    conv_off = False #whether lid has thickened sufficiently to switch off convection
    #output variables
    Xs = np.ones([m])*Xs_0 #core sulfur fraction
    Ra = np.zeros([m])
    Racrit = np.zeros([m])
    RaH = np.zeros([m])
    RanoH = np.zeros([m])
    eta = np.zeros([m])
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
    Rem = np.zeros([m])
    B = np.zeros([m])
    buoyr = np.zeros([2,m])
    qcore = np.zeros([4,m])
    tsolve = np.zeros([m])
        
    #Step 0. Calculate time, get two separate temperature arrays
    # assume mantle and core are isothermal at point of differentiation - valid due to strong heating
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
    buoyr_new = [0,0]
    
    ##################    Initial step    #########################
    # Step 1. Calculate conductive profile for mantle
    T_new_mantle = Tm_cond_calc(tsolve_new,dt,T0_mantle,sparse_mat_m)
    dTdt_mantle_new = (T_new_mantle-T0_mantle)[1] #default use temp change at base above CMB
    Tm_conv_new = 0 #default convective temp
    Fs_new = -km*(T_new_mantle[-1]-T_new_mantle[-2])/dr
    Flid_new = Fs_new #by default lid flux is same as surface
    h = al_heating(tsolve_new)
    Ur_new = rhom*Vm*h/abs(Fs_new*As) #calculate Urey ratio
    
    # Step 2. Is the mantle convecting? Calculate stagnant lid thickness, base thickness and Rayleigh number
    Ra_new, d0_new, RaH_new, RanoH_new, eta_new = rayleigh_calc(tsolve_new,T0_mantle[1],Ur_new,default) #use temp at base of mantle 
    Racrit_new = rayleigh_crit(T0_mantle[1])   
    nlid_cells = round(d0_new/dr)
    if nlid_cells ==0:
        lid_start = nmantle_cells -2
    else:
        lid_start = nmantle_cells - nlid_cells - 1  #index in temp array where stagnant lid starts
    
    # in initial step use balance of eqn 23 and 24 to find first estimate of Tcmb
    Tcmb_new = (km/kc*T0_mantle[1]+T0_core[-2])/(km/kc+1)
    T_new_mantle[0] = Tcmb_new
    Fcmb_new = -km*(T0_mantle[1]-Tcmb_new)/dr # CMB heat flux eqn 23 in Dodds 2020
    
    #CMB and mantle isothermal at differentiation so assume no bottom b.l.
    stratification_old = True 
    lid_thickness = d0_new 
    dl_new = 0 #no cmb b.l.
    nbase_cells = 0 #no cmb b.l.
    
    if np.round((lid_thickness/(r-rc)),1)<1: #mantle convecting, 
 
        mantle_conv = True               
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
    Rem_new = 0 # by default
    B_new = 0
    
    Fad_new = kc*T0_core[-2]*alpha_c*gc/cpc
    # Step 4. Is the core solidifying? 
    # is the core solidifying?
    Tliquidus = fe_fes_liquidus_bw(Xs_0,Pc)
    if T0_core[-2] < Tliquidus: #core solidifies outside in - convecting so isothermal beneath CMB
        core_conv = False
        stratification_old = False #no stratification when solidifying
        if Xs_0>= Xs_eutectic:
            dTcdt = 0 # whilst undergoing eutectic solidification there is no temp change
            dfdt = -(Fcmb_new*Acmb-qr(tsolve_new))/(4*np.pi*rc**3*f0**2*Lc*rhoc)                    
            f_new = f0 + dfdt*dt
            Xs_new = Xs_0 #sulfur concentration unchanged in eutectic solidification
            #find new convective temperature
            Tc_conv_new = T0_core[0] 
            T_new_core[:] = Tc_conv_new #replace everything with the convective temperature
            Rem_new = 0
            B_new = 0
        else:             
            dTcdt, f_new, Rem_new, B_new, buoyr_new[0], buoyr_new[1], qcore_new = dTcdt_calc_solid(tsolve_new,Fcmb_new, T0_core, f0, Xs_0, dt) 
            #find new convective temperature
            Tc_conv_new = T0_core[0] + dTcdt*dt 
            T_new_core[:] = Tc_conv_new #replace everything with the convective temperature
            Xs_new = Xs_0/(f_new**3) #update sulfur content
            
    # Is the core convecting  without solidification?
    elif Fcmb_new > 0: #Fcmb > 0, core convects, don't check for stratification as haven't recalculated Tcmb yet
        core_conv = True
        nbl_cells = round(dc_new/dr)
        bl_start = ncore_cells - nbl_cells - 1 #index in temp array where lid starts
        dTcdt, Rem_new, B_new, buoyr_new[0], buoyr_new[1], qcore_new = dTcdt_calc(tsolve_new,Fcmb_new, T0_core, f0,Xs_0)
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
    Fcmb_new = -km*(T_new_mantle[1]-Tcmb_new)/dr # CMB heat flux eqn 23 in Dodds 2020
    Tc_new = T_new_core[0] #temperature of core is always taken at centre
    
    #Check for initial stratification
    if np.any(T_new_core[:-1]-Tcmb_new < temp_tol):
        stratification_new = True
        if np.all(T_new_core[:-1]-Tcmb_new < temp_tol):
            min_unstable_new = i_core -1
        else:
            b_ind = np.where(T_new_core[:-1]-Tcmb_new >= temp_tol)[0] #indices of unstable layer as array
            min_unstable_new = b_ind[0]
            dl_new = 0 #reset dl as 0 when core is stratified
    else:
        stratifcation_new = False
        min_unstable_new = 0
        
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
    stratification_old = stratification_new
    Ur_old = Ur_new
    Ra_old = Ra_new
    RaH_old = RaH_new
    RanoH_old = RanoH_new
    Racrit_old = Racrit_new
    Flid_old = Flid_new
    Fs_old = Fs_new
    Fad_old = Fad_new
    Fcmb_old = Fcmb_new
    tsolve_old = tsolve_new
    
    while tsolve_new < tend:

        #Step -1. Check for unphysicality
        if f_old > f0:
            raise ValueError("Inner core exceeds core size")
        if np.any(T_old_mantle > 2000):
            raise ValueError("Mantle overheating! Temperature exceeds 2000K.")
        
        #Step 0. Calculate time
        tsolve_new = tsolve_old + dt
        
        # Step 1. Calculate conductive profile for mantle
        T_new_mantle = Tm_cond_calc(tsolve_new,dt,T_old_mantle,sparse_mat_m)
        Fs_new = -km*(T_new_mantle[-1]-T_new_mantle[-2])/dr
        Flid_new = Fs_new #by default lid flux is same as surface
        dTdt_mantle_new = (T_new_mantle-T_old_mantle)[1] #default use temp change at base just above CMB
        
        #Step 2. Is the mantle convecting?
        if conv_off == False: #check convection hasn't switched off
            if stratification_old == False:
                lid_thickness = d0_old + dl_old #only include cmb b.l. when core not stratified
                nbase_cells = round(dl_old/dr)
            else:
                lid_thickness = d0_old 
                dl_new = 0 #no cmb b.l.
                nbase_cells = 0 #no cmb b.l.
        
        if conv_off == False: #convection still could be on
            if np.round((lid_thickness/(r-rc)),1)<1: #mantle convecting, use temp below stagnant lid
                #calculate number of cells in stagnant lid 
                nlid_cells = round(d0_old/dr)
                if nlid_cells ==0:
                    lid_start = nmantle_cells -2
                else:
                    lid_start = nmantle_cells - nlid_cells - 1  #index in temp array where lid starts
                if Tm_conv_old ==0: #if was conducting use base temp
                    Tm_conv_old = T_old_mantle[1] #use mantle base temp
                Ra_new, d0_new, RaH_new, RanoH_new, eta_new = rayleigh_calc(tsolve_new,Tm_conv_old,Ur_old,default) #Use the temperature just below the lid 
                Racrit_new = rayleigh_crit(Tm_conv_old)
                #calculate temp change
                mantle_conv = True               
                dTdt_mantle_conv = dTmdt_calc(tsolve_old,Tm_conv_old,d0_old,Flid_old,Fcmb_old) #convective dTdt - use temp below lid
                dTdt_mantle_new = dTdt_mantle_conv
                Tm_conv_new = Tm_conv_old + dTdt_mantle_new*dt #temperature of convecting region 
                T_new_mantle[nbase_cells+1:lid_start+1] = Tm_conv_new
                T_new_mantle[-1] = Ts #pin surface to 200K in case d0 < 1 cell thick
                if d0_new < dr:
                    Flid_new = -km*(Ts-Tm_conv_new)/d0_new #flux from convecting region to stagnant lid
                    Fs_new = Flid_new #lid determines flux out of surface
                else:
                    Flid_new = -km*(T_new_mantle[lid_start+1]-T_new_mantle[lid_start])/dr
                
            else: #mantle not convecting, use basal mantle temp
                Ra_new, d0_new, RaH_new, RanoH_new, eta_new = rayleigh_calc(tsolve_new,T_old_mantle[1],Ur_old,default) #use temp at base of mantle 
                Racrit_new = rayleigh_crit(T_old_mantle[1])            
                
                if tsolve_new/Myr > 5: #turn off convection if mantle no longer heating up
                    conv_off = True
                    fcond_t = tsolve_new/Myr
                    print(f'Lid is too thick convection turned off at {tsolve_new/Myr:.2f} Myr')
                mantle_conv = False
                Tm_conv_new = 0

        else: #lids too thick convection turned off
            mantle_conv = False
            Tm_conv_new = 0
        
        # Calculate Urey ratio with correct heat flux
        h = al_heating(tsolve_new)
        Ur_new = rhom*Vm*h/abs(Fs_new*As)
        
        #store values
        Tm_mid_new = T_new_mantle[round(nmantle_cells/2)] # temperature at the mid mantle
        Tm_surf_new = T_new_mantle[-2] #temperature one cell above surface pinned to 200K
          
        # Step 3. Calculate conductive profile for the core
        T_new_core = Tc_cond_calc(tsolve_new,dt,T_old_core,sparse_mat_c,True)
        Tc_conv_new = 0 #by default, overwrite if core convects
        Rem_new = 0 #by default assume no compositional convection
        B_new = 0
        buoyr_new = [0, 0]
        qcore_new = np.zeros([4]) #different heat sources not relevant when conducting
        f_new = f_old #by default overwrite if solidifies
        Xs_new = Xs_old
        min_unstable_new = min_unstable_old #continuity of mixed layer thickness by default
        
        # Step 4. Is the core solidifying? 
        # is the core solidifying?
        Tliquidus = round(fe_fes_liquidus_bw(Xs_old,Pc)) #round as Buono & Walker data to nearest integer
        if T_old_core[-2] < Tliquidus: #core solidifies from outside in- convecting so isothermal beneath CMB
            core_conv = False #core convects but b.l. thickness set by rho not T so use conductive Fcmb
            stratification_new = False #thermal stratification is removed by solidification
            if Xs_old>= Xs_eutectic:
                dTcdt = 0 # whilst undergoing eutectic solidification there is no temp change
                dfdt = -(Fcmb_old*Acmb-qr(tsolve_new))/(4*np.pi*rc**3*f_old**2*Lc*rhoc)                    
                f_new = f_old + dfdt*dt
                Xs_new = Xs_old #sulfur concentration unchanged in eutectic solidification
                #find new convective temperature
                Tc_conv_new = T_old_core[0] 
                T_new_core[:] = Tc_conv_new #replace everything with the convective temperature
                
            else:              
                dTcdt, f_new, Rem_new, B_new, buoyr_new[0], buoyr_new[1], qcore_new = dTcdt_calc_solid(tsolve_new,Fcmb_old, T_old_core, f_old, Xs_old, dt) 
                #find new convective temperature
                Tc_conv_new = T_old_core[0] + dTcdt*dt 
                T_new_core[:] = Tc_conv_new #replace everything with the convective temperature
                Xs_new = Xs_0/(f_new**3) #update sulfur content
         
        # Is the core convecting  without solidification?
        elif (Fcmb_old > 0) and (stratification_old==False): #super adiabatic and no stratification, core convects
            core_conv = True
            stratification_new = False
            min_unstable_new = 0
            nbl_cells = round(dc_old/dr)
            bl_start = ncore_cells - nbl_cells - 1 #index in temp array where lid starts

            dTcdt, Rem_new, B_new, buoyr_new[0], buoyr_new[1], qcore_new = dTcdt_calc(tsolve_new,Fcmb_old, T_old_core, f_old,Xs_old)
            
            #find new convective temperature
            Tc_conv_new = T_old_core[bl_start-1] + dTcdt*dt 
            
            T_new_core[:bl_start] = Tc_conv_new #replace everything with the convective temperature up to b.l
            
        else: # don't have whole core convection, 
            #check if there is thermal stratification
            core_conv = False 
            
            if stratification_old == True:  #there is thermal stratification 
                
                if np.all(T_old_core[:-1]-Tcmb_old < temp_tol): #stratification stable
                    stratification_new = True
                    min_unstable_new = i_core - 1
                        # scenario 1 - just conduction in the core
                        # use already calculated conductive profile and keep core in current state

                else: #scenario 2 - erosion of stratification, convective layer at top of core
                    core_conv = True
                    b_ind = np.where(T_old_core[:-1]-Tcmb_old >= temp_tol)[0] #indices of unstable layer as array
                    min_unstable_new = b_ind[0]
                    dTcdt, Rem_new, B_new, buoyr_new[0], buoyr_new[1], qcore_new = dTcdt_calc(tsolve_new,Fcmb_old, T_old_core, f_old, Xs_old, stratification = [True, min_unstable_old])
                    Tc_conv_new = T_old_core[min_unstable_old]+dTcdt*dt #replace convecting layer from last timestep with new temp - in later steps use i-1 and i
                    T_new_core[min_unstable_old:-1] = Tc_conv_new
                    #now perform volume average over unstable layer
                    Tlayer = volume_average(T_new_core, b_ind,dr)
                    T_new_core[min_unstable_new:-1] = Tlayer #replace unstable layer with average temp
                    if min_unstable_new ==0:
                        stratification_new = False #the whole core is now convecting in the next step
                    else:
                        stratification_new = True
            else: #there was no stratification in previous step, the core is not convecting
                if np.all(T_old_core[:-1]-Tcmb_old < temp_tol): #is there now stratification?
                    stratification_new = True #add stratification for next step
                    min_unstable_new = i_core - 1
                else: #core is hotter than the mantle and the core is not thermally convecting 
                    stratification_new = False #keep conductive profile
                    min_unstable_new = 0
                 
        Tc_new = T_new_core[0] #temperature of core is always taken at centre
        Fad_new = kc*T_new_core[-2]*alpha_c*gc/cpc   

        # Step 5. Find the CMB temperature, how you do this depends on if core and mantle are convecting
        if stratification_new == False:
            if mantle_conv == True:
                if core_conv == True:
                    # eqn 26 and 29 in Dodds 2020 - use bl from previous timestep
                    # check if non zero if just switched to convection b.l. at previous timestep = 0
                    if dl_old == 0 and dc_old == 0:
                        dc_old = delta_c(Tc_conv_new,Tcmb_old) #find approximate core cmb b.l. thickness
                        dl_old = delta_l(Tm_conv_new,Tcmb_old) #find approximate mantle cmb b.l. thickness
                        factor = (kc*dl_old)/(km*dc_old)
                        
                    elif dl_old != 0 and dc_old == 0:
                         dc_old = delta_c(Tc_conv_new,Tcmb_old) #find approximate core cmb b.l. thickness
                         factor = (kc*dl_old)/(km*dc_old)
                         
                    elif dl_old == 0 and dc_old != 0:
                         dl_old = delta_l(Tm_conv_new,Tcmb_old) #find approximate mantle cmb b.l. thickness
                         factor = (kc*dl_old)/(km*dc_old)
                         
                    else: #both convective b.l. are non zero
                        factor = (kc*dl_old)/(km*dc_old)
                        
                    Tcmb_new = (Tm_conv_new + factor*Tc_conv_new)/(1+factor) 
                    dc_new = delta_c(Tc_conv_new,Tcmb_new) #find core cmb b.l. thickness
                    dl_new = delta_l(Tm_conv_new,Tcmb_new) #find mantle cmb b.l. thickness
                    Fcmb_new = -km*(Tm_conv_new-Tcmb_new)/dl_new
                    
                else: #balance mantle b.l. flux and diffusive flux from core  
                    factor = (kc*dl_old)/(km*dr)
                    Tcmb_new = (Tm_conv_new + factor*T_new_core[-2])/(1+factor)
                    dl_new = delta_l(Tm_conv_new,Tcmb_new) #find mantle cmb b.l. thickness
                    Fcmb_new = -km*(T_new_mantle[nbase_cells+1]-Tcmb_new)/dl_new # CMB heat flux eqn 29 in Dodds 2020 - until core starts to convect heat transport is modelled as diffusive
        
            else:
                if core_conv == True:
                    if dc_old == 0: #if core wasn't convecting in previous step
                         dc_old = delta_c(Tc_conv_new,Tcmb_old) #find approximate core cmb b.l. thickness
                    factor = (kc*dr)/(km*dc_old)
                    Tcmb_new = (T_new_mantle[1]+ factor*Tc_conv_new)/(1+factor)
                    dc_new = delta_c(Tc_conv_new,Tcmb_new) #find core cmb b.l. thickness
                    Fcmb_new = -km*(T_new_mantle[1]-Tcmb_new)/dr # CMB heat flux eqn 23 in Dodds 2020
                else:                
                    Tcmb_new = (T_new_mantle[1]+kc/km*T_new_core[-2])/(1+kc/km)
                    Fcmb_new = -km*(T_new_mantle[1]-Tcmb_new)/dr # CMB heat flux eqn 23 in Dodds 2020  
        else:                
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
        stratification_old = stratification_new
        Ur_old = Ur_new
        Ra_old = Ra_new
        RaH_old = RaH_new
        RanoH_old = RanoH_new
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
            Racrit[save_ind] = Racrit_old
            eta[save_ind] = eta_new
            Fs[save_ind] = Fs_old
            Flid[save_ind] = Flid_old
            Fad[save_ind] = Fad_old
            Fcmb[save_ind] = Fcmb_old
            Rem[save_ind] = Rem_new
            B[save_ind] = B_new
            buoyr[:,save_ind] = buoyr_new
            qcore[:,save_ind] = qcore_new
            tsolve[save_ind] = tsolve_old
        
        if i%int((tstart-tend)/(dt*out_interval))==0: #every 1/out_interval fraction of the run print the time
            if automated == False: #don't print on automated runs
                print('t={:.2f}Myr'.format(tsolve_new/Myr)) #useful to track progress of simulation
 
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
                Racrit[save_ind] = Racrit_old
                eta[save_ind] = eta_new
                Fs[save_ind] = Fs_old
                Flid[save_ind] = Flid_new
                Fad[save_ind] = Fad_old
                Fcmb[save_ind] = Fcmb_old
                Rem[save_ind] = Rem_new
                B[save_ind] = B_new
                buoyr[:,save_ind] = buoyr_new
                qcore[:,save_ind] = qcore_new
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
            Racrit = Racrit[:save_ind+1]
            eta = eta[:save_ind+1]
            Fs = Fs[:save_ind+1]
            Flid = Flid[:save_ind+1]
            Fad = Fad[:save_ind+1]
            Fcmb = Fcmb[:save_ind+1]
            Rem = Rem[:save_ind+1]
            B = B[:save_ind+1]
            buoyr = buoyr[:,:save_ind+1]
            qcore = qcore[:,:save_ind+1]
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
    eta = eta[:save_ind+1]
    Fs = Fs[:save_ind+1]
    Flid = Flid[:save_ind+1]
    Fad = Fad[:save_ind+1]
    Fcmb = Fcmb[:save_ind+1]
    Rem = Rem[:save_ind+1]
    B = B[:save_ind+1]
    buoyr = buoyr[:,:save_ind+1]
    qcore = qcore[:,:save_ind+1]
    tsolve = tsolve[:save_ind+1]
            
    return Tc, Tc_conv, Tcmb, Tm_mid, Tm_conv, Tm_surf, Tprofile, f, Xs, dl, dc, d0, min_unstable, Ur, Ra, RaH, RanoH, Racrit, eta, Fs, Flid, Fad, Fcmb, Rem, B, buoyr, qcore, tsolve, fcond_t