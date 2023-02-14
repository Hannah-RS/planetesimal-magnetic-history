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
import scipy.sparse as sp
import scipy.optimize as sco
from parameters import Ts, Myr, Rac, B, dr, out_interval, km, kc, alpha_m, alpha_c, rc, rhoc, rhom, eta_c, g, gc, cpc, Xs_0, default, kappa, kappa_c, c1, gamma, Xs_eutectic, Acmb, Lc

#import required functions
from Tm_cond import T_cond_calc
from dTmdt_def import dTmdt_calc
from dTcdt_def import dTcdt_calc 
from Rayleigh_def import Rayleigh_calc
from viscosity_def import viscosity
from cmb_bl import delta_l, delta_c
from fe_fes_liquidus import fe_fes_liquidus 
from stratification import volume_average

def thermal_evolution(tstart,tend,dt,T0,f0,sparse_mat_c,sparse_mat_m):
    """
    

    Parameters
    ----------
    tstart : float  
        start time
    tend : float
        end time
    dt : float
        time step
    T0 : array
        initial temperature profile
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
        thickness of mantle CMB boundary layer (0 if mantle not convecting)
    dc: array
        thickness of core CMB boundary layer (0 if core not thermally convecting)
    d0: array
        stagnant lid thickness for convecting mantle
    Ra: array
        Rayleigh number for convecting mantle
    Fs: array
        surface heat flux [W m^-2]
    Fad: array
        adiabatic CMB heat flux [W m^-2]
    Fcmb: array
        CMB heat flux [W m^-2]
    tsolve: array
        time points corresponding to each of the values above
    cond_i: integer
        index in the array when the mantle switched from conduction to convection, nan if didn't switch
         
    """

    #initialise arrays for output
    #p = int((tend-tstart)/(out_interval)) #only output temp profiles every 10 Myr
    m = int((tend-tstart)/dt)
    #ratio = int(m/p) #use for calculating when to save temp profiles
    ratio =100
    p = 1000
    i_save=0
    n_cells = len(T0) #number of cells
    i_core = round(n_cells/2)-1 # index in array of last core cell (-1 as indexing starts at 0)
    cond = 0 #flag for when first switch to conductive regime
    cond_i = 'nan' #set as string and then reset if switches to conduction
    mantle_conv = False #flag for mantle convection
    core_conv = False #flag for core convection
    min_unstable_old = i_core-1 #smallest index of cells in the core that are convectively unstable - as a minimum it is the one below the CMB
    
    #output variables
    Xs = np.ones([m])*Xs_0 #core sulfur fraction
    Ra = np.zeros([m])
    d0 = np.zeros([m])
    dl = np.zeros([m])
    dc = np.zeros([m])
    Tprofile= np.zeros([p,n_cells])
    Tc = np.zeros([m])
    Tc_conv = np.zeros([m])
    Tcmb = np.zeros([m])
    Tm_conv = np.zeros([m])
    Tm_mid = np.zeros([m])
    Tm_surf = np.zeros([m])
    f = np.ones([m])*f0 #set as initial core size by default and only override if core starts to change
    Fs = np.zeros([m])
    Fad = np.zeros([m])
    Fcmb = np.zeros([m])
    tsolve = np.zeros([m])

    
    #Step 0. Calculate time, get two separate temperature arrays
    # the last cell of the core array is the same as the first cell of the mantle array
    tsolve[0] = tstart + dt
    T0_core = T0[:i_core+1] #include last core cell
    T0_mantle = T0[i_core:]
    nmantle_cells = len(T0_mantle)
    ncore_cells = len(T0_core)
    
    ##################    Initial step    #########################
    # Step 1. Calculate conductive profile for mantle
    T_new_mantle = T_cond_calc(tsolve[0],dt,T0_mantle,sparse_mat_m,True)
    Fs[0] = -km*(T_new_mantle[-1]-T_new_mantle[-2])/dr
    
    # Step 2. Is the mantle convecting? Calculate stagnant lid thickness, base thickness and Rayleigh number
    Ra[0], d0[0] = Rayleigh_calc(T0_mantle[1],default) #use temp at base of mantle 
    nlid_cells = round(d0[0]/dr)
    if nlid_cells ==0:
        lid_start = nmantle_cells -2
    else:
        lid_start = nmantle_cells - nlid_cells - 1 #index in temp array where lid starts
    
    
    # in initial step use balance of eqn 23 and 24 to find Tcmb
    Tcmb[0] = (km/kc*T0_mantle[0]+T0_core[0])/(km/kc+1)
    T_new_mantle[0] = Tcmb[0]
    Fcmb[0] = -km*(T0_mantle[0]-Tcmb[0])/dr # CMB heat flux eqn 23 in Dodds 2020
    
    if Ra[0] < Rac:# not convecting

        Tm_conv[0] = 0 # convective mantle temperature is 0 if mantle not convecting
       
    else: #mantle is convecting replace mantle below stagnant lid with isothermal convective profile
        dl[0] = delta_l(T0_mantle[1],T0_mantle[0])
        nbase_cells = round(dl[0]/dr)
        Tm_conv[0] = T0_mantle[nbase_cells] + dTmdt_calc(tsolve[0],Fs[0],Fcmb[0])*dt
        T_new_mantle[nbase_cells:lid_start+1] = Tm_conv[0]
        Fs[0] = -km*(Ts-Tm_conv[0])/d0[0] #replace Fs with convective version
        
    #store values   
    Tm_mid[0] = T_new_mantle[round(nmantle_cells/2)] # temperature at the mid mantle
    Tm_surf[0] = T_new_mantle[-2] #temperature one cell above surface pinned to 200K
                
    # Step 3. Calculate conductive profile for the core
    T_new_core = T_cond_calc(tsolve[0],dt,T0_core,sparse_mat_c,False)
        
    # Step 4. Is the core convecting? 
    # check if heat flux is super adiabatic 
    
    Fad[0] = kc*T0_core[-2]*alpha_c*gc/cpc
    
    if Fcmb[0] > Fad[0]: #super adiabatic, core convects

        dc[0] = delta_c(T0_core[-2],T0_mantle[0]) #second input is CMB temp

        # is the core solidifying?
        Tliquidus = fe_fes_liquidus(Xs_0)
        if np.any(T0_core < Tliquidus) == True: #core solidifies
            dTcdt = dTcdt_calc(tsolve[0], Fcmb[0], T0_core, f0, solidification = True) #save dTcdt seperately as need for f
            dfdt = -B*dTcdt/(T0_core[-2]*f0)
            f[0] = f0 + dfdt*dt
            Xs[0] = (1-(f[0]**3))*Xs_0 #update sulfur content
        
        else: # core not solidifying
            dTcdt = dTcdt_calc(tsolve[0], Fcmb[0], T0_core, f0, solidification = False) 
            f[0]=f0
            Xs[0]=Xs_0
        
        #find number of cells which are solid
        nic_cells = round(f[0]*rc/dr)
        Tc_conv[0] = T0_core[-2] + dTcdt*dt
        T_new_core[:] = Tc_conv[0] #replace everything inside the core with the convective temperature
        
              
    else: # don't have whole core convection, 
        #check if there is thermal stratification
        
        if Tcmb[0] > np.any(T0_core): #there is thermal stratification
            
            if Tcmb[0] > np.all(T0_core):
                    # scenario 1 - just conduction in the core
                    pass # use already calculated condctive profile, don't do anything
            else: #scenario 2 - erosion of stratification, convective layer at top of core
                dc[0] = delta_c(T0_core[-2],T0_mantle[0]) #second input is CMB temp
                b_ind = np.where(Tcmb[0] <= T0_core) #indices of unstable layer
                min_unstable_new = b_ind[0]

                Tc_conv[0] = T0_core[min_unstable_old]
                # is the core solidifying?
                Tliquidus = fe_fes_liquidus(Xs_0)
                if np.any(T0_core) < Tliquidus: #core solidifies
                    raise NotImplementedError('Core solidification erodes thermal stratification - write this code!')
                
                else: # core not solidifying
                    dTcdt = dTcdt_calc(Fcmb[0], T0_core, f0, solidification = False, stratification = [True, min_unstable_old])
                    f[0]=f0
                    Xs[0]=Xs_0
                    
                    Tc_conv[0] = Tc_conv[0]+dTcdt*dt #replace convecting layer from last timestep with new temp - in later steps use i-1 and i
                    T_new_core[min_unstable_old:-1] = Tc_conv[0]
                    T_new_core[-1] = Tcmb[0]
                    
                    #now perform volume average over unstable layer
                    Tlayer = volume_average(T_new_core, b_ind)
                    T_new_core[min_unstable_new:-1] = Tlayer #replace unstable layer with average temp
                    min_unstable_old = min_unstable_new #replace for next step
                    
        
        else: #there is no stratification and the core is not thermally convecting 
            
            """Below here needs rewriting for conductive core solidification/compositional convection"""
            # is the core solidifying?
            Tliquidus = fe_fes_liquidus(Xs_0)
            if T0_core[-2] < Tliquidus: #core solidifies
                dTcdt = dTcdt_calc(Fcmb[0], T0_core, f0, solidification = True) #save dTcdt seperately as need for f
                dfdt = -B*dTcdt/(T0_core[lnb-1]*f0)
                f[0] = f0 + dfdt*dt
                Xs[0] = (1-(f[0])**3)*Xs_0 #update sulfur content
                #Tc_conv[0] = T0_core[lnb-1] + dTcdt*dt
                #T_new_core[:lnb+1] = Tc_conv[0] #replace everything above the solid core
            else: # core not solidifying or convecting keep conductive profile
                pass 
            
            

            
        
    # Step 5. Recalculate the CMB temperature so now it is determined by the end of the step as it will be in subsequent steps
    #use balance of eqn 23 and 24 to find Tcmb
    Tcmb[0] = (km/kc*T_new_mantle[1]+T_new_core[-2])/(km/kc+1)
    T_new_mantle[0] = Tcmb[0]
    T_new_core[-1] = Tcmb[0]
    Fcmb[0] = -km*(T_new_mantle[0]-Tcmb[0])/dr # CMB heat flux eqn 23 in Dodds 2020
    Tc[0] = T_new_core[0] #temperature of core is always taken at centre
    
    # Step 6. Replace old array with new ready for next step
    T_old_core = T_new_core
    T_old_mantle = T_new_mantle
 
    for i in range(1,m):

        #Step 0. Calculate time
        tsolve[i] = tsolve[i-1] + dt
        
        # Step 1. Calculate conductive profile for mantle
        T_new_mantle = T_cond_calc(tsolve[i],dt,T_old_mantle,sparse_mat_m,False)
        Fs[i] = -km*(T_new_mantle[-1]-T_new_mantle[-2])/dr
        
        # Step 2. Is the mantle convecting? Calculate stagnant lid thickness, base thickness and Rayleigh number
        Ra[i], d0[i] = Rayleigh_calc(T_old_mantle[1],default) #use temp at base of mantle 
            
        
        if Ra[i] < Rac or cond==1: #once Rayleigh number subcritical don't want to use that criterion anymore
            mantle_conv = False
            if cond == 0: #check if first time it is conductive i.e. the switch
                cond_i = i
                cond = 1
            else: #mantle already conducting
                pass
            Tm_conv[i] = 0 # convective mantle temperature is 0 if mantle not convecting

        else: #mantle is convecting replace mantle below stagnant lid with isothermal convective profile 
            mantle_conv = True
            nlid_cells = round(d0[i]/dr)
            
            if nlid_cells ==0:
                lid_start = nmantle_cells -2
            else:
                lid_start = nmantle_cells - nlid_cells - 1 #index in temp array where lid starts
            if Tm_conv[i-1]!=0: #mantle already convecting
                Tm_old = Tm_conv[i-1]
            else: #mantle just started convecting
                Tm_old = T_old_mantle[lid_start-1] # take temperature below stagnant lid as mantle temp
            
            
            
            nbase_cells = round(dl[i]/dr)
            Tm_conv[i] = Tm_old + dTmdt_calc(tsolve[i-1],Fs[i-1],Fcmb[i-1])*dt #temperature of convecting region 
            T_new_mantle[nbase_cells:lid_start] = Tm_conv[i]
            Fs[i] = -km*(Ts-Tm_conv[i])/d0[i]
            
            
        #store values
        Tm_mid[i] = T_new_mantle[round(nmantle_cells/2)] # temperature at the mid mantle
        Tm_surf[i] = T_new_mantle[-2] #temperature one cell above surface pinned to 200K
          
        # Step 3. Calculate conductive profile for the core
        T_new_core = T_cond_calc(tsolve[i],dt,T_old_core,sparse_mat_c,False)
   
        # Step 4. Is the core convecting? 
        # check if heat flux is super adiabatic 
        
        if Fcmb[i-1] > Fad[i-1]: #super adiabatic, core convects
            core_conv = True
            
            nbl_cells = round(dc[i-1]/dr)
            bl_start = ncore_cells - nbl_cells - 1 #index in temp array where lid starts

            # is the core solidifying?
            Tliquidus = fe_fes_liquidus(Xs[i-1])
            if np.any(T_old_core < Tliquidus) == True: #core solidifies
                if Xs[i-1]>= Xs_eutectic:
                    dTcdt = 0 # whilst undergoing eutectic solidification there is no temp change
                    dfdt = Fcmb[i-1]*Acmb/(4*np.pi*rc**3*f[i-1]**2*Lc*rhoc)                    
                    f[i] = f[i-1] + dfdt*dt
                    Xs[i] = Xs[i-1] #sulfur concentration unchanged in eutectic solidification
                else:
                    dTcdt = dTcdt_calc(Fcmb[i-1], T_old_core, f[i-1], solidification = True) #save dTcdt seperately as need for f
                    dfdt = -B*dTcdt/(T_old_core[-2]*f[i-1])
                    f[i] = f[i-1] + dfdt*dt
                    Xs[i] = (1-(f[i]**3))*Xs_0 #update sulfur content
            
            else: # core not solidifying
                dTcdt = dTcdt_calc(Fcmb[i-1], T_old_core, f[i-1], solidification = False)
                f[i]=f[i-1] #keep the current core size
                Xs[i]=Xs[i-1]
            
            #find new convective temperature
            Tc_conv[i] = T_old_core[bl_start-1] + dTcdt*dt 
            #print(f"Tc_conv[i] is {Tc_conv[i]:.2f}")
            T_new_core[:] = Tc_conv[i] #replace everything with the convective temperature
            
            
            
        else: # don't have whole core convection, 
            #check if there is thermal stratification
            core_conv = False 
            if Tcmb[i-1] > np.any(T_old_core): #there is thermal stratification
                
                if Tcmb[i-1] > np.all(T_old_core):
                        # scenario 1 - just conduction in the core
                        # use already calculated condctive profile and keep core in current state
                        f[i] = f[i-1]
                        Xs[i] = Xs[i-1]
                else: #scenario 2 - erosion of stratification, convective layer at top of core
                    core_conv = True
                    b_ind = np.where(Tcmb[i-1] <= T_old_core) #indices of unstable layer
                    
                    min_unstable_new = b_ind[0]

                    
                    
                    # is the core solidifying?
                    Tliquidus = fe_fes_liquidus(Xs[i-1])
                    if np.any(T_old_core) < Tliquidus: #core solidifies
                        raise NotImplementedError('Core solidification erodes thermal stratification - write this code!')
                    
                    else: # core not solidifying
                        dTcdt = dTcdt_calc(Fcmb[i-1], T_old_core, f[i-1], solidification = False, stratification = [True, min_unstable_old])
                        f[i]=f[i-1]
                        Xs[i]=Xs[i-1]
                        
                        Tc_conv[i] = T_old_core[min_unstable_old]+dTcdt*dt #replace convecting layer from last timestep with new temp - in later steps use i-1 and i
                        T_new_core[min_unstable_old:-1] = Tc_conv[i]
                                                 
                        #now perform volume average over unstable layer
                        Tlayer = volume_average(T_new_core, b_ind)
                        T_new_core[min_unstable_new:-1] = Tlayer #replace unstable layer with average temp
                        
                        min_unstable_old = min_unstable_new #replace for next step
                           
            
            else: #there is no stratification and the core is not thermally convecting 
                
                """Below here needs rewriting for conductive core solidification/compositional convection"""
                # is the core solidifying?
                Tliquidus = fe_fes_liquidus(Xs[i-1])
                if T_old_core[-2] < Tliquidus: #core solidifies
                    if Xs[i-1]>= Xs_eutectic: #eutectic solidification
                        dTcdt = 0 # whilst undergoing eutectic solidification there is no temp change
                        dfdt = Fcmb[i-1]*Acmb/(4*np.pi*rc**3*f[i-1]**2*Lc*rhoc)                    
                        f[i] = f[i-1] + dfdt*dt
                        Xs[i] = Xs[i-1] #sulfur concentration unchanged in eutectic solidification
                    else:
                        dTcdt = dTcdt_calc(Fcmb[i-1], T_old_core, f[i-1], solidification = True) #save dTcdt seperately as need for f
                        dfdt = -B*dTcdt/(T_old_core[-2]*f[i-1])
                        f[i] = f[i-1] + dfdt*dt
                        Xs[i] = (1-(f[i])**3)*Xs_0 #update sulfur content
                    
                else: # core not solidifying or convecting keep conductive profile
                    f[i] = f[i-1]
                    Xs[i] = Xs[i-1]
                 
        Tc[i] = T_new_core[0] #temperature of core is always taken at centre
        Fad[i] = kc*T_new_core[-2]*alpha_c*gc/cpc   

        # Step 5. Find the CMB temperature, how you do this depends on if core and mantle are convecting
        if mantle_conv == True:
            if core_conv == True:
                # eqn 26 and 29 in Dodds 2020 - use bl from previous timestep
                # check if non zero if just switched to convection b.l. at previous timestep = 0
                if dl[i-1] == 0 and dc[i-1] == 0:
                    dc[i-1] = delta_c(Tc_conv[i],Tcmb[i-1]) #find approximate core cmb b.l. thickness
                    dl[i-1] = delta_l(Tm_conv[i],Tcmb[i-1]) #find approximate mantle cmb b.l. thickness
                    factor = (kc*dl[i-1])/(km*dc[i-1])
                    
                elif dl[i-1] != 0 and dc[i-1] == 0:
                     dc[i-1] = delta_c(Tc_conv[i],Tcmb[i-1]) #find approximate core cmb b.l. thickness
                     factor = (kc*dl[i-1])/(km*dc[i-1])
                     
                elif dl[i-1] == 0 and dc[i-1] != 0:
                     dl[i-1] = delta_l(Tm_conv[i],Tcmb[i-1]) #find approximate mantle cmb b.l. thickness
                     factor = (kc*dl[i-1])/(km*dc[i-1])
                     
                else: #both convective b.l. are non zero
                    factor = (kc*dl[i-1])/(km*dc[i-1])
                    
                Tcmb[i] = (Tm_conv[i] + factor*Tc_conv[i])/(1+factor) 
                dc[i] = delta_c(Tc_conv[i],Tcmb[i]) #find core cmb b.l. thickness
                dl[i] = delta_l(Tm_conv[i],Tcmb[i]) #find mantle cmb b.l. thickness
                Fcmb[i] = -km*(Tm_conv[i]-Tcmb[i])/dr#dl[i]
            else: #eqn 23 = 24  
                
                Tcmb[i] = (T_new_mantle[1]+kc/km*T_new_core[-2])/(1+kc/km)
                dl[i] = delta_l(Tm_conv[i],Tcmb[i]) #find mantle cmb b.l. thickness
                Fcmb[i] = -km*(T_new_mantle[1]-Tcmb[i])/dr # CMB heat flux eqn 23 in Dodds 2020 - until core starts to convect heat transport is modelled as diffusive
        else:
            if core_conv == True:
                #mantle sets CMB temp
                #Tcmb[i] = T_new_mantle[0]
                if dc[i-1] == 0:
                    dc[i-1] = delta_c(Tc_conv[i],Tcmb[i-1]) #find approximate core cmb b.l. thickness
                factor = (kc*dr)/(km*dc[i-1])
                Tcmb[i] = (T_new_mantle[1]+factor*T_new_core[-2])/(1+factor)
                dc[i] = delta_c(Tc_conv[i],Tcmb[i]) #find core cmb b.l. thickness
                Fcmb[i] = -km*(T_new_mantle[1]-Tcmb[i])/dr # CMB heat flux eqn 23 in Dodds 2020
            else: # eqn 23 = 24
                Tcmb[i] = (T_new_mantle[1]+kc/km*T_new_core[-2])/(1+kc/km)
                Fcmb[i] = -km*(T_new_mantle[1]-Tcmb[i])/dr # CMB heat flux eqn 23 in Dodds 2020

        #replace CMB nodes
        T_new_mantle[0] = Tcmb[i]
        T_new_core[-1] = Tcmb[i]
 
        # Step 6. Replace old array with new ready for next step
        T_old_core = T_new_core
        T_old_mantle = T_new_mantle

        #write Tprofile to file if appropriate
        if i%ratio== 0: #if multiple of 10Myr then save the profile - will need to check if this works as tsolve might not be integer multiples of 10 ever
            print('t={:.2f}Myr'.format(tsolve[i]/Myr)) #useful to track progress of simulation
            if i_save >= p: # array is full, pass so don't throw an error
                pass
            else:
                T_old = np.hstack([T_old_core,T_old_mantle[1:]])
                Tprofile[i_save,:] = T_old
                i_save = i_save + 1 #increment so saves in next space
        else: 
            pass
 
        if f[i]>=1: #stop integration if core is solid, truncate arrays to only return non-zero values
            
            
            Tc = Tc[:i+1]
            Tc_conv = Tc_conv[:i+1]
            Tcmb = Tcmb[:i+1]
            Tm_mid = Tm_mid[:i+1]
            Tm_conv = Tm_conv[:i+1]
            Tm_surf = Tm_surf[:i+1]
            Tprofile = Tprofile[:i_save+1]
            f=f[:i+1]
            Xs = Xs[:i+1]
            dl = dl[:i+1]
            dc = dc[:i+1]
            d0 = d0[:i+1]
            Ra = Ra[:i+1] 
            Fs = Fs[:i+1]
            Fad = Fad[:i+1]
            Fcmb = Fcmb[:i+1]
            tsolve = tsolve[:i+1]
            
            break
        else: 
            pass
              
  
    return Tc, Tc_conv, Tcmb, Tm_mid, Tm_conv, Tm_surf, Tprofile, f, Xs, dl, dc, d0, Ra,  Fs, Fad, Fcmb, tsolve, cond_i