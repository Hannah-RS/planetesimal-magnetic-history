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
def thermal_evolution(tstart,tend,dt,T0,f0):
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

    Returns
    -------
    Ra: array
        Rayleigh number for convecting mantle
    d0: stagnant lid thickness for convecting mantle
    Tprofile: array
        temperature profile as a function of time, output every 10Myr
    Tc: array
        core temperature
    Tm_base: array
        mantle temp one node above CMB
    Tm_surf: array
        mantle temp one node below surface
    f: array
        fractional inner core radius as function of time
    tsolve: array
        time points corresponding to each of the values above
    cond_i: integer
        index in the array when the mantle switched from conduction to convection, nan if didn't switch
         
    """
    import numpy as np
    from parameters import Myr, Rac, B, Tsolidus, dr
    
    #initialise arrays for output
    out_interval = 0.01*Myr
    p = int((tend-tstart)/(out_interval)) #only output temp profiles every 10 Myr
    m = int((tend-tstart)/dt)
    ratio = int(m/p) #use for calculating when to save temp profiles
    i_save=0
    n_cells = len(T0) #number of cells
    i_core = int(n_cells/2) # index in array of last core cell
    cond = 0 #flag for when first switch to conductive regime
    
    #output variables
    Ra = np.zeros([m])
    d0 = np.zeros([m])
    Tprofile= np.zeros([p,n_cells])
    Tc = np.zeros([m])
    Tm_base = np.zeros([m])
    Tm_surf = np.zeros([m])
    f = np.zeros([m])
    tsolve = np.zeros([m])
    
    #import required functions
    from Tm_cond import Tm_cond_calc
    from dTmdt_def import dTmdt_calc
    from dTcdt_def import dTcdt_calc #convective CMB heat flux
    from dTcdt_def2 import dTcdt_calc2  #conductive CMB heat flux
    
    #Step 0. Calculate time
    tsolve[0] = tstart + dt
    
    # Step 1. Calculate conductive profile for whole body
    T_new = Tm_cond_calc(dt,T0)
    
    # Step 2. Calculate stagnant lid thickness and Rayleigh number
    from Rayleigh_def import Rayleigh_calc
    Ra[0], d0[0] = Rayleigh_calc(T0[i_core+1]) #use temp at base of mantle 
    nlid_cells = int(d0[0]/dr)
    lid_start = n_cells - nlid_cells - 1 #index in temp array where lid starts

    if Ra[0] < Rac:
        Tm_base[0] = T_new[i_core+1]
        Tm_surf[0] = T_new[-2]
        # don't replace any of the mantle temperature, just replace the core
        # Step 4. Replace core with isothermal profile
        dTcdt = dTcdt_calc2(tsolve[0],T0[i_core+1],T0[i_core],f0) #save dTcdt seperately as need for f
        Tc[0] = T0[i_core]+dTcdt*dt
        T_new[:i_core+1] = Tc[0] #replace core including core cell
    else: 
        # Step 3. Replace mantle below stagnant lid with isothermal convective profile
        Tm = T0[i_core+1] + dTmdt_calc(tsolve[0],T0[i_core+1],T0[i_core])*dt
        T_new[:lid_start] = Tm
        
        Tm_base[0] = T_new[i_core+1]
        Tm_surf[0] = T_new[-2]
        
        # Step 4. Replace core with isothermal profile
        dTcdt = dTcdt_calc(tsolve[0],T0[i_core+1],T0[i_core],f0,d0[0])
        Tc[0] = T0[i_core]+ dTcdt*dt
        T_new[:i_core+1] = Tc[0] #replace core including core cell
    
    # Step 5. Calculate f
    if Tc[0] > Tsolidus:
        f[0]=f0  #core not solidifying - based on 3-9 wt% S and phase diagram in Scheinberg 2016 
    else:  
    
        dfdt = -B*dTcdt/(T0[i_core]*f0)
        f[0] = f0 + dfdt*dt
    
    T_old = T_new 

    
    for i in range(1,m):
      
        #Step 0. Calculate time
        tsolve[i] = tsolve[i-1] + dt
        
        # Step 1. Calculate conductive profile for whole body
        T_new = Tm_cond_calc(dt,T_old)
        
        # Step 2. Calculate stagnant lid thickness and Rayleigh number
        from Rayleigh_def import Rayleigh_calc
        Ra[i], d0[i] = Rayleigh_calc(T_old[i_core+1]) #use temp at base of mantle 
        nlid_cells = int(d0[i]/dr)
        lid_start = n_cells - nlid_cells - 1 #index in temp array where lid starts

        if Ra[i] < Rac or cond==1: #once Rayleigh number subcritical don't want to use that criterion anymore
            if cond == 0: #check if first time it is conductive i.e. the switch
                cond_i = i
                cond = 1
            else: 
                pass
            
            Tm_base[i] = T_new[i_core+1]
            Tm_surf[i] = T_new[-2]
            # don't replace any of the mantle temperature, just replace the core
            # Step 4. Replace core with isothermal profile
            dTcdt = dTcdt_calc2(tsolve[i-1],T_old[i_core+1],T_old[i_core],f[i-1]) #save dTcdt seperately as need for f
            Tc[i] = T_old[i_core]+dTcdt*dt
            T_new[:i_core+1] = Tc[i] #replace core including core cell
        else: 
            # Step 3. Replace mantle below stagnant lid with isothermal convective profile
            Tm = T_old[i_core+1] + dTmdt_calc(tsolve[i-1],T_old[i_core+1],T_old[i_core])*dt
            T_new[:lid_start] = Tm
            
            Tm_base[i] = T_new[i_core+1]
            Tm_surf[i] = T_new[-2]
            
            # Step 4. Replace core with isothermal profile
            dTcdt = dTcdt_calc(tsolve[i],T_old[i_core+1],T_old[i_core],f[i],d0[i])
            Tc[i] = T_old[i_core]+ dTcdt*dt
            T_new[:i_core+1] = Tc[i] #replace core including core cell
        
        # Step 5. Calculate f
        if Tc[i] > Tsolidus:
            f[i]=f0  #core not solidifying - based on 3-9 wt% S and phase diagram in Scheinberg 2016 
        else:  
        
            dfdt = -B*dTcdt/(Tc[i-1]*f[i-1])
            f[i] = f[i-1] + dfdt*dt
        
        T_old = T_new 
        
        #write Tprofile to file if appropriate
        if i%ratio== 0: #if multiple of 10Myr then save the profile - will need to check if this works as tsolve might not be integer multiples of 10 ever
            print('t={:.0f}Myr'.format(tsolve[i]/Myr)) #useful to track progress of simulation
            if i_save >= p: # array is full, pass so don't throw an error
                pass
            else:
                Tprofile[i_save,:] = T_old
                i_save = i_save +1 #increment so saves in next space
        else:
            pass
        
        if f[i]>=1: #stop integration if core is solid, truncate arrays to only return non-zero values
            
            Ra = Ra[:i]
            d0 = d0[:i]
            Tprofile = Tprofile[:i,:]
            Tc = Tc[:i]
            Tm_base = Tm_base[:i]
            Tm_surf = Tm_surf[:i]           
            f=f[:i]
            tsolve = tsolve[:i]
            Tprofile = Tprofile[:i_save]
            break
        else:
            pass
              
        
    return Ra, d0, Tprofile, Tc, Tm_base, Tm_surf, f, tsolve, cond_i