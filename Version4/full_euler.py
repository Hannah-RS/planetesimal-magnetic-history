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
from parameters import Myr, Rac, B, dr, out_interval, km, kc, alpha_c, rhoc, eta_c, gc, cpc, Xs_0

#import required functions
from Tm_cond import T_cond_calc
from dTmdt_def import dTmdt_calc
from dTcdt_def import dTcdt_calc #convective CMB heat flux
from dTcdt_def2 import dTcdt_calc2  #conductive CMB heat flux
from Rayleigh_def import Rayleigh_calc
from cmb_bl import delta_l, delta_c
from fe_fes_liquidus import fe_fes_liquidus 

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
        core temperature, measured at inner core boundary (or centre) [K]
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
    bl: array
        thickness of mantle CMB boundary layer
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
    p = int((tend-tstart)/(out_interval)) #only output temp profiles every 10 Myr
    m = int((tend-tstart)/dt)
    ratio = int(m/p) #use for calculating when to save temp profiles
    i_save=0
    n_cells = len(T0) #number of cells
    i_core = round(n_cells/2)-1 # index in array of last core cell (-1 as indexing starts at 0)
    cond = 0 #flag for when first switch to conductive regime
    cond_i = np.nan #set as nan and then reset if switches to conduction
    
    #output variables
    Xs = np.zeros([m]) #core sulfur fraction
    Ra = np.zeros([m])
    d0 = np.zeros([m])
    bl = np.zeros([m])
    Tprofile= np.zeros([p,n_cells])
    Tc = np.zeros([m])
    Tcmb = np.zeros([m])
    Tm_conv = np.zeros([m])
    Tm_mid = np.zeros([m])
    Tm_surf = np.zeros([m])
    f = np.zeros([m])
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
    
    # Step 1. Calculate conductive profile for mantle
    T_new_mantle = T_cond_calc(dt,T0_mantle,sparse_mat_m)
    Fs[0] = -km*(T_new_mantle[-1]-T_new_mantle[-2])/dr
    
    # Step 2. Is the mantle convecting? Calculate stagnant lid thickness, base thickness and Rayleigh number
    Ra[0], d0[0] = Rayleigh_calc(T0_mantle[1]) #use temp at base of mantle 
    nlid_cells = round(d0[0]/dr)
    lid_start = nmantle_cells - nlid_cells - 1 #index in temp array where lid starts
    base = delta_l(T0_mantle[1],T0_mantle[0])
    nbase_cells = round(base/dr)
    
    if Ra[0] < Rac:# not convecting
       Tm_conv[0] = 0 # convective mantle temperature is 0 if mantle not convecting
       pass 
    else: #mantle is convecting replace mantle below stagnant lid with isothermal convective profile
        Tm = T0_mantle[nbase_cells] + dTmdt_calc(tsolve[0],T0_mantle[nbase_cells],T0_mantle[0])*dt
        T_new_mantle[nbase_cells:lid_start] = Tm
        Tm_conv[0] = Tm #temperature of convecting region
        
    #store values
    Tcmb[0] = T_new_mantle[0]
    Tm_mid[0] = T_new_mantle[round(nmantle_cells/2)] # temperature at the mid mantle
    Tm_surf[0] = T_new_mantle[-2] #temperature one cell above surface pinned to 200K
                
    # Step 3. Calculate conductive profile for the core
    T_new_core = T_cond_calc(dt,T0_core,sparse_mat_c)
        
    # Step 4. Is the core convecting? 
    # check if heat flux is super adiabatic 
    Fcmb[0] = -km*(T_new_mantle[1]-T_new_mantle[0])/dr # CMB heat flux eqn 23 in Dodds 2020
    Fad[0] = kc*T0_core[-2]*alpha_c*gc/cpc
    
    if Fcmb[0] > Fad[0]: #super adiabatic, core convects
        bl[0] = delta_c(T0_core[-2],T0_mantle[0]) #second input is CMB temp
        nbl_cells = round(bl[0]/dr)
        bl_start = ncore_cells - nlid_cells - 1 #index in temp array where lid starts
        
        # is the core solidifying?
        Tliquidus = fe_fes_liquidus(Xs_0)
        if T0_core[-2] < Tliquidus: #core solidifies
            dTcdt = dTcdt_calc2(tsolve[0],T0_mantle[1],T0_core[bl_start-1],f0) #save dTcdt seperately as need for f
            dfdt = -B*dTcdt/(T0_core[-2]*f0)
            f[0] = f0 + dfdt*dt
            Xs[0] = (1-f[0])**(-3)*Xs_0 #update sulfur content
        
        else: # core not solidifying
            dTcdt = dTcdt_calc2(tsolve[0],T0_mantle[1],T0_core[bl_start-1],f0, Qg = False, Ql = False) #Qg, Ql are false as no core solidification
            f[0]=f0
            Xs[0]=Xs_0
        
        #find number of cells which are solid
        nic_cells = round(f[0]/dr)
        T_new_core[nic_cells:bl_start] = T0_core[bl_start-1] + dTcdt*dt #replace everything above the solid core
         
        
        
    else: # don't have whole core convection, 
        #check if there is thermal stratification
        #calculate dT/dr
        dTdr = np.gradient(T_new_core,dr)
        if np.all(dTdr)>0: #whole core is thermally stratified
            # is the core solidifying?
            Tliquidus = fe_fes_liquidus(Xs_0)
            if np.any(T0_core) < Tliquidus: #core solidifies
                raise NotImplementedError('Purely conductive core solidification has not been developed.')
            pass
        else: #only part of the core is stably stratified
            lnb = np.max(np.where(dTdr <= 0)) #index for level of neutral buoyancy
            
            # is the core solidifying?
            Tliquidus = fe_fes_liquidus(Xs_0)
            if T0_core[lnb-1] < Tliquidus: #core solidifies
                dTcdt = dTcdt_calc2(tsolve[0],T0_mantle[1],T0_core[lnb-1],f0) #save dTcdt seperately as need for f
                dfdt = -B*dTcdt/(T0_core[lnb-1]*f0)
                f[0] = f0 + dfdt*dt
                Xs[0] = (1-f[0])**(-3)*Xs_0 #update sulfur content
            
            else: # core not solidifying
                dTcdt = dTcdt_calc2(tsolve[0],T0_mantle[1],T0_core[bl_start-1],f0, Qg = False, Ql = False) #Qg, Ql are false as no core solidification
                f[0]=f0
                Xs[0]=Xs_0
            
            #find number of cells which are solid
            nic_cells = round(f[0]/dr)
            T_new_core[nic_cells:lnb] = T0_core[lnb-1] + dTcdt*dt #replace everything above the solid core
            
        
    # Step 5. Set the CMB temperature in the core array to be the same as the mantle, mantle determines core cooling
    T_new_core[-1] = T_new_mantle[0]
    Tc[0] = T_new_core[nic_cells] #temperature of core is always taken at inner core boundary (or centre)
    
    # Step 6. Replace old array with new ready for next step
    T_old_core = T_new_core
    T_old_mantle = T_new_mantle
 
    print('Moving to loop')
    
    for i in range(1,m):

        #Step 0. Calculate time
        tsolve[i] = tsolve[i-1] + dt
        
        # Step 1. Calculate conductive profile for mantle
        T_new_mantle = T_cond_calc(dt,T_old_mantle,sparse_mat_m)
        Fs[i] = -km*(T_new_mantle[-1]-T_new_mantle[-2])/dr
        
        # Step 2. Is the mantle convecting? Calculate stagnant lid thickness, base thickness and Rayleigh number
        Ra[i], d0[i] = Rayleigh_calc(T_old_mantle[1]) #use temp at base of mantle 
        nlid_cells = round(d0[i]/dr)
        lid_start = nmantle_cells - nlid_cells - 1 #index in temp array where lid starts
        base = delta_l(T_old_mantle[1],T_old_mantle[0])
        nbase_cells = round(base/dr)
        
        if Ra[i] < Rac or cond==1: #once Rayleigh number subcritical don't want to use that criterion anymore
            if cond == 0: #check if first time it is conductive i.e. the switch
                cond_i = i
                cond = 1
            else: #mantle already conducting
                pass
            
            Tm_conv[i] = 0 # convective mantle temperature is 0 if mantle not convecting

        else: #mantle is convecting replace mantle below stagnant lid with isothermal convective profile
            Tm = T_old_mantle[nbase_cells] + dTmdt_calc(tsolve[i],T_old_mantle[nbase_cells],T_old_mantle[0])*dt
            T_new_mantle[nbase_cells:lid_start] = Tm
            Tm_conv[i] = Tm #temperature of convecting region
            
        #store values
        Tcmb[i] = T_new_mantle[0]
        Tm_mid[i] = T_new_mantle[round(nmantle_cells/2)] # temperature at the mid mantle
        Tm_surf[i] = T_new_mantle[-2] #temperature one cell above surface pinned to 200K
                    
        # Step 3. Calculate conductive profile for the core
        T_new_core = T_cond_calc(dt,T_old_core,sparse_mat_c)
            
        # Step 4. Is the core convecting? 
        # check if heat flux is super adiabatic 
        Fcmb[i] = -km*(T_new_mantle[1]-T_new_mantle[0])/dr # CMB heat flux eqn 23 in Dodds 2020
        Fad[i] = kc*T_old_core[-2]*alpha_c*gc/cpc
        
        if Fcmb[i] > Fad[i]: #super adiabatic, core convects
            bl[i] = delta_c(T_old_core[-2],T_old_mantle[0]) #second input is CMB temp
            nbl_cells = round(bl[i]/dr)
            bl_start = ncore_cells - nlid_cells - 1 #index in temp array where lid starts
            
            # is the core solidifying?
            Tliquidus = fe_fes_liquidus(Xs[i-1])
            if T_old_core[-2] < Tliquidus: #core solidifies
                dTcdt = dTcdt_calc2(tsolve[i],T_old_mantle[1],T_old_core[bl_start-1],f[i-1]) #save dTcdt seperately as need for f
                dfdt = -B*dTcdt/(T_old_core[-2]*f[i-1])
                f[i] = f[i-1] + dfdt*dt
                Xs[i] = (1-f[i])**(-3)*Xs_0 #update sulfur content
            
            else: # core not solidifying
                dTcdt = dTcdt_calc2(tsolve[i],T_old_mantle[1],T_old_core[bl_start-1],f0, Qg = False, Ql = False) #Qg, Ql are false as no core solidification
                f[i]=f0
                Xs[i]=Xs_0
            
            #find number of cells which are solid
            nic_cells = round(f[i]/dr)
            T_new_core[nic_cells:bl_start] = T_old_core[bl_start-1] + dTcdt*dt #replace everything above the solid core
             
            
            
        else: # don't have whole core convection, 
            #check if there is thermal stratification
            #calculate dT/dr
            dTdr = np.gradient(T_new_core,dr)
            if np.all(dTdr)>0: #whole core is thermally stratified
                # is the core solidifying?
                Tliquidus = fe_fes_liquidus(Xs[i])
                if np.any(T_old_core) < Tliquidus: #core solidifies
                    raise NotImplementedError('Purely conductive core solidification has not been developed.')
                pass
            else: #only part of the core is stably stratified
                lnb = np.max(np.where(dTdr <= 0)) #index for level of neutral buoyancy
                
                # is the core solidifying?
                Tliquidus = fe_fes_liquidus(Xs[i-1])
                if T_old_core[lnb-1] < Tliquidus: #core solidifies
                    dTcdt = dTcdt_calc2(tsolve[i],T_old_mantle[1],T_old_core[lnb-1],f[i-1]) #save dTcdt seperately as need for f
                    print(dTcdt)
                    dfdt = -B*dTcdt/(T_old_core[lnb-1]*f[i-1])
                    f[i] = f[i-1] + dfdt*dt
                    Xs[i] = (1-f[i])**(-3)*Xs_0 #update sulfur content
                
                else: # core not solidifying
                    dTcdt = dTcdt_calc2(tsolve[i],T_old_mantle[1],T_old_core[bl_start-1],f0, Qg = False, Ql = False) #Qg, Ql are false as no core solidification
                    f[i]=f0
                    Xs[i]=Xs_0
                
                #find number of cells which are solid
                nic_cells = round(f[i]/dr)
                T_new_core[nic_cells:lnb] = T_old_core[lnb-1] + dTcdt*dt #replace everything above the solid core
                
            
        # Step 5. Set the CMB temperature in the core array to be the same as the mantle, mantle determines core cooling
        T_new_core[-1] = T_new_mantle[0]
        Tc[i] = T_new_core[nic_cells] #temperature of core is always taken at inner core boundary (or centre)
        
        # Step 6. Replace old array with new ready for next step
        T_old_core = T_new_core
        T_old_mantle = T_new_mantle

        #write Tprofile to file if appropriate
        if i%ratio== 0: #if multiple of 10Myr then save the profile - will need to check if this works as tsolve might not be integer multiples of 10 ever
            print('t={:.1f}Myr'.format(tsolve[i]/Myr)) #useful to track progress of simulation
            if i_save >= p: # array is full, pass so don't throw an error
                pass
            else:
                T_old = np.hstack([T_old_core,T_old_mantle[1:]])
                Tprofile[i_save,:] = T_old
                i_save = i_save + 1 #increment so saves in next space
        else: 
            pass
        print(f[i])
        if f[i]>=1: #stop integration if core is solid, truncate arrays to only return non-zero values
            
            
            Tc = Tc[:i+1]
            Tcmb = Tcmb[:i+1]
            Tm_mid = Tm_mid[:i+1]
            Tm_conv = Tm_conv[:i+1]
            Tm_surf = Tm_surf[:i+1]
            Tprofile = Tprofile[:i_save+1]
            f=f[:i+1]
            Xs = Xs[:i+1]
            bl = bl[:i+1]
            d0 = d0[:i+1]
            Ra = Ra[:i+1] 
            Fs = Fs[:i+1]
            Fad = Fad[:i+1]
            Fcmb = Fcmb[:i+1]
            tsolve = tsolve[:i+1]
            
            break
        else: 
            pass
              
        
    return Tc, Tcmb, Tm_mid, Tm_conv, Tm_surf, Tprofile, f, Xs, bl, d0, Ra,  Fs, Fad, Fcmb, tsolve, cond_i