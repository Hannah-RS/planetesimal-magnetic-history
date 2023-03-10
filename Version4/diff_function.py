#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function for differentiating an asteroid. Based on the process in Dodds et. al. (2021)
"""
from stencil import cond_stencil_general
from fe_fes_liquidus import fe_fes_liquidus
from Rayleigh_def import Rayleigh_differentiate
from heating import AlFe_heating
from scipy import sparse as sp
from cp_func import cp_calc_arr, cp_calc_int, cp_calc_eut_arr, cp_calc_eut_int
import numpy as np
from parameters import kc, km, ka, rhoa, rhoc, rhom, XFe_a, Xs_0, Xs_eutectic, cpa, Lc, Myr, Ts_fe, Tl_fe, Tml, Tms, Ts
def differentiation(Tint,tacc,r,dr,dt):
    """
    

    Parameters
    ----------
    Tint : float
        array of initial temperatures
    tacc: float
        accretion time after CAIs [s]
    r : float
        asteroid size [m]
    dr : float
        cell spacing [m]
    dt : float
        timesetp [s]

    Returns
    -------
    Tdiff: float
        final temperature profile after differentiation [K]
    Tdiff_profile: float
        temperature profiles at each step in differentiation [K]
    Xfe: float
        proportion of iron in each cell which is melted in differentiation [0 to 1]
    Xsi: float
        proportion of silicate in each cell which is melted in differentiation [0 to 1]
    cp : float
        effective specific heat capacity of each cell
    t_diff: float
        array of timesteps during differentiation [s]

    """
    sparse_mat = sp.dia_matrix(cond_stencil_general(r,dr))
    ncells = int(r/dr)
    dTphase_fe = Tl_fe - Ts_fe
    dTphase_si = Tml - Tms

    #Initial step
    # Create arrays - column is one timestep
    Xfe = np.zeros([ncells,1]) #fraction of iron melted
    Xsi = np.zeros([ncells,1]) #fraction of silicate melted
    cp = np.zeros([ncells,1]) #specific heat capacity of each cell
    T = np.zeros([ncells,1]) #temperature
    Ra = np.ones([1]) #Rayleigh number
    d0 = np.ones([1]) #stagnant lid thickness
    Ra_crit = np.ones([1]) # critical Rayleigh number
    convect = np.ones([1]) # is anything convecting
    
    t = np.asarray([tacc])
    
    if Xs_0 != Xs_eutectic:
        #Initial step 
        # check for convection
        Ra[0], d0[0], Ra_crit[0], convect[0] = Rayleigh_differentiate(t[0],T[0,0])
        
        #calculate radiogenic heating
        H = np.array([AlFe_heating(t[0])])
        Tk = ka*Tint
        #Calculate rhs 1/r^2dt/dr(r^2dt/dr)
        rhs = sparse_mat.dot(Tk) + H*rhoa
        cp[:,0] = cp_calc_arr(Tint,True) #calculate cp
        
        #calculate temperature change
        dTdt = rhs/(rhoa*cp[:,0])
        T[:-1,0] = Tint[:-1] + dt*dTdt[:-1]
        T[-1,0] = Tint[-1] #pin top cell to 200
        
        if convect[0] == True: #replace convecting portion with convecting profile
            
            nlid_cells = round(d0[0]/dr)
            if nlid_cells ==0:
                lid_start = ncells -2
            else:
                lid_start = ncells - nlid_cells - 1 #index in temp array where lid starts
            cp[:lid_start,0] = cp_calc_int(Tint[0],True) #replace heat capacity below lid
            cp[-1,0] = cpa            
            Fs = -ka*(Ts-Tint[lid_start])/d0[0]
            
            dTdt = (rhoa*H-Fs)/(rhoa*cp[0,0])
            T[:lid_start,0] = Tint[:lid_start] + dTdt*dt
            T[-1,0] = Ts
        
        #calculate melting
        #iron
        Xfe[T[:,0]<Ts_fe,0] = 0 #subsolidus
        Xfe[((T[:,0]>=Ts_fe) & (T[:,0]<Tl_fe)),0] = (T[((T[:,0]>=Ts_fe) & (T[:,0]<Tl_fe)),0]-Ts_fe)/dTphase_fe #melting
        Xfe[T[:,0]>=Tl_fe,0] = 1 #above liquidus
        
        #silicate
        Xsi[T[:,0]<Tms,0] = 0 #subsolidus
        Xsi[((T[:,0]>=Tms) & (T[:,0]<Tml)),0] = (T[((T[:,0]>=Tms) & (T[:,0]<Tml)),0]-Tms)/dTphase_si #melting
        Xsi[T[:,0]>=Tml,0] = 1 #above liquidus
        
        #now loop
        i = 1
        cond_i = 0 #convective switch
        #while (convect[i-1]==False) or (Xsi[int(ncells/2),i-1]<0.05): #whilst any part except the top cell is not convecting and less than 5% silicate melted
        #while t[i-1]<1.16*Myr: 
        #while convect[i-1]==False:
        #while Xsi[int(ncells/2),i-1]<0.2625: #Kathryn criteria for case 1
        while Xsi[int(ncells/2),i-1]<0.05:
            
            app_array = np.zeros([ncells,1])
            T = np.append(T,app_array,1)
            Xfe = np.append(Xfe, app_array, 1)
            Xsi = np.append(Xsi, app_array, 1)
            cp = np.append(cp, app_array, 1)
            Ra = np.append(Ra, 0)
            d0 = np.append(d0, 0)
            Ra_crit = np.append(Ra_crit, 0)
            convect = np.append(convect, 0)
            
            if cond_i == 0:
                t = np.append(t,t[i-1]+dt)
            else: 
                t = np.append(t,t[i-1]+0.01*dt) #adaptive timestep smaller when convecting
                
            Ra[i], d0[i], Ra_crit[i], convect[i] = Rayleigh_differentiate(t[i],T[0,i-1])
            
            #calculate radiogenic heating
            H = np.append(H,AlFe_heating(t[i]))
            Tk = ka*T[:,i-1]
            #Calculate rhs 1/r^2dt/dr(r^2dt/dr)
            rhs = sparse_mat.dot(Tk) + H[i]*rhoa
            cp[:,i] = cp_calc_arr(T[:,i-1],True) #calculate cp

            #calculate temperature change
            dTdt = rhs/(rhoa*cp[:,i])
            T[:-1,i] = T[:-1,i-1] + dt*dTdt[:-1]
            T[-1,i] = Ts
            
            if convect[i] == True or cond_i==1: #overwrite convecting portion
                if cond_i == 0:
                    cond_i =1 
                    print('Convecting, changing timestep')
                    
                nlid_cells = round(d0[i]/dr)
                if nlid_cells ==0:
                    lid_start = ncells -2
                else:
                    lid_start = ncells - nlid_cells - 1 #index in temp array where lid starts
                
                cp[:lid_start,i] = cp_calc_int(T[0,i-1],True)
                cp[-1,i] = cpa
                Fs = -ka*(Ts-T[lid_start,i-1])/d0[i]
                dTdt = (rhoa*H[i]-Fs)/(rhoa*cp[0,i])
                T[:lid_start,i] = T[:lid_start,i-1] + dTdt*0.01*dt 
                T[-1,i] = Ts   
                
            #calculate melting
            #iron
            Xfe[T[:,i]<Ts_fe,i] = 0 #subsolidus
            Xfe[((T[:,i]>=Ts_fe) & (T[:,i]<Tl_fe)),i] = (T[((T[:,i]>=Ts_fe) & (T[:,i]<Tl_fe)),i]-Ts_fe)/dTphase_fe #melting
            Xfe[T[:,i]>=Tl_fe,i] = 1 #above liquidus
            
            #silicate
            Xsi[T[:,i]<Tms,i] = 0 #subsolidus
            Xsi[((T[:,i]>=Tms) & (T[:,i]<Tml)),i] = (T[((T[:,i]>=Tms) & (T[:,i]<Tml)),i]-Tms)/dTphase_si #melting
            Xsi[T[:,i]>=Tml,i] = 1 #above liquidus
            
            i = i+1
            
            #relabel for returning
            Tdiff = T
            t_diff = t
    else:
        Tdiff, Xfe, Xsi, cp, Ra, Ra_crit, convect, t_diff, H = differentiation_eutectic(Tint,tacc,r,dr,dt)
              
    return Tdiff, Xfe, Xsi, cp, Ra, Ra_crit, convect, t_diff, H

def differentiation_eutectic(Tint,tacc,r,dr,dt):
    """
    

    Parameters
    ----------
    Tint : float
        array of initial temperatures
    tacc: float
        accretion time after CAIs [s]
    r : float
        asteroid size [m]
    dr : float
        cell spacing [m]
    dt : float
        timesetp [s]

    Returns
    -------
    Tdiff: float
        final temperature profile after differentiation [K]
    Tdiff_profile: float
        temperature profiles at each step in differentiation [K]
    Xfe: float
        proportion of iron in each cell which is melted in differentiation [0 to 1]
    Xsi: float
        proportion of silicate in each cell which is melted in differentiation [0 to 1]
    cp : float
        effective specific heat capacity of each cell
    t_diff: float
        array of timesteps during differentiation [s]

    """
    sparse_mat = sp.dia_matrix(cond_stencil_general(r,dr))
    ncells = int(r/dr)
    nmid = int(ncells/2)
    dTphase_fe = Tl_fe - Ts_fe
    dTphase_si = Tml - Tms

    #Initial step
    # Create arrays - column is one timestep
    Xfe = np.zeros([ncells,1]) #fraction of iron melted
    Xsi = np.zeros([ncells,1]) #fraction of silicate melted
    cp = np.zeros([ncells,1]) #specific heat capacity of each cell
    T = np.zeros([ncells,1]) #temperature
    Ra = np.ones([1]) #Rayleigh number
    d0 = np.ones([1]) #stagnant lid thickness
    Ra_crit = np.ones([1]) # critical Rayleigh number
    convect = np.ones([1]) # is anything convecting
    
    t = np.asarray([tacc])
    
    #Initial step 
    
    
    # check for convection
    Ra[0], d0[0], Ra_crit[0], convect[0] = Rayleigh_differentiate(t[0],T[0,0])
    
    #calculate radiogenic heating
    H = np.array([AlFe_heating(t[0])])
    
    #calculate convective profile
    Tk = ka*Tint
    #Calculate rhs 1/r^2dt/dr(r^2dt/dr)
    rhs = sparse_mat.dot(Tk) + H*rhoa
    #calculate temperature change
    cp[:,0] = cp_calc_eut_arr(Tint,True) #calculate c
    dTdt = rhs/(rhoa*cp[:,0])
    T[:-1,0] = Tint[:-1] + dt*dTdt[:-1] 
    T[-1,0] = Tint[-1] #pin top cell to 200K
    
    #check for melting
    if np.any((np.int_(Tint)>=Ts_fe) & (Xfe[:,0]<1)): #once solidus is crossed initiate melting
        melt = np.where((np.int_(Tint)>=Ts_fe) & (Xfe[:,0]<1))
        Xfe[melt,0] = rhs[melt]/(rhoa*XFe_a*Lc)*dt
        T[melt,0] = Tint[melt] #overwrite melting cells so their temp is constant    
   
    
    if convect[0] == True: #overwrite conductive profile
        nlid_cells = round(d0[0]/dr)
        if nlid_cells ==0:
            lid_start = ncells -2
        else:
            lid_start = ncells - nlid_cells - 1 #index in temp array where lid starts
                
        Fs = -ka*(Ts-Tint[lid_start])/d0[0]
        if int(Tint[lid_start-1]) >= Ts_fe and (Xfe[lid_start-1,0]<1): #no temp change only melting
            Xfe[:lid_start,0] = (rhoa*H-Fs)/(rhoa*XFe_a*Lc)*dt
        else:
            cp[:lid_start,0] = cp_calc_eut_int(Tint[0],True)
            cp[lid_start:,0] = cp_calc_eut_arr(Tint[lid_start:],True)
            dTdt = (rhoa*H-Fs)/(rhoa*cp[0,0])
            T[:lid_start,0] = Tint[:-1] + dTdt*dt
            T[-1,0] = Ts  

    #calculate silicate melting  
    Xsi[T[:,0]<Tms,0] = 0 #subsolidus
    Xsi[((T[:,0]>=Tms) & (T[:,0]<Tml)),0] = (T[((T[:,0]>=Tms) & (T[:,0]<Tml)),0]-Tms)/dTphase_si #melting
    Xsi[T[:,0]>=Tml,0] = 1 #above liquidus
    
    #now loop
    i = 1
    cond_i = 0 #convective switch
    #while (convect[i-1]==False) or (Xsi[int(ncells/2),i-1]<0.05): #whilst any part except the top cell is not convecting and less than 5% silicate melted
    #while t[i-1]<1.16*Myr: 
    #while convect[i-1]==False:
    #while Xsi[int(ncells/2),i-1]<0.2625: #Kathryn criteria for case 1
    while Xsi[int(ncells/2),i-1]<0.05:
        
        app_array = np.zeros([ncells,1])
        T = np.append(T,app_array,1)
        Xfe = np.append(Xfe, app_array, 1) #automatically append previous timestep
        Xsi = np.append(Xsi, app_array, 1)
        cp = np.append(cp, app_array, 1)
        Ra = np.append(Ra, 0)
        d0 = np.append(d0, 0)
        Ra_crit = np.append(Ra_crit, 0)
        convect = np.append(convect, 0)
        
        if cond_i == 0:
            t = np.append(t,t[i-1]+dt)
        else: 
            t = np.append(t,t[i-1]+0.01*dt) #adaptive timestep smaller when convecting
            
        Ra[i], d0[i], Ra_crit[i], convect[i] = Rayleigh_differentiate(t[i],T[0,i-1])
        
        #calculate radiogenic heating
        H = np.append(H,AlFe_heating(t[i]))
        Tk = ka*T[:,i-1]
        #Calculate rhs 1/r^2dt/dr(r^2dt/dr)
        rhs = sparse_mat.dot(Tk) + H[i]*rhoa
        cp[:,i] = cp_calc_eut_arr(T[:,i-1],True) #calculate cp
                
        #calculate temperature change
        dTdt = rhs/(rhoa*cp[:,i])
        T[:-1,i] = T[:-1,i-1] + dt*dTdt[:-1]
        T[-1,i] = Ts
        Xfe[:,i] = Xfe[:,i-1] #continuity of Xfe between timesteps  
        
        if np.any((np.int_(T)[:,i-1]>=Ts_fe) & (Xfe[:,i-1]<1)): #see if there are any melting cells
            melt = np.where((np.int_(T)[:,i-1]>=Ts_fe) & (Xfe[:,i-1]<1))
            Xfe[melt,i] = Xfe[melt,i-1]+rhs[melt]/(rhoa*XFe_a*Lc)*dt
            T[melt,i] = T[melt,i-1] #melting region has constant temperature
        
        if convect[i] == True or cond_i==1: #overwrite convecting portion
            print('Convecting')
            nlid_cells = round(d0[i]/dr)
            if nlid_cells ==0:
                lid_start = ncells -2
            else:
                lid_start = ncells - nlid_cells - 1 #index in temp array where lid starts
            cp[:lid_start,i] = cp_calc_eut_int(T[lid_start-1,i-1],True)
                        
            if cond_i == 0: #first time it starts convecting
                print('Convecting, changing timestep')             
                cond_i =1
                
            Fs = -ka*(Ts-T[lid_start,i-1])/d0[i]
            if int(T[lid_start-1,i-1]) >= Ts_fe and (Xfe[lid_start-1,i-1]<1): #no temp change only melting
                Xfe[:lid_start,i] = Xfe[:lid_start,i-1]+(rhoa*H-Fs)/(rhoa*XFe_a*Lc)*dt
            else:
                cp[:lid_start,i] = cp_calc_eut_int(T[:lid_start,i-1],True)
                cp[lid_start:,i] = cp_calc_eut_arr(T[lid_start:,i-1],True)
                dTdt = (rhoa*H[i]-Fs)/(rhoa*cp[0,i])
                T[:lid_start,i] = T[:lid_start,i-1] + dTdt*0.01*dt
                T[-1,i] = Ts 
                
     
        #calculate silicate melting        
        Xsi[T[:,i]<Tms,i] = 0 #subsolidus
        Xsi[((T[:,i]>=Tms) & (T[:,i]<Tml)),i] = (T[((T[:,i]>=Tms) & (T[:,i]<Tml)),i]-Tms)/dTphase_si #melting
        Xsi[T[:,i]>=Tml,i] = 1 #above liquidus
        
        i = i+1
        
        #relabel for returning
        Tdiff = T
        t_diff = t

    return Tdiff, Xfe, Xsi, cp, Ra, Ra_crit, convect, t_diff, H