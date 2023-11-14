#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function for differentiating an asteroid. Based on the process in Dodds et. al. (2021)
"""
from stencil import cond_stencil_general
from Rayleigh_def import Rayleigh_differentiate
from heating import AlFe_heating
from dTmdt_def import dTadt_calc
from scipy import sparse as sp
from cp_func import cp_calc_arr, cp_calc_int, cp_calc_eut_arr, cp_calc_eut_int
import numpy as np
from parameters import  ka, rhoa, XFe_a, Xs_0, Xs_eutectic, cpa, Lc, Ts_fe, Tl_fe, Tml, Tms, Ts, As, V, Rac, rcmf, n_cells
def differentiation(Tint,tacc,r,dr,dt):
    """
    

    Parameters
    ----------
    Tint : float
        array of initial temperatures [K]
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
        temperature profiles at each step in differentiation [K]
    Xfe: float
        proportion of iron in each cell which is melted in differentiation [0 to 1]
    Xsi: float
        proportion of silicate in each cell which is melted in differentiation [0 to 1]
    cp : float
        effective specific heat capacity of each cell [J /kg /K]
    Ra: float
        Rayleigh number for body
    Ra_crit: float
        critical Rayleigh number for body
    convect : bool
        whether the body is convecting
    d0 : float
        stagnant lid thickness [m]
    t_diff: float
        array of timesteps during differentiation [s]
    H : float
        radiogenic heating [W/kg]

    """
    sparse_mat = sp.dia_matrix(cond_stencil_general(r,dr))
    dTphase_fe = Tl_fe - Ts_fe
    dTphase_si = Tml - Tms

    #Initial step
    # Create arrays - column is one timestep
    Xfe = np.zeros([n_cells,1]) #fraction of iron melted
    Xsi = np.zeros([n_cells,1]) #fraction of silicate melted
    cp = np.zeros([n_cells,1]) #specific heat capacity of each cell
    T = np.zeros([n_cells,1]) #temperature
    Ra = np.ones([1]) #Rayleigh number
    d0 = np.ones([1]) #stagnant lid thickness
    Ra_crit = np.ones([1]) # critical Rayleigh number
    convect = np.ones([1]) # is anything convecting
    
    t = np.asarray([tacc])
    
    if Xs_0 != Xs_eutectic:
        #Initial step 
        # don't check for convection as accretes isothermal
        Ra[0] = 0
        d0[0] = r
        Ra_crit[0] = Rac
        convect[0] = False
        
        #calculate radiogenic heating
        H = np.array([AlFe_heating(t[0])])
        Tk = ka*Tint
        #Calculate rhs 1/r^2dt/dr(r^2dt/dr)
        rhs = sparse_mat.dot(Tk) + H*rhoa
        cp[:,0] = cp_calc_arr(Tint,True) #calculate cp
        
        #calculate temperature change
        dTdt = rhs/(rhoa*cp[:,0])
        dTdt_old = dTdt[0] #take temp change at centre for Rayeligh-Roberts number
        T[:-1,0] = Tint[:-1] + dt*dTdt[:-1]
        T[-1,0] = Tint[-1] #pin top cell to 200
        Flid_old = -ka*(Ts - T[-2,0])/dr
        
        Ur = rhoa*H*V/abs(Flid_old*As) #calculate Urey ratio
        
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
        
        while Xsi[int(n_cells/2),i-1]<rcmf: #assume differentiation occurs at rcmf
            
            #add to existing arrays ready for new timestep
            app_array = np.zeros([n_cells,1])
            T = np.append(T,app_array,1)
            Xfe = np.append(Xfe, app_array, 1)
            Xsi = np.append(Xsi, app_array, 1)
            cp = np.append(cp, app_array, 1)
            Ra = np.append(Ra, 0)
            d0 = np.append(d0, 0)
            Ra_crit = np.append(Ra_crit, 0)
            convect = np.append(convect, 0)
            
            t = np.append(t,t[i-1]+dt)
            
            Ra[i], d0[i], Ra_crit[i], convect[i] = Rayleigh_differentiate(t[i],T[0,i-1], Ur)
            #calculate radiogenic heating
            H = np.append(H,AlFe_heating(t[i]))
            Tk = ka*T[:,i-1]
            #Calculate rhs 1/r^2dt/dr(r^2dt/dr)
            rhs = sparse_mat.dot(Tk) + H[i]*rhoa
            cp[:,i] = cp_calc_arr(T[:,i-1],True) #calculate cp

            #calculate temperature change
            dTdt = rhs/(rhoa*cp[:,i])
            dTdt_new = dTdt[0]
            T[:-1,i] = T[:-1,i-1] + dt*dTdt[:-1]
            T[-1,i] = Ts
            Flid_new = -ka*(Ts-T[-2,i])/dr #default surface flux
            
            if convect[i-1] == True: #overwrite convecting portion                   
                nlid_cells = round(d0[i]/dr)
                if nlid_cells ==0:
                    lid_start = n_cells -2
                else:
                    lid_start = n_cells - nlid_cells - 1 #index in temp array where lid starts
                
                cp[:lid_start,i] = cp_calc_int(T[0,i-1],True)
                cp[-1,i] = cpa
                dTdt_new = dTadt_calc(t[i-1],T[lid_start-1,i-1],d0[i-1],Flid_old,False)
                T[:lid_start,i] = T[:lid_start,i-1] + dTdt_new*dt 
                T[-1,i] = Ts 
                
                if d0[i] < dr:
                    Flid_new = -ka*(Ts-T[lid_start,i])/d0[i] #if less than grid thickness choose d0 so don't overestimate thickness
                else:
                    Flid_new = -ka*(T[lid_start+1,i]-T[lid_start,i])/dr 
                
            #calculate Urey ratio
            Fs = -ka*(Ts-T[-2,i])/dr
            Ur = rhoa*V*H[i]/(abs(Fs*As))
            #relabel for next step
            Flid_old = Flid_new
            dTdt_old = dTdt_new
            
            #calculate melting
            #iron
            Xfe[T[:,i]<Ts_fe,i] = 0 #subsolidus
            Xfe[((T[:,i]>=Ts_fe) & (T[:,i]<Tl_fe)),i] = (T[((T[:,i]>=Ts_fe) & (T[:,i]<Tl_fe)),i]-Ts_fe)/dTphase_fe #melting
            Xfe[T[:,i]>=Tl_fe,i] = 1 #above liquidus
            
            #silicate
            Xsi[T[:,i]<Tms,i] = 0 #subsolidus
            Xsi[((T[:,i]>=Tms) & (T[:,i]<Tml)),i] = (T[((T[:,i]>=Tms) & (T[:,i]<Tml)),i]-Tms)/dTphase_si #melting
            Xsi[T[:,i]>=Tml,i] = 1 #above liquidus
            if np.any(Xsi[:,i]>rcmf):
                print(Xsi[int(n_cells/2),i])
            i = i+1
            
        #relabel for returning
        Tdiff = T
        t_diff = t
    else:
        Tdiff, Xfe, Xsi, cp, Ra, Ra_crit, convect, d0, t_diff, H = differentiation_eutectic(Tint,tacc,r,dr,dt)
              
    return Tdiff, Xfe, Xsi, cp, Ra, Ra_crit, convect, d0, t_diff, H

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
    Ra: float
        Rayleigh number for body
    Ra_crit: float
        critical Rayleigh number for body
    convect : bool
        whether the body is convecting
    d0 : float
        stagnant lid thickness [m]
    t_diff: float
        array of timesteps during differentiation [s]
    H : float
        radiogenic heating [W/kg]
    t_diff: float
        array of timesteps during differentiation [s]

    """
    sparse_mat = sp.dia_matrix(cond_stencil_general(r,dr))
    dTphase_si = Tml - Tms

    #Initial step
    # Create arrays - column is one timestep
    Xfe = np.zeros([n_cells,1]) #fraction of iron melted
    Xsi = np.zeros([n_cells,1]) #fraction of silicate melted
    cp = np.zeros([n_cells,1]) #specific heat capacity of each cell
    T = np.zeros([n_cells,1]) #temperature
    Ra = np.ones([1]) #Rayleigh number
    d0 = np.ones([1]) #stagnant lid thickness
    Ra_crit = np.ones([1]) # critical Rayleigh number
    convect = np.ones([1]) # is anything convecting
    
    t = np.asarray([tacc])
    
    #Initial step 
        
    # don't check for convection as accretes isothermal   
    Ra[0] = 0
    d0[0] = r
    Ra_crit[0] = Rac
    convect[0] = False
    #calculate radiogenic heating
    H = np.array([AlFe_heating(t[0])])
    
    #calculate convective profile
    Tk = ka*Tint
    #Calculate rhs 1/r^2dt/dr(r^2dt/dr)
    rhs = sparse_mat.dot(Tk) + H*rhoa
    #calculate temperature change
    cp[:,0] = cp_calc_eut_arr(Tint,True) #calculate c
    dTdt = rhs/(rhoa*cp[:,0])
    dTdt_old = dTdt[0] #take temp change at centre for Rayeligh-Roberts number
    T[:-1,0] = Tint[:-1] + dt*dTdt[:-1] 
    T[-1,0] = Tint[-1] #pin top cell to 200K
    Flid_old = -ka*(Ts - T[-2,0])/dr
    Ur = rhoa*H[0]*V/abs(Flid_old*As) #calculate Urey ratio
    
    #check for melting
    if np.any((np.int_(Tint)>=Ts_fe) & (Xfe[:,0]<1)): #once solidus is crossed initiate melting
        melt = np.where((np.int_(Tint)>=Ts_fe) & (Xfe[:,0]<1))
        Xfe[melt,0] = rhs[melt]/(rhoa*XFe_a*Lc)*dt
        T[melt,0] = Tint[melt] #overwrite melting cells so their temp is constant      

    #calculate silicate melting  
    Xsi[T[:,0]<Tms,0] = 0 #subsolidus
    Xsi[((T[:,0]>=Tms) & (T[:,0]<Tml)),0] = (T[((T[:,0]>=Tms) & (T[:,0]<Tml)),0]-Tms)/dTphase_si #melting
    Xsi[T[:,0]>=Tml,0] = 1 #above liquidus
    
    #now loop
    i = 1
    
    while Xsi[int(n_cells/2),i-1]<rcmf: #differentiation occurs at rcmf
        
        app_array = np.zeros([n_cells,1])
        T = np.append(T,app_array,1)
        Xfe = np.append(Xfe, app_array, 1) #automatically append previous timestep
        Xsi = np.append(Xsi, app_array, 1)
        cp = np.append(cp, app_array, 1)
        Ra = np.append(Ra, 0)
        d0 = np.append(d0, 0)
        Ra_crit = np.append(Ra_crit, 0)
        convect = np.append(convect, 0)
        
        t = np.append(t,t[i-1]+dt)
        
        Ra[i], d0[i], Ra_crit[i], convect[i] = Rayleigh_differentiate(t[i],T[0,i-1],Ur)
        #calculate radiogenic heating
        H = np.append(H,AlFe_heating(t[i]))
        Tk = ka*T[:,i-1]
        #Calculate rhs 1/r^2dt/dr(r^2dt/dr)
        rhs = sparse_mat.dot(Tk) + H[i]*rhoa
        cp[:,i] = cp_calc_eut_arr(T[:,i-1],True) #calculate cp
                
        #calculate temperature change
        dTdt = rhs/(rhoa*cp[:,i])
        dTdt_new = dTdt[0]
        T[:-1,i] = T[:-1,i-1] + dt*dTdt[:-1]
        T[-1,i] = Ts
        Flid_new = -ka*(Ts-T[-2,i])/dr #default surface flux
        Fs = -ka*(Ts-T[-2,i])/dr #surface flux for Ur ratio
        Xfe[:,i] = Xfe[:,i-1] #continuity of Xfe between timesteps       
        if np.any((T[:,i-1]>=Ts_fe) & (Xfe[:,i-1]<1)): #see if there are any melting cells
            melt = np.where((T[:,i-1]>=Ts_fe) & (Xfe[:,i-1]<1)&(dTdt>=0)) #only melting where heating up
            Xfe[melt,i] = Xfe[melt,i-1]+rhs[melt]/(rhoa*XFe_a*Lc)*dt
            Xfe[-1,i]=0 #surface node is unmelted by default
            dTdt_new = 0 #no temp change
            T[melt,i] = T[melt,i-1] #melting region has constant temperature
        if convect[i-1] == True: #overwrite convecting portion
            nlid_cells = round(d0[i]/dr)
            if nlid_cells ==0:
                lid_start = n_cells -2
            else:
                lid_start = n_cells - nlid_cells - 1 #index in temp array where lid starts
            cp[:lid_start,i] = cp_calc_eut_int(T[lid_start-1,i-1],True)
            if (int(T[lid_start-1,i-1]) >= Ts_fe) & (Xfe[lid_start-1,i-1]<1): #no temp change only melting
                Vocean = 4/3*np.pi*(r-d0[i])**3
                Alid = 4*np.pi*(r-d0[i])**2
                Xfe[:lid_start,i] = Xfe[:lid_start,i-1]+(Vocean*rhoa*H[i]-Flid_old*Alid)/(rhoa*XFe_a*Lc*Vocean)*dt
                Xfe[-1,i]=0 #surface node is unmelted by default
                if d0[i] < dr:
                    Flid_new = -ka*(Ts-T[lid_start,i])/d0[i] #if less than grid thickness choose d0 so don't overestimate thickness
                else:
                    Flid_new = -ka*(T[lid_start+1,i]-T[lid_start,i])/dr   
            else:
                cp[:lid_start,i] = cp_calc_eut_int(T[0,i-1],True)
                cp[-1,i] = cpa
                dTdt_new = dTadt_calc(t[i-1],T[lid_start-1,i-1],d0[i-1],Flid_old,True)
                T[:lid_start,i] = T[:lid_start,i-1] + dTdt_new*dt 
                T[-1,i] = Ts 
                if d0[i] < dr:
                    Flid_new = -ka*(Ts-T[lid_start,i])/d0[i] #if less than grid thickness choose d0 so don't overestimate thickness
                else:
                    Flid_new = -ka*(T[lid_start+1,i]-T[lid_start,i])/dr      
        
        #calculate Urey ratio
        Ur = rhoa*V*H[i]/(abs(Fs*As))
        #relabel for next step
        Flid_old = Flid_new
        dTdt_old = dTdt_new
            
        #calculate silicate melting        
        Xsi[T[:,i]<Tms,i] = 0 #subsolidus
        Xsi[((T[:,i]>=Tms) & (T[:,i]<Tml)),i] = (T[((T[:,i]>=Tms) & (T[:,i]<Tml)),i]-Tms)/dTphase_si #melting
        Xsi[T[:,i]>=Tml,i] = 1 #above liquidus
                
        i = i+1
           
    #relabel for returning
    Tdiff = T
    t_diff = t

    return Tdiff, Xfe, Xsi, cp, Ra, Ra_crit, convect, d0, t_diff, H