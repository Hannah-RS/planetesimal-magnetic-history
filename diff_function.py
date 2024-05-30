#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from stencil import cond_stencil_general
from rayleigh_def import rayleigh_differentiate
from heating import alfe_heating
from dTmdt_def import dTadt_calc
from scipy import sparse as sp
from cp_func import cp_calc_arr, cp_calc_int, cp_calc_eut_arr, cp_calc_eut_int
import numpy as np
from parameters import  ka, rhoa, XFe_a, Xs_0, Xs_eutectic, cpa, Lc, Ts_fe, Tl_fe, Tml, Tms, Ts, As, V, Rac, rcmf, n_cells
def differentiation(Tint,tacc,r,dr,dt):
    """
    Thermal evolution of planetesimal prior to differentiation as described in
    Section 2.4 of Sanderson et. al. (2024)

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
    Xfe = np.zeros([1,n_cells]) #fraction of iron melted
    Xsi = np.zeros([1,n_cells]) #fraction of silicate melted
    cp = np.zeros([1,n_cells]) #specific heat capacity of each cell
    T = np.zeros([1,n_cells]) #temperature
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
        H = np.array([alfe_heating(t[0])])
        Tk = ka*Tint
        #Calculate rhs 1/r^2dt/dr(r^2dt/dr)
        rhs = sparse_mat.dot(Tk) + H*rhoa
        cp[0,:] = cp_calc_arr(Tint,True) #calculate cp
        
        #calculate temperature change
        dTdt = rhs/(rhoa*cp[:,0])
        dTdt_old = dTdt[0] #take temp change at centre for Rayeligh-Roberts number
        T[0,:-1] = Tint[:-1] + dt*dTdt[:-1]
        T[0,-1] = Tint[-1] #pin top cell to 200
        Flid_old = -ka*(Ts - T[0,-2])/dr
        
        Ur = rhoa*H*V/abs(Flid_old*As) #calculate Urey ratio
        
        #calculate melting
        #iron
        Xfe[0,T[0,:]<Ts_fe] = 0 #subsolidus
        Xfe[0,((T[0,:]>=Ts_fe) & (T[0,:]<Tl_fe))] = (T[0,((T[0,:]>=Ts_fe) & (T[0,:]<Tl_fe))]-Ts_fe)/dTphase_fe #melting
        Xfe[0,T[0,:]>=Tl_fe] = 1 #above liquidus
        
        #silicate
        Xsi[0,T[0,:]<Tms] = 0 #subsolidus
        Xsi[0,((T[0,:]>=Tms) & (T[0,:]<Tml))] = (T[0,((T[0,:]>=Tms) & (T[0,:]<Tml))]-Tms)/dTphase_si #melting
        Xsi[0,T[0,:]>=Tml] = 1 #above liquidus
        
        #now loop
        i = 1
        
        while Xsi[i-1,int(n_cells/2)]<rcmf: #assume differentiation occurs at rcmf
            
            #add to existing arrays ready for new timestep
            app_array = np.zeros([1,n_cells])
            T = np.append(T,app_array,0)
            Xfe = np.append(Xfe, app_array, 0)
            Xsi = np.append(Xsi, app_array, 0)
            cp = np.append(cp, app_array, 0)
            Ra = np.append(Ra, 0)
            d0 = np.append(d0, 0)
            Ra_crit = np.append(Ra_crit, 0)
            convect = np.append(convect, 0)
            
            t = np.append(t,t[i-1]+dt)
            
            Ra[i], d0[i], Ra_crit[i], convect[i] = rayleigh_differentiate(t[i],T[i-1,0], Ur)
            #calculate radiogenic heating
            H = np.append(H,alfe_heating(t[i]))
            Tk = ka*T[i-1,:]
            #Calculate rhs 1/r^2dt/dr(r^2dt/dr)
            rhs = sparse_mat.dot(Tk) + H[i]*rhoa
            cp[i,:] = cp_calc_arr(T[i-1,:],True) #calculate cp

            #calculate temperature change
            dTdt = rhs/(rhoa*cp[i,:])
            dTdt_new = dTdt[0]
            T[i,:-1] = T[i-1,:-1] + dt*dTdt[:-1]
            T[i,-1] = Ts
            Flid_new = -ka*(Ts-T[i,-2])/dr #default surface flux
            
            if convect[i-1] == True: #overwrite convecting portion                   
                nlid_cells = round(d0[i]/dr)
                if nlid_cells ==0:
                    lid_start = n_cells -2
                else:
                    lid_start = n_cells - nlid_cells - 1 #index in temp array where lid starts
                
                cp[i,:lid_start] = cp_calc_int(T[i-1,0],True)
                cp[i,-1] = cpa
                dTdt_new = dTadt_calc(t[i-1],T[i-1,lid_start-1],d0[i-1],Flid_old,False)
                T[i,:lid_start] = T[i-1,:lid_start] + dTdt_new*dt 
                T[i,-1] = Ts 
                
                if d0[i] < dr:
                    Flid_new = -ka*(Ts-T[i,lid_start])/d0[i] #if less than grid thickness choose d0 so don't overestimate thickness
                else:
                    Flid_new = -ka*(T[i,lid_start+1]-T[i,lid_start])/dr 
                
            #calculate Urey ratio
            Fs = -ka*(Ts-T[i,-2])/dr
            Ur = rhoa*V*H[i]/(abs(Fs*As))
            #relabel for next step
            Flid_old = Flid_new
            dTdt_old = dTdt_new
            
            #calculate melting
            #iron
            Xfe[i,T[i,:]<Ts_fe] = 0 #subsolidus
            Xfe[i,((T[i,:]>=Ts_fe) & (T[i,:]<Tl_fe))] = (T[i,((T[i,:]>=Ts_fe) & (T[i,:]<Tl_fe))]-Ts_fe)/dTphase_fe #melting
            Xfe[i,T[i,:]>=Tl_fe] = 1 #above liquidus
            
            #silicate
            Xsi[i,T[i,:]<Tms] = 0 #subsolidus
            Xsi[i,((T[i,:]>=Tms) & (T[i,:]<Tml))] = (T[i,((T[i,:]>=Tms) & (T[i,:]<Tml))]-Tms)/dTphase_si #melting
            Xsi[i,T[i,:]>=Tml] = 1 #above liquidus

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
        final temperature profile after differentiation [K]
    Tdiff_profile: float
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
    t_diff: float
        array of timesteps during differentiation [s]

    """
    sparse_mat = sp.dia_matrix(cond_stencil_general(r,dr))
    dTphase_si = Tml - Tms

    #Initial step
    # Create arrays - column is one timestep
    Xfe = np.zeros([1,n_cells]) #fraction of iron melted
    Xsi = np.zeros([1,n_cells]) #fraction of silicate melted
    cp = np.zeros([1,n_cells]) #specific heat capacity of each cell
    T = np.zeros([1,n_cells]) #temperature
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
    H = np.array([alfe_heating(t[0])])
    
    #calculate convective profile
    Tk = ka*Tint
    #Calculate rhs 1/r^2dt/dr(r^2dt/dr)
    rhs = sparse_mat.dot(Tk) + H*rhoa
    #calculate temperature change
    cp[0,:] = cp_calc_eut_arr(Tint,True) #calculate c
    dTdt = rhs/(rhoa*cp[:,0])
    dTdt_old = dTdt[0] #take temp change at centre for Rayeligh-Roberts number
    T[0,:-1] = Tint[:-1] + dt*dTdt[:-1] 
    T[0,-1] = Tint[-1] #pin top cell to 200K
    Flid_old = -ka*(Ts - T[0,-2])/dr
    Ur = rhoa*H[0]*V/abs(Flid_old*As) #calculate Urey ratio
    
    #check for melting
    if np.any((np.int_(Tint)>=Ts_fe) & (Xfe[0,:]<1)): #once solidus is crossed initiate melting
        melt = np.where((np.int_(Tint)>=Ts_fe) & (Xfe[0,:]<1))
        Xfe[0,melt] = rhs[melt]/(rhoa*XFe_a*Lc)*dt
        T[0,melt] = Tint[melt] #overwrite melting cells so their temp is constant      

    #calculate silicate melting  
    Xsi[0,T[0,:]<Tms] = 0 #subsolidus
    Xsi[0,((T[0,:]>=Tms) & (T[0,:]<Tml))] = (T[0,((T[0,:]>=Tms) & (T[0,:]<Tml))]-Tms)/dTphase_si #melting
    Xsi[0,T[0,:]>=Tml] = 1 #above liquidus
    
    #now loop
    i = 1
    
    while Xsi[i-1,int(n_cells/2)]<rcmf: #differentiation occurs at rcmf
        
        app_array = np.zeros([1,n_cells])
        T = np.append(T,app_array,0)
        Xfe = np.append(Xfe, app_array, 0)
        Xsi = np.append(Xsi, app_array, 0)
        cp = np.append(cp, app_array, 0)
        Ra = np.append(Ra, 0)
        d0 = np.append(d0, 0)
        Ra_crit = np.append(Ra_crit, 0)
        convect = np.append(convect, 0)
        
        t = np.append(t,t[i-1]+dt)
        
        Ra[i], d0[i], Ra_crit[i], convect[i] = rayleigh_differentiate(t[i],T[i-1,0],Ur)
        #calculate radiogenic heating
        H = np.append(H,alfe_heating(t[i]))
        Tk = ka*T[i-1,:]
        #Calculate rhs 1/r^2dt/dr(r^2dt/dr)
        rhs = sparse_mat.dot(Tk) + H[i]*rhoa
        cp[i,:] = cp_calc_eut_arr(T[i-1,:],True) #calculate cp
                
        #calculate temperature change
        dTdt = rhs/(rhoa*cp[i,:])
        dTdt_new = dTdt[0]
        T[i,:-1] = T[i-1,:-1] + dt*dTdt[:-1]
        T[i,-1] = Ts
        Flid_new = -ka*(Ts-T[i,-2])/dr #default surface flux
        Fs = -ka*(Ts-T[i,-2])/dr #surface flux for Ur ratio
        Xfe[i,:] = Xfe[i-1,:] #continuity of Xfe between timesteps       
        if np.any((T[i-1,:]>=Ts_fe) & (Xfe[i-1,:]<1)): #see if there are any melting cells
            melt = np.where((T[i-1,:]>=Ts_fe) & (Xfe[i-1,:]<1)&(dTdt>=0)) #only melting where heating up
            Xfe[i,melt] = Xfe[i-1,melt]+rhs[melt]/(rhoa*XFe_a*Lc)*dt
            Xfe[i,-1]=0 #surface node is unmelted by default
            dTdt_new = 0 #no temp change
            T[i,melt] = T[i-1,melt] #melting region has constant temperature
        if convect[i-1] == True: #overwrite convecting portion
            nlid_cells = round(d0[i]/dr)
            if nlid_cells ==0:
                lid_start = n_cells -2
            else:
                lid_start = n_cells - nlid_cells - 1 #index in temp array where lid starts
            cp[i,:lid_start] = cp_calc_eut_int(T[i-1,lid_start-1],True)
            if (int(T[i-1,lid_start-1]) >= Ts_fe) & (Xfe[i-1,lid_start-1]<1): #no temp change only melting
                Vocean = 4/3*np.pi*(r-d0[i])**3
                Alid = 4*np.pi*(r-d0[i])**2
                Xfe[i,:lid_start] = Xfe[i-1,:lid_start]+(Vocean*rhoa*H[i]-Flid_old*Alid)/(rhoa*XFe_a*Lc*Vocean)*dt
                Xfe[i,-1]=0 #surface node is unmelted by default
                if d0[i] < dr:
                    Flid_new = -ka*(Ts-T[i,lid_start])/d0[i] #if less than grid thickness choose d0 so don't overestimate thickness
                else:
                    Flid_new = -ka*(T[i,lid_start+1]-T[i,lid_start])/dr   
            else:
                cp[i,:lid_start] = cp_calc_eut_int(T[i-1,0],True)
                cp[i,-1] = cpa
                dTdt_new = dTadt_calc(t[i-1],T[i-1,lid_start-1],d0[i-1],Flid_old,True)
                T[i,:lid_start] = T[i-1,:lid_start] + dTdt_new*dt 
                T[i,-1] = Ts 
                if d0[i] < dr:
                    Flid_new = -ka*(Ts-T[i,lid_start])/d0[i] #if less than grid thickness choose d0 so don't overestimate thickness
                else:
                    Flid_new = -ka*(T[i,lid_start+1]-T[i,lid_start])/dr      
        
        #calculate Urey ratio
        Ur = rhoa*V*H[i]/(abs(Fs*As))
        #relabel for next step
        Flid_old = Flid_new
        dTdt_old = dTdt_new
            
        #calculate silicate melting        
        Xsi[i,T[i,:]<Tms] = 0 #subsolidus
        Xsi[i,((T[i,:]>=Tms) & (T[i,:]<Tml))] = (T[i,((T[i,:]>=Tms) & (T[i,:]<Tml))]-Tms)/dTphase_si #melting
        Xsi[i,T[i,:]>=Tml] = 1 #above liquidus
                
        i = i+1
           
    #relabel for returning
    Tdiff = T
    t_diff = t

    return Tdiff, Xfe, Xsi, cp, Ra, Ra_crit, convect, d0, t_diff, H