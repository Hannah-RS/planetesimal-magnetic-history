#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for solving thermal evolution for the entire evolution of the body
This version uses one mega script to see if it improves performance - a fun experiment.
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
    import scipy.sparse as sp
    from parameters import Myr, Rac, B, Tsolidus, dr, out_interval, rhom, cpm_p, Vm, As, Ts, h0, Al0, XAl, thalf_al, Acmb, km, gamma, alpha_m, g, r, rc, kappa, rhoc, rc, drho, D, Delta, G, cpc, Lc, de_fe, de_k, thalf_fe, thalf_k, f60, ppm_k
    
    #initialise arrays for output
    p = int((tend-tstart)/(out_interval)) #only output temp profiles every 10 Myr
    m = int((tend-tstart)/dt)
    ratio = int(m/p) #use for calculating when to save temp profiles
    i_save=0
    n_cells = len(T0) #number of cells
    i_core = int(n_cells/2) # index in array of last core cell
    cond = 0 #flag for when first switch to conductive regime
    cond_i = np.nan #set as nan and then reset if switches to conduction
    
    Mc = 4/3*np.pi*rc**3*rhoc # mass of core [kg]
    
    #output variables
    Ra = np.zeros([m])
    d0 = np.zeros([m])
    Tprofile= np.zeros([p,n_cells])
    Tc = np.zeros([m])
    Tm_base = np.zeros([m])
    Tm_surf = np.zeros([m])
    f = np.zeros([m])
    tsolve = np.zeros([m])
    
    #define required functions
    
    ################## F_calc ####################################
    def F_calc(f):
        return (1/5+2/15*f**5-(f**2)/3)*f/(1-f**3)
    
    #####################  Tm_cond_calc ###################################
    def Tm_cond_calc(dt,T):
        """
        Expression for dTm/dt from rearranging 13 in supplementary materials of Bryson (2019)
        Calculates rate of change of conductive profile for the whole body and new temperature profile
        Parameters
        ----------
        dt : float
            time step [s]
        T: array
            temperature profile, first value is r=0, last value is surface

        Returns
        -------
        Tm_new :  array
                new temperature profile

        """

             
        # calculate dTdt
        #import stencil
        r_mat = sp.load_npz('Stencil.npz')

            
        #calculate dTdt for conduction
        dTdt = r_mat.dot(T)
        
        Tnew = dTdt*dt + T
            
        return Tnew
    
   ####################### dTmdt_calc ###################################
   def dTmdt_calc(t,Tm,Tc,Qr=True,Qs=True,Qg=True,Ql=True):
       """
       Expression for dTm/dt from rearranging eqn 9 in supplementary materials of Bryson (2019)
       Parameters
       ----------
       t : float
           time, s
       Tm: float
           base mantle temperature
       Tc : float
           core temperature  
       Qr : boolean, optional
           power contribution from radiogenic heating. The default is True.
       Qs : boolean, optional
           power contribution from secular cooling The default is True.
       Qg : boolean, optional
           power contribution from compositional buoyancy. The default is True.
       Ql : boolean, optional
           power contribution from latent heat. The default is True.

       Returns
       -------
       dTmdt : float
               rate of change of mantle temperature 

       """


       
       #calculate lid thickness
       Ra, d0 = Rayleigh_calc(Tm)
       
       #calculate radiogenic heating 
       h = h0*Al0*XAl*np.exp(-np.log(2)*t/thalf_al) 
       rad = h*rhom*Vm #radiogenic heating contribution
       
       #surface flux
       Fs = km*(Tm-Ts)/d0
       
       #CMB heat flux
       Fcmb = Fs*gamma*(Tc-Tm)
       
       return 1/(rhom*cpm_p*Vm)*(rad-Fs*As+Fcmb*Acmb)
   
   ##################### dTcdt_calc - convective CMB heat flux ################################

   def dTcdt_calc(t,Tm,Tc,f,d0,Qr=True,Qs=True,Qg=True,Ql=True):
       """
       Expression for dTc/dt using Equation 70 from Vol 8, Treatise on Geophysics (Nimmo 2007) and Table 1 from Nimmo, F. (2009). Qnergetics of asteroid dynamos and the role of compositional convection. 
       Can toggle on and off effects of gravitational potential energy release, radiogenic heating, secular cooling, latent heat. 

       dTc/dt=(Qcmb-Qr)/(Qst+Qgt+Qlt) where Qst, Qgt, QLt are the expressions for Qs, Qg, Ql divided by dTc/dt


       Parameters
       ----------
       t : float
           time, s
       Tm: float
           base of mantle temperature
       Tc : float
           core temperature
       f : float
           fractional inner core radius 
       d0: float
           stagnant lid thickness
       Qr : boolean, optional
           power contribution from radiogenic heating. The default is True.
       Qs : boolean, optional
           power contribution from secular cooling The default is True.
       Qg : boolean, optional
           power contribution from compositional buoyancy. The default is True.
       Ql : boolean, optional
           power contribution from latent heat. The default is True.

       Returns
       -------
       dTcdt : float
               rate of change of core temperature (if negative core is cooling)

       """

       #check at least one term containing dTc/dt is non zero
       if Qs == False and Qg == False and Ql == False:
           raise ValueError('At least one source of power production due to core cooling must be non-zero')
       else: 
           den = 0 # create a variable which will be the value of the denominator, will be changed later as the criterion has passed
           

           
       #calculate CMB heat flux
       Fcmb = (Tc-Tm)*gamma*km*(Tm-Ts)/d0
       Qcmb = Acmb*Fcmb
       
       num = Qcmb #create a variable which has the value of the numerator, assign Qcmb


       if Qr == True:
           h_fe = de_fe*f60*np.exp(-t/thalf_fe)/thalf_fe #internal heat generation rate from iron [J /kg /s]
           h_k = de_k*ppm_k*np.exp(-t/thalf_k)/thalf_k #internal heat generation rate from potassium
           
           Qr = Mc*(h_fe+h_k)
           num = num - Qr #add radiogenic contribution
       else: pass
       
       
       if Qs == True:
            Qst = -Mc*cpc*(1+2/5*rc**2/D**2)
            den = den + Qst #add secular heating contribution
       else: pass
           
       if Qg == True:
           Mc = 4/3*np.pi*rc**3*rhoc # mass of core [kg]
           F = F_calc(f)           
           Qgt = -3*np.pi*G*rhoc*Mc*F*drho/(Delta-1)*D**2/Tc
           den = den + Qgt #add compositional buoyancy contribution
       else: pass   
       
       if Ql == True:
           Qlt = -3/2*Mc*f*Lc*D**2/(rc**2*Tc*(Delta-1))
           den = den + Qlt #add latent heat contribution
       else: pass
       
       #now have added all terms calculate dTc/dt  
           
       return num/den
    #################  dTcdt_calc2 - conductive CMB heat flux #############################
    
    def dTcdt_calc2(t,Tm,Tc,f,Qr=True,Qs=True,Qg=True,Ql=True):
        """
        Expression for dTc/dt using Equation 70 from Vol 8, Treatise on Geophysics (Nimmo 2007) and Table 1 from Nimmo, F. (2009). Energetics of asteroid dynamos and the role of compositional convection. 
        Can toggle on and off effects of gravitational potential energy release, radiogenic heating, secular cooling, latent heat. 

        dTc/dt=(Qcmb-Qr)/(Qst+Qgt+Qlt) where Qst, Qgt, QLt are the expressions for Qs, Qg, Ql divided by dTc/dt

        This version is for the conductive step, so takes slightly different arguments to the one in solve ivp
        Fcmb is also k(Tc-Tm)/dr now - eqn 14 in Supplementary materials of Bryson (2019)
        Parameters
        ----------
        t: float
            time
        Tm: float
            temperature at the base of the mantle [K]
        Tc : float
            core temperature
        f : float
            fractional inner core radius 
        Qr : boolean, optional
            power contribution from radiogenic heating. The default is True.
        Qs : boolean, optional
            power contribution from secular cooling The default is True.
        Qg : boolean, optional
            power contribution from compositional buoyancy. The default is True.
        Ql : boolean, optional
            power contribution from latent heat. The default is True.

        Returns
        -------
        dTcdt : float
                rate of change of core temperature (if negative core is cooling)

        """

        #check at least one term containing dTc/dt is non zero
        if Qs == False and Qg == False and Ql == False:
            raise ValueError('At least one source of power production due to core cooling must be non-zero')
        else: 
            den = 0 # create a variable which will be the value of the denominator, will be changed later as the criterion has passed
           

            
        #calculate CMB heat flux - eqn 26 in Dodds (2021) or 14 in bryson supplementary materials as flux is equal

        Fcmb = km*(Tc-Tm)/dr
        Qcmb = Acmb*Fcmb
        
        num = Qcmb #create a variable which has the value of the numerator, assign Qcmb


        if Qr == True:
            h_fe = de_fe*f60*np.exp(-t/thalf_fe)/thalf_fe #internal heat generation rate from iron [J /kg /s]
            h_k = de_k*ppm_k*np.exp(-t/thalf_k)/thalf_k #internal heat generation rate from potassium
            
            Qr = Mc*(h_fe+h_k)
            num = num - Qr #add radiogenic contribution
        else: pass
        
        
        if Qs == True:
            Qst = -Mc*cpc*(1+2/5*rc**2/D**2)
            den = den + Qst #add secular heating contribution
        else: pass
            
        if Qg == True:
            
            F = F_calc(f)           
            Qgt = -3*np.pi*G*rhoc*Mc*F*drho/(Delta-1)*D**2/Tc
            den = den + Qgt #add compositional buoyancy contribution
        else: pass   
        
        if Ql == True:
            Qlt = -3/2*Mc*f*Lc*D**2/(rc**2*Tc*(Delta-1))
            den = den + Qlt #add latent heat contribution
        else: pass
        
        #now have added all terms calculate dTc/dt  
            
        return num/den
    
    ##################### Rayleigh_calc ##################################
    
   
    def Rayleigh_calc(Tm):
        """
        Script for calculating the Rayleigh number, expression for d0 comes from substituting in the expression for Rayleigh number
        so that length of the domain cancels

        Parameters
        ----------
        Tm : float
            mantle temperature

        Returns
        -------
        Rayleigh number, stagnant lid thickness

        """
           

        log10_eta = 64 - Tm/29 - 5*np.tanh((Tm-1625)/15)
        eta= 10**log10_eta
        d0 = (gamma/8)**(4/3)*(Tm-Ts)*((Rac*kappa*eta)/(rhom*g*alpha_m))**(1/3) #upper bl
        drh = d0/(gamma*(Tm-Ts)) #lower bl
        Ram= rhom*g*alpha_m*(Tm-Ts)*(r-rc-(d0-drh))**3/(kappa*eta) 
        
        return Ram, d0
    
    ############### Start integration ################################
    #Step 0. Calculate time
    tsolve[0] = tstart + dt
    
    # Step 1. Calculate conductive profile for whole body
    T_new = Tm_cond_calc(dt,T0)
    
    # Step 2. Calculate stagnant lid thickness and Rayleigh number
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
            print('t={:.1f}Myr'.format(tsolve[i]/Myr)) #useful to track progress of simulation
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