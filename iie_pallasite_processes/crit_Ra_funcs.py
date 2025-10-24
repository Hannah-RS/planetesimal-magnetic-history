import numpy as np
nu = 1e-6  # kinematic viscosity Aubert et. al. 2009 [m^2/s]
alphac = 9.2e-5 # thermal expansivity [K^-1] Sanderson et. al. 2025
omega = 1.75e-4 #10hr period [s^-1] Sanderson et. al. 2025
kappac = 5e-6 # thermal diffusivity [m^2/s] Aubert et. al. 2009
kc = 30 # thermal conductivity [Wm^-1K^-1] Dodd et. al. 2021
G = 6.67e-11 # gravitational constant [m^3/kg/s^2]
dr = 500 #grid spacing [m]
#%% Define functions
def Rac_calc(omega, d):
    """
    Calculate critical Rayleigh number for a rotating system
    
    Parameters
    ----------
    omega : float
        Rotational velocity
    d : float
        liquid core thickness [m]
    Returns
    -------
    Rac : float
        critical Rayleigh number
    """
    
    Ek = nu / (omega * d**2) #Ekman number
    Rac = 11*Ek**(-4/3) #4/3 or -4/3
    #Rac = 1e20
    return Rac

def Ra_calc(dTdr,d,gc):
    """
    Calculate Rayleigh number for a fixed flux condition
    Parameters
    ----------
    dTdr : float
        temperature gradient [K/m]
    d : float
        liquid core thickness [m]
    gc : float
        gravitational field strength at CMB [m/s^2]
    Returns
    -------
    Ra : float
        Rayleigh number
    """
    Ra = alphac*gc*abs(dTdr)*d**4/(kappac*nu)
    return Ra

def rarac_check(Bdat,rind,t,temp,rarac,racrit,tsolid_start,rplot):
    """
    Check if the dynamo is on when each position cools through 593K
    and whether the core is solidifying at each time
    Parameters
    ----------
    Bdat : array
        Magnetic field strength from data [muT]
    rind : array
        Indices within the mantle for each meteorite
    t : array
        Model output time [Myr]
    temp : array
        Model mantle temperature profile [K]
    rarac : array
        Ra/Rac for core
    racrit : float
        Critical Ra/Rac value
    tsolid_start : float
        Onset of core solidification [Myr]
    rplot : array
        Radial grid in mantle [km]
    Returns
    -------
    f : bool
        True if the dynamo is on for the correct meteorites
    depth : array
        Depths for each meteorite (upper and lower bounds) [km]
    tdata : array
        Date for each magnetic remanence (upper and lower bounds) [Myr]
    c : list
        List of booleans for each radiometric date,
        True if the core is solidifying at that time
    """
    c = [] #core solidification check
    fcheck = np.zeros([5]) #counter for dynamo on/off
    depth = np.zeros([5,2]) #depths for each meteorite (upper and lower bounds) [km]
    tdata = np.zeros([5,2]) #dates for each magnetic remanence (upper and lower bounds) [Myr]
    for i, rval in enumerate(rind[:,0]):
        frcheck = np.zeros([2]) #counter for either time for remanance working
        if np.any(temp[:,rval]<=593): #check if cools below 593K
            #find shallowest depth for which dynamo is on/off as expected
            pmax = rind[i,0]-rind[i,-1] #max depth range
            p = 0
            while p<pmax:
                if np.any(temp[:,rval-p]<=593): #if this depth cools below 593K
                    tval = np.where(temp[:,rval-p]<=593)[0][0] #find first time cool below 593K
                    if ((rarac[tval] >= racrit) & (Bdat[i] > 0)) | ((rarac[tval] < racrit) & (Bdat[i] <=0)): 
                        #dynamo is on/off and should be
                        frcheck[0] = 1
                        #save depths, times, field strength
                        depth[i,0] = rplot[-1] - rplot[rind[i,0]-p] 
                        tdata[i,0] = t[np.where(temp[:,rind[i,0]-p]<=593)[0][0]]
                        break
                    else: #if dynamo criteria fails, go one position deeper and repeat
                        p += 1 
                else: #below this depth, no longer cools below 593K
                    break          
                
        #do for upper bound
        if frcheck[0] == 1: #check shallow limit works
            #find upper bound on depth that does cool below 593K and for which dynamo behaves correctly
            p = 0
            while p<pmax:
                while (np.any(temp[:,(rind[i,1]+p)]<=593)==False): #first index that is cool enough
                    p += 1
                #check if dynamo is on/off
                tval = np.where(temp[:,(rind[i,1]+p)]<=593)[0][0]
                if ((rarac[tval] >= racrit) & (Bdat[i] > 0)) | ((rarac[tval] < racrit) & (Bdat[i] <=0)): 
                    #dynamo is on/off and should be
                    frcheck[1] = 1
                    depth[i,1] = rplot[-1] - rplot[rind[i,1]+p]
                    tdata[i,1] = t[np.where(temp[:,(rind[i,1]+p)]<=593)[0][0]]
                    break
                else: #go one index higher and repeat
                    p += 1

            
        if np.any(frcheck == 1): #if either time works
            fcheck[i] = 1 
            if (tsolid_start<=tdata[i,1]): #check if core is solidifying
                c.append(True)
            else:
                c.append(False)
        else:
            c.append(False)

    if np.all(fcheck==1): #if all meteorites work
        f = True
    else:
        f = False
        #overwrite any saved data
        depth = np.zeros([5,2])
        tdata = np.zeros([5,2])
        
    return f, depth, tdata, c