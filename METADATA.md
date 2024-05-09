This file describes all the model outputs.

# Npz files
Each of the individual arrays in the npz file can be accessed using the following code where `var` is the name of the array you want to access.
```
npzfile = np.load('run_n.npz')
var = npzfile['var']
```
Units are given in [] brackets.

## run_<n>_diff.npz
The following are all time series arrays saved at frequency `save_interval_d`.

Tdiff: array
        temperature profiles at each step in differentiation [K]
Xfe: array
    proportion of iron in each cell which is melted in differentiation [0 to 1]
Xsi: array
    proportion of silicate in each cell which is melted in differentiation [0 to 1]
cp : array
    effective specific heat capacity of each cell [J /kg /K]
Ra: array
    Rayleigh number for body
Ra_crit: array
    critical Rayleigh number for body
convect : array
    whether the body is convecting
d0 : array
    stagnant lid thickness [m]
t_diff: array
    array of timesteps during differentiation [s]
H : array
    radiogenic heating [W/kg]
        
## run_<n>.npz

The following are all time series arrays saved at frequency `save_interval_t`. 

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
    full temperature profile [K]
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
Racrit: array
    critical Rayleigh number for convecting mantle
Flux : array
    Fs Flux[0,:]: surface heat flux [W m^-2]
    Fad Flux[1,:]: adiabatic CMB heat flux [W m^-2]
    Fcmb Flux[2,:]: CMB heat flux [W m^-2]
    Frad Flux[3,:]: mantle radiogenic heat flux [W m^-2]
    Flid Flux[4,:]: heat flux across base of stagnant lid [W m^-2]
Rem : array
    magnetic Reynolds number
B : array
    dipole magnetic field strength at the surface [T]
buoyr : array
        compositional buoyr[0,:] and thermal buoyr[1,:] buoyancy fluxes [kg/s] 
tsolve: array
    time points corresponding to each of the values above [s]

# Csv files
These files can be loaded using numpy or pandas. The function `load_run_results` from the `load_info.py` module automates loading in results.  

## run_info.csv or auto_params.csv
run : run number
r : planetesimal radius [m]
rcr : fractional core radius 
default : chosen viscosity model - always set to `'vary'`
rcmf : critical melt fraction
eta0 : reference viscosity [Pas]
beta : Arrhenius slope [$K^{-1}$]
w : width of viscosity approximation region [K] - default is 5K
etal : liquid viscosity [Pas]
alpha_n : melt weakening exponent
Xs_0 : initial core sulfur content [wt %] - minimum possible value can be calculated in `Analysis/minimum_xs.py`
Fe0 : $^{60}Fe/^{56}Fe$ in accreted material
t_acc_m : accretion time [Ma after CAI formation]
t_end_m : model end time [Ma after CAI formation] - set this to be long enough for the core to solidify
dr : grid spacing [m]
dt : time step [core conductive timestep]
icfrac : core solidification endmember ($m_{frac}$)

## run_results.csv
All values in this csv are floats. Units are given in the second row of the csv and must be skipped when reading in the csv.

run : run number
tsolid : time of core solidification [Ma after CAI formation]
int_time : total computation time [s]
diff_time : time of planetesimal differentiation
diff_T : temperature of differentiation [K]
peakT : peak mantle temperature [K]
tmax : time of peak mantle temperature [Ma after CAI formation]
peak_coreT : peak core temperature [K]
tcoremax : time of peak core temperature [Ma after CAI formation] 
tstrat_start : time of formation of core thermal stratification [Ma after CAI formation]
tstrat_remove : time of beginning of erosion of core thermal stratification [Ma after CAI formation]
strat_end : time of complete removal of core thermal stratification [Ma after CAI formation]
fcond_t : end of mantle convection [Ma after CAI formation]
fcond_T : temperature at base of the mantle at the cessation of mantle convection [Ma after CAI formation]
tsolid_start : time of onset of core solidification [Ma after CAI formation]
max_Rtherm : maximum magnetic Reynolds number during thermal dynamo generation
max_Rthermt : time of maximum magnetic Reynolds number during thermal dynamo generation [Ma after CAI formation]
max_Btherm : maximum magnetic field strength during thermal dynamo generation [T]
max_Bthermt : time of maximum magnetic field strength during thermal dynamo generation [Ma after CAI formation]
max_Rcomp : maximum magnetic Reynolds number during compostional dynamo generation
max_Bcomp : maximum magnetic field strength during compostional dynamo generation [T]
Bn1 : number of dynamo generation periods for Rem > 10
magon_1, magoff_1, magon_2, magoff_2, magon_3, magoff_3 : on and off times for the 1st, 2nd and 3rd dynamo generation periods for Rem > 10  [Ma after CAI formation]
Bn2 : number of dynamo generation periods for Rem > 40
magon_4, magoff_4, magon_5, magoff_5 :  on and off times for the 1st and 2nd dynamo generation periods for Rem > 40 [Ma after CAI formation]
Bn3 : number of dynamo generation periods for Rem > 100
magon_6, magoff_6, magon_7, magoff_7 on and off times for the 1stand 2nd dynamo generation periods for Rem > 100 [Ma after CAI formation]
