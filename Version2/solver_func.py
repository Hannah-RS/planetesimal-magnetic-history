#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functionalised form of solver.py script for solving the dTcdt ODE defined in dTcdt_def
"""
def solver_func(run,t_end_m,Ek=True,Er=True,Es=True,Eg=True,El=True):
    """
    

    Parameters
    ----------
    run: int
         run number
    t_end_m: float,
        time for end of integration in Myr
    Ek : boolean, optional
        entropy contribution from core conductivity. The default is True. If True the contribution is included
    Er : boolean, optional
        entropy contribution from radiogenic heating. The default is True.
    Es : boolean, optional
        entropy contribution from secular cooling The default is True.
    Eg : boolean, optional
        entropy contribution from compositional buoyancy. The default is True.
    El : boolean, optional
        entropy contribution from latent heat. The default is True.

    Returns
    -------
    None. Saves a file with parameters of the run

    """
    # import modules
    import scipy.integrate as scint
    import numpy as np
    
    # import function to be integrated and event function to stop integration
    from dTcdt_def import dTcdt_calc
    from stop_def import stop
    
    stop.terminal = True # make the event terminal so it will stop the integration
    
    #import time constants and initial conditions
    from parameters import year, Myr, Tc0, f0, B
    
    # define the parameters for the integration and run number
    t_end=t_end_m*Myr
    step_m=0.0001*t_end
    
    # perform integration 
    yout=scint.solve_ivp(dTcdt_calc, [0,t_end], np.array([Tc0]), args= (Ek,Er,Es,Eg,El), max_step=step_m, events=stop)
    
    # assign the output to variables 
    Tc_solve=yout['y'][0,:]
    t_solve=yout['t']
    
    # calculate rate of change of core temperature
    #make sure to match argument here with optional arguments in integration so you include the same terms both times
    dTcdt = dTcdt_calc(t_solve,Tc_solve)
    
    # calculate fractional inner core radius
    from f_def import f_calc
    f_out=f_calc(Tc_solve,Tc0,f0,B)
    
    # calculate rate of change of fractional inner core radius
    from dfdt_def import dfdt_calc
    dfdt=dfdt_calc(f_out,Tc_solve,dTcdt)
    
    # calculate magentic Reynolds number
    from Rem_func import Rem_calc
    Rem = Rem_calc(f_out[:-1], dfdt[:-1]) #truncate array as last value in f is so close to 1 that 1-f produces a runtimewarning
    
    #save variables to file
    np.savez('Results/run_{}'.format(run), Tc=Tc_solve, t=t_solve, f=f_out, dTcdt=dTcdt, dfdt=dfdt, Rem=Rem) 
    
    return None
